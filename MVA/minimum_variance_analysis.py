"""
mva.py - minimum variance analysis (MVA) of the magnetic field, with errors.

    from mva import mva
    eigen, bs = mva(b1, b2, b3, normalization=True, smooth=True, n=2000)

Works in any right-handed frame (RTN, GSE, ...): b1, b2, b3 are the field
components in that frame, and all output vectors are in the same frame.
Data must be evenly sampled; keep gaps as NaN (they are skipped).

Eigenvectors (columns of eigen['vect'], eigenvalues ascending):
  i  minimum       short window across the front: local normal
                   -> normalization=False (B.n is constant only for B itself)
  j  intermediate  rope crossed not too far from its axis: the axis
                   -> normalization=True, smooth=True
  k  maximum       k = i x j

Normalization (only if normalization=True):
  smooth=True   divide B by a cubic fit of |B| over the interval: removes the
                large-scale |B| profile (e.g. expansion) but cannot follow the
                noise, so it does not shift the axis.
  smooth=False  divide by |B| at each point. Not recommended with noise: it
                can shift the axis by tens of degrees (test_noise.py).

Errors (bootstrap, n > 0) cover random noise only. Errors from the impact
distance, expansion and rope structure are often larger (test_expansion.py).
"""
import numpy as np
from scipy.signal import lfilter

PAIRS = {"ij": (0, 1), "ik": (0, 2), "jk": (1, 2)}


def _normalize(B, pos, normalization, smooth):
    if not normalization:
        return B
    m = np.linalg.norm(B, axis=1)
    if smooth:
        m = np.polynomial.Polynomial.fit(pos, m, 3)(pos)
    return B / m[:, None]


def _eig(X, ref):
    lam, V = np.linalg.eigh(np.cov(X.T))
    lam = np.maximum(lam, 0.0)                       # round-off can give tiny negatives
    i, j = V[:, 0], V[:, 1]
    if X.mean(axis=0) @ j < 0:                       # axis along the field
        j = -j
    if i @ ref < 0:
        i = -i
    return lam, np.column_stack((i, j, np.cross(i, j)))


def _angles(v):
    """theta: elevation from the 1-2 plane; phi: azimuth from axis 1 toward axis 2 (deg)."""
    return (np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1))),
            np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360)


def mva(b1, b2, b3, normalization=False, smooth=True, n=2000, ref=(1, 0, 0), seed=0):
    """
    normalization : divide B by |B| before MVA (use for the axis)
    smooth        : with normalization, use a cubic fit of |B| (recommended)
    n             : number of bootstrap copies (0: no bootstrap, bs = None).
                    More copies give smoother statistics, not smaller errors.
    ref           : i is turned so that i . ref > 0. Give the anti-sunward
                    direction of your frame: (1,0,0) in RTN, (-1,0,0) in GSE.

    eigen (the MVA result):
      value         eigenvalues (i, j, k), ascending
      vect          3x3, columns i, j, k
      ratio_ji      value_j / value_i,  ratio_kj  value_k / value_j
      theta, phi    direction angles (deg) of 'i', 'j', 'k'

    bs (the errors, 1 sigma, deg):
      dphi          spread of the turn of one eigenvector toward another:
                    'ij' (i <-> j), 'ik' (i <-> k), 'jk' (j <-> k)
                    axis j: dphi['ij'] toward i, dphi['jk'] toward k
                    normal i: dphi['ij'] toward j, dphi['ik'] toward k
      cone68        68 % of the copies lie within this angle: 'i', 'j', 'k'
      theta_err, phi_err   errors of theta and phi: 'i', 'j', 'k'
                    (phi_err is meaningless when |theta| > ~70 deg)
      dphi_analytic Khrabrov & Sonnerup (1998) eq. 24 with k_eff, same keys as dphi
      r, k_eff      noise correlation between neighbouring points, and the
                    number of independent points K (1-r)/(1+r)
      vect          (n, 3, 3) eigenvectors of every copy (same sense as eigen)

    Bootstrap (only the noise is resampled):
      1. Fit each component with a polynomial (degree min(7, K//10)): the
         signal S. The rest, R = B - S, is the noise.
      2. r = lag-1 correlation of R; independent steps e = R[k] - r R[k-1].
      3. Copy = S + noise rebuilt from randomly drawn steps with the same r.
         Each copy goes through the same normalization/smoothing and MVA as
         the real data. The spread of the n results is the error.
    """
    B = np.column_stack((b1, b2, b3)).astype(float)
    pos = np.arange(len(B))
    ok = np.all(np.isfinite(B), axis=1)
    B, pos = B[ok], pos[ok]
    ref = np.asarray(ref, float)

    lam, E = _eig(_normalize(B, pos, normalization, smooth), ref)
    theta, phi = _angles(E.T)
    with np.errstate(divide="ignore"):
        eigen = dict(value=lam, vect=E, K=len(B),
                     ratio_ji=lam[1] / lam[0], ratio_kj=lam[2] / lam[1],
                     theta=dict(zip("ijk", theta)), phi=dict(zip("ijk", phi)))
    if not n:
        return eigen, None

    # 1-2: split the data into signal and noise
    K = len(B)
    S = np.column_stack([np.polynomial.Polynomial.fit(pos, c, min(7, max(1, K // 10)))(pos)
                         for c in B.T])
    R = B - S
    r = float(np.clip(np.sum(R[:-1] * R[1:]) / np.sum(R * R), 0.0, 0.99))
    steps = R[1:] - r * R[:-1]
    steps -= steps.mean(axis=0)

    # 3: copies, each through the same pipeline as the real data
    rng = np.random.default_rng(seed)
    burn = 100                                        # let the rebuilt noise settle
    vect = np.empty((n, 3, 3))
    for s in range(n):
        noise = lfilter([1.0], [1.0, -r], steps[rng.integers(0, K - 1, K + burn)], axis=0)[burn:]
        Es = _eig(_normalize(S + noise, pos, normalization, smooth), ref)[1]
        vect[s] = Es * np.sign(np.sum(Es * E, axis=0))  # same sense as E

    # statistics of the copies
    tilt = np.degrees(np.arcsin(np.clip(np.einsum("sia,ib->sab", vect, E), -1, 1)))
    th, ph = _angles(vect.transpose(0, 2, 1))
    ph = (ph - phi + 180) % 360 - 180                 # azimuth change, no 0/360 jump
    k_eff = K * (1 - r) / (1 + r)
    l3 = lam[0]
    bs = dict(
        dphi={p: float(tilt[:, a, b].std(ddof=1)) for p, (a, b) in PAIRS.items()},
        cone68={c: float(np.percentile(90 - tilt[:, a, a], 68.3)) for a, c in enumerate("ijk")},
        theta_err={c: float(th[:, a].std(ddof=1)) for a, c in enumerate("ijk")},
        phi_err={c: float(ph[:, a].std(ddof=1)) for a, c in enumerate("ijk")},
        dphi_analytic={p: float(np.degrees(np.sqrt(l3 * (lam[a] + lam[b] - l3)
                                                   / ((k_eff - 1) * (lam[a] - lam[b]) ** 2))))
                       for p, (a, b) in PAIRS.items()},
        r=r, k_eff=k_eff, vect=vect)
    return eigen, bs


def plot_bootstrap(eigen, bs, vector="j", fname=None):
    """
    Left : the error cone seen end-on. Each dot is one copy; its position is
           how far that copy's vector turned toward each of the other two
           eigenvectors (the 'tilts'). Ellipse: 1 sigma (dphi). Circle: cone68.
    Right: the same copies as direction angles theta and phi.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Ellipse
    c = "ijk".index(vector)
    o1, o2 = [d for d in range(3) if d != c]
    E, V = eigen["vect"], bs["vect"][:, :, c]
    t1, t2 = (np.degrees(np.arcsin(np.clip(V @ E[:, d], -1, 1))) for d in (o1, o2))
    s1, s2 = (bs["dphi"]["".join(sorted(vector + "ijk"[d]))] for d in (o1, o2))
    th, ph = _angles(V)
    ph0 = eigen["phi"][vector]
    ph = ph0 + (ph - ph0 + 180) % 360 - 180

    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.8))
    a.scatter(t2, t1, s=3, color="#A8BCDB", lw=0)
    a.add_patch(Ellipse((0, 0), 2 * s2, 2 * s1, fill=False, color="#4C72B0", lw=2, label="1 σ (dphi)"))
    a.add_patch(Circle((0, 0), bs["cone68"][vector], fill=False, color="#DD8452", lw=2,
                       ls="--", label=f"cone68 = {bs['cone68'][vector]:.1f}°"))
    a.plot(0, 0, "k+", ms=12, label=f"MVA {vector}")
    lim = 1.1 * max(np.percentile(np.hypot(t1, t2), 99.5), bs["cone68"][vector])
    a.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal",
          xlabel=f"tilt toward {'ijk'[o2]}  (deg)", ylabel=f"tilt toward {'ijk'[o1]}  (deg)",
          title=f"{vector}: error cone seen end-on")
    a.legend(frameon=False, fontsize=9, loc="upper right")

    b.scatter(ph, th, s=3, color="#A8BCDB", lw=0)
    b.errorbar(ph0, eigen["theta"][vector], xerr=bs["phi_err"][vector],
               yerr=bs["theta_err"][vector], fmt="k+", ms=12, capsize=4, lw=1.5)
    b.set(xlabel="φ  (deg)", ylabel="θ  (deg)",
          title=f"θ = {eigen['theta'][vector]:.1f} ± {bs['theta_err'][vector]:.1f}°,   "
                f"φ = {ph0:.1f} ± {bs['phi_err'][vector]:.1f}°")
    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    if fname:
        fig.savefig(fname, dpi=150)
    return fig
