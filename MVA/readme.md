# MVA for magnetic clouds

Minimum variance analysis (MVA) of in-situ magnetic field data, made for magnetic clouds and other flux ropes in CMEs. It finds
- the **flux-rope axis** (intermediate eigenvector), and/or
- the **local normal of the cloud front** (minimum eigenvector),
and estimates their errors with a bootstrap that resamples only the noise.

## Quick start

```python
from mva import mva, plot_bootstrap

# flux-rope axis: the whole cloud interval
eigen, bs = mva(b1, b2, b3, normalization=True, smooth=True, n=2000)
axis = eigen["vect"][:, 1]
print(f"theta = {eigen['theta']['j']:.1f} ± {bs['theta_err']['j']:.1f} deg")
print(f"phi   = {eigen['phi']['j']:.1f} ± {bs['phi_err']['j']:.1f} deg")
print(f"68 % error cone = {bs['cone68']['j']:.1f} deg")
plot_bootstrap(eigen, bs, "j")

# local normal of the front: a short window across the front, raw B
eigen, bs = mva(b1[a:b], b2[a:b], b3[a:b], normalization=False, n=2000)
normal = eigen["vect"][:, 0]
```

The data must be **evenly sampled**. Keep data gaps as `NaN`; they are skipped.

## Parameters

| parameter | default | meaning |
|---|---|---|
| `normalization` | `False` | Divide B by \|B\| before MVA. `True` for the axis, `False` for the front normal. |
| `smooth` | `True` | Only used with `normalization=True`. `True`: divide by a cubic fit of \|B\| (recommended). `False`: divide by \|B\| at each point (not recommended when there is noise). |
| `n` | `2000` | Number of bootstrap copies. `0` skips the bootstrap (`bs` is `None`). More copies give smoother statistics, not smaller errors. |
| `ref` | `(1, 0, 0)` | The minimum eigenvector i is turned so that i · ref > 0. Use the anti-sunward direction of your frame: `(1, 0, 0)` in RTN, `(-1, 0, 0)` in GSE. |
| `seed` | `0` | Random seed of the bootstrap. |

## Output

### `eigen`: the MVA result

| key | meaning |
|---|---|
| `value` | eigenvalues (i, j, k), smallest first |
| `vect` | 3×3 array; columns are the unit vectors i, j, k |
| `ratio_ji`, `ratio_kj` | eigenvalue ratios λj/λi and λk/λj |
| `theta`, `phi` | direction angles (deg) of each vector, keys `'i'`, `'j'`, `'k'` |
| `K` | number of data points used |

The three eigenvectors:

| vector | variance | meaning | sign |
|---|---|---|---|
| i | minimum | whole cloud: roughly the crossing direction. Short window across the front: the local normal. | i · ref > 0 |
| j | intermediate | the flux-rope axis, if the spacecraft passes not too far from the axis | along the field, ⟨B · j⟩ > 0 |
| k | maximum | completes the frame | k = i × j |

The angles are measured in the input frame: θ is the elevation from the
plane of axes 1 and 2, and φ is the azimuth from axis 1 toward axis 2
(0–360°). In RTN, θ is measured from the R–T plane and φ from R toward T. In
GSE they are the usual latitude and longitude.

### `bs`: the errors (1σ, in degrees)

| key | meaning |
|---|---|
| `dphi['ij']`, `dphi['ik']`, `dphi['jk']` | spread of the turn of one eigenvector toward another (see below) |
| `cone68` | 68 % of the bootstrap copies lie within this angle of the MVA result; keys `'i'`, `'j'`, `'k'` |
| `theta_err`, `phi_err` | errors of θ and φ; keys `'i'`, `'j'`, `'k'` |
| `dphi_analytic` | the same as `dphi`, from the formula of Khrabrov & Sonnerup (1998), using `k_eff` |
| `r` | correlation between neighbouring noise values |
| `k_eff` | number of independent points, K(1 − r)/(1 + r) |
| `vect` | (n, 3, 3) eigenvectors of every bootstrap copy |

## Which settings to use

**For the axis, use `normalization=True, smooth=True`.** The field strength
\|B\| changes across a cloud, for example because the cloud expands while the
spacecraft is inside it. With raw B, this \|B\| profile tilts the MVA axis.
Dividing by \|B\| removes it. Dividing by a *cubic fit* of \|B\|, instead of
\|B\| at each point, keeps the noise from shifting the axis: a cubic can follow
the slow \|B\| profile but not the noise.

**For the front normal, use `normalization=False`.** The method relies on the
normal component B · n being constant across the front (from ∇ · B = 0). That
holds for B itself, not for B/\|B\|.

## How the errors are estimated

### The bootstrap

If we had many crossings of the same structure, each with different noise, we
could run MVA on each one and see how much the result changes. We only have one
crossing, so the bootstrap makes fake copies from the data itself:

1. **Split the data into signal and noise.** A smooth polynomial fit of each
   component (degree up to 7) is the signal. What is left over is the noise.
2. **Measure how the noise behaves.** `r` is the correlation between
   neighbouring noise values. High-cadence data often have r ≈ 0.5–0.9.
3. **Make a copy.** Shuffle the noise, rebuild it with the same `r`, and add it
   back to the signal.
4. **Run each copy through exactly the same steps as the real data**: the same
   normalization (the cubic fit of \|B\| is redone for every copy) and MVA.
5. **Repeat `n` times.** The spread of the results is the error.

Only the noise is resampled. The plain bootstrap (Kawano & Higuchi 1995)
resamples whole data points, so it also reshuffles the rope's own field. That
works for the minimum direction, but it makes the axis error toward k two to four
times too large (see `test_bootstrap.py`).

### What `dphi` and "tilt" mean

Each bootstrap copy gives an axis that points in a slightly different
direction. Any small change in the direction of j can be split into two
parts: a turn **toward i** and a turn **toward k**. These are the "tilts".

- For the **axis j**, `dphi['ij']` is the error toward i and `dphi['jk']` is
  the error toward k. Together they form an error ellipse.
- For the **normal i**, `dphi['ij']` is the error toward j and `dphi['ik']`
  is the error toward k.

A turn of j toward i is the same rotation as a turn of i toward j, so one
number, `dphi['ij']`, serves both vectors.

`cone68` sums this up in one number, and `theta_err` / `phi_err` give the same
errors as direction angles: the code computes θ and φ for every copy and takes
their spread.

![Bootstrap example](example_bootstrap.png)

*Left: each dot is one bootstrap copy of the axis, placed by its tilt toward k
and toward i. Blue ellipse: 1σ (`dphi`). Orange circle: `cone68`. Right: the
same copies as θ and φ.*

## Checks before trusting a result

- **Eigenvalue ratios.** `ratio_ji` and `ratio_kj` should be well above 1. If a
  ratio is close to 1, the two directions cannot be told apart. For a front
  normal, `ratio_ji` ≥ 10 is a common rule of thumb (Sonnerup & Scheible 1998).
- **Histograms.** In `plot_bootstrap`, the cloud of copies should be compact and
  roughly elliptical. If it is not, trust the plot, not the numbers.
- **Noise correlation.** If `r` is large (above about 0.5) on a short window,
  average the data first so that neighbouring points are nearly independent.
- **Near the pole.** When \|θ\| is above about 70°, φ has no clear meaning;
  quote `cone68` instead of `phi_err`.
- **Robustness.** Compare with `normalization=False`, and move the start and end
  of the interval a little. If the axis changes a lot, it is not well
  determined.


## References

- Kawano, H., & Higuchi, T. (1995). The bootstrap method in space physics: Error estimation for the minimum variance analysis. GRL, 22, 307.
- Khrabrov, A. V., & Sonnerup, B. U. Ö. (1998). Error estimates for minimum variance analysis. JGR, 103, 6641.
- Sonnerup, B. U. Ö., & Cahill, L. J. (1967). Magnetopause structure and attitude from Explorer 12 observations. JGR, 72, 171.
- Sonnerup, B. U. Ö., & Scheible, M. (1998). Minimum and maximum variance analysis. In Analysis Methods for Multi-Spacecraft Data, ISSI SR-001, 185.
- Xiao, C. J., et al. (2004). Inferring of flux rope orientation with the minimum variance analysis technique. JGR, 109, A11218.
