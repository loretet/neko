# Verification of the AMD SGS model

This case checks the CPU implementation of the anisotropic minimum-dissipation
(AMD) model (`model: "amd"` in the `les_model` simulation component) against
analytic values, including the anisotropic filter-width tensor
`G` of `les_model_compute_delta_tensor`.

The fluid and the temperature are frozen (`"freeze": true`) and set to
polynomial fields: a velocity with a constant gradient tensor `A`,
`u_i = A_ij x_j`, and `theta = 300 + 6 x - 4 y + 0.05 z + 0.3 z^2`. The spectral
element derivatives of these fields are exact, so the eddy viscosity `nut` and
the eddy diffusivity `temperature_alphat` computed by Neko must agree with the
AMD formulas (Abkar & Moin 2017, Eqs. 17-18), evaluated analytically in the user
file, to round-off. The user file prints the maximum value and the maximum
absolute error of both fields in the log, in lines starting with `AMD TEST`.

Generate the mesh and build the case with

```
genmeshbox 0 2 0 2 0 1 4 4 4 .false. .false. .false.
makeneko amd_verification.f90
```

The velocity field is selected with the environment variable `AMD_TEST_MODE`:

| `AMD_TEST_MODE` | velocity field                  | expected result               |
|-----------------|---------------------------------|-------------------------------|
| `base` (default)| general gradient tensor `A`     | `nut`, `alphat` > 0           |
| `shear`         | laminar shear, `u = 0.8 z`      | `nut = 0`                     |
| `rotation`      | solid-body rotation about `z`   | `nut = 0`                     |
| `negA`          | `-A` (reversed flow)            | `nut = 0`, `alphat = 0` (clipping) |
| `skew`          | `A` on a sheared mesh, `x -> x + 0.3 y + 0.2 z` | tests the off-diagonal filter widths |
| `neutral`       | `A`, run with `amd_verification_neutral.case` | no buoyancy term |
| `nomean`        | `A`, run with `amd_verification_nomean.case` | no horizontal-mean subtraction |

For example

```
AMD_TEST_MODE=shear mpirun -np 2 -x AMD_TEST_MODE ./neko amd_verification.case
```

In all cases the reported errors should be of the order of machine precision.
Remember to rerun `makeneko` after rebuilding Neko, since the case executable
is linked against the Neko library at build time.
