# OpenFOAM Numerics

Custom gradient schemes for OpenFOAM that adapt to mesh quality. On well-shaped cells the schemes favour accurate least-squares or Gauss reconstructions; on skewed, non-orthogonal, or poorly conditioned cells they blend toward a more robust alternative.

The library ships as `libnumerics`, plus a small preprocessing utility that can precompute blending weights once per case.

---

## Features

| Scheme | Role |
| --- | --- |
| `blendedGradScheme` | Blend any two gradient schemes using a precomputed `gradSchemeBlending` field |
| `hybridLeastSquares` | Least-squares gradients weighted by the blending field |
| `hybridLeastSquares1` | Variant of hybrid least-squares with corrected boundary treatment |
| `hybridFaceLeastSquares` | Face-oriented hybrid least-squares using the same blending field |
| `hybridLsqGG` | Runtime hybrid of least-squares and Green–Gauss from cell-quality helpers |
| `hybridGaussLS` | Blend Gauss and least-squares from non-orthogonality / skewness |
| `hybridIterGaussLS` | Smooth blend of iterative Gauss and least-squares |

> Which schemes are compiled depends on the branch (see below).

---

## Branches

| Branch | Target | Schemes built into `libnumerics` |
| --- | --- | --- |
| `ESI-v2412` | OpenFOAM.com v2412 (ESI) | `blendedGradScheme`, `hybridLeastSquares`, `hybridLeastSquares1`, `hybridFaceLeastSquares` |
| `OpenFoam_org` | OpenFOAM.org | `blendedGradScheme`, `hybridLeastSquares` |
| `TCAE` | TCAE / related stacks | `blendedGradScheme` |
| `main` | General / experimental | `blendedGradScheme`, `hybridGaussLS`, `hybridIterGaussLS`, `hybridLsqGG` |

Check out the branch that matches your OpenFOAM installation before compiling.

---

## Requirements

- A working OpenFOAM environment (`wmSet` / `source` your OpenFOAM `bashrc` so `wmake`, `$FOAM_USER_LIBBIN`, and `$FOAM_USER_APPBIN` are set)
- Compatible compilers for that OpenFOAM build

---

## Build

```bash
cd OpenFoamNumerics
./Allwmake
```

This compiles:

- **Library** → `$FOAM_USER_LIBBIN/libnumerics.so`
- **Utility** → `$FOAM_USER_APPBIN/computeGradSchemeBlending`

Clean with:

```bash
./Allwclean
```

---

## Typical workflow

### 1. Load the library

In `system/controlDict`:

```c++
libs
(
    "libnumerics.so"
);
```

### 2. (Optional) Precompute blending weights

Copy the template dictionary into your case:

```bash
cp gradSchemeBlendingDict path/to/case/system/gradSchemeBlendingDict
```

Adjust quality thresholds as needed, then from the case directory:

```bash
computeGradSchemeBlending
```

This writes a volume field `gradSchemeBlending` (values in `[0, 1]`):

- `0` → prefer the background / least-squares side  
- `1` → prefer the blended / Gauss side  

Criteria are driven by mesh quality (aspect ratio, LSQ tensor conditioning, flat/chevron cells, non-orthogonality, skewness — depending on dictionary entries).

### 3. Select a scheme in `fvSchemes`

**Blend two arbitrary schemes** (dictionary form):

```c++
gradSchemes
{
    default
    {
        type    blendedGradScheme;
        scheme1 leastSquares;
        scheme2 Gauss linear;
    }
}
```

**Positional form:**

```c++
gradSchemes
{
    default     blendedGradScheme leastSquares Gauss linear;
}
```

**Hybrid least-squares (needs `gradSchemeBlending`):**

```c++
gradSchemes
{
    default     hybridLeastSquares;
    // or:      hybridLeastSquares1;
    // or:      hybridFaceLeastSquares;
}
```

**Runtime hybrids** (available on `main` / matching branches — no precomputed field required):

```c++
gradSchemes
{
    grad(U)     hybridGaussLS linear 0.5 0.3;
    //          hybridIterGaussLS linear 0.15 0.6;
}
```

---

## Repository layout

```text
OpenFoamNumerics/
├── Allwmake / Allwclean
├── computeGradSchemeBlending/   # preprocessing utility
├── gradSchemes/                 # scheme sources → libnumerics
├── gradSchemeBlendingDict       # template for system/
└── README.md
```

---

## Notes

- Precomputed blending is computed **once** from the mesh; remesh or change topology → rerun `computeGradSchemeBlending`.
- If `gradSchemeBlending` is missing, hybrid least-squares variants fall back to standard least-squares behaviour where implemented.
- Scheme availability always follows the **checked-out branch** and its `gradSchemes/Make/files`.

---

## License

Distributed under the terms of the GNU General Public License (GPL), consistent with OpenFOAM.
