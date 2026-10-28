# Einstein Field Equations from the Einstein-Hilbert Action (xAct)

Mathematica / xAct code that derives the Einstein field equations by varying the action

S = ∫ d⁴x √(−g) [ R/(2κ) + L_matter ],  κ = 8πG (c = 1), signature (−,+,+,+)

with **no cosmological constant**. Matter sectors included: minimally coupled scalar field with potential V(φ), and Maxwell field. The code also computes the energy-momentum tensor T_ab = −(2/√−g) δ(√−g L_m)/δg^{ab}, the matter field equations, and consistency checks (Bianchi identity).

## Repository contents

| File | Description |
|------|-------------|
| `EFE_variation_xAct.nb` | Notebook version (recommended), run cell by cell |
| `EFE_variation_xAct.wl` | Plain script version |
| `README.md` | This file |

## Requirements

| Requirement | Version / note |
|-------------|----------------|
| Wolfram Mathematica | 10 or newer (12+ recommended), or the free Wolfram Engine for developers |
| xAct package | latest release from <http://www.xact.es> (uses `xTensor`, which loads `xCore` and `xPerm` automatically) |
| OS | Windows, macOS, or Linux (xPerm binary is bundled with xAct) |

xAct is not included in this repo (it is GPL-licensed and distributed separately). Install it yourself as below.

## Setup

### 1. Install xAct (one time)

1. Download the latest xAct zip from <http://www.xact.es>.
2. In Mathematica, find your applications folder:
   ```mathematica
   FileNameJoin[{$UserBaseDirectory, "Applications"}]
   ```
3. Unzip so that a folder named `xAct` sits directly inside it, for example `~/.Mathematica/Applications/xAct/` on Linux.
4. Restart Mathematica and test:
   ```mathematica
   << xAct`xTensor`
   ```
   A license banner should print with no errors.

### 2. Get the code

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
```

### 3. Run

**Notebook (recommended):** open `EFE_variation_xAct.nb` in Mathematica and press `Shift+Enter` on each input cell from top to bottom, in one kernel session. The first cell (`Quit[]`) is optional and only for a fresh kernel. Do not use *Evaluate Notebook*, as it would run `Quit[]` mid-run.

**Script from the command line** (needs `wolframscript`):

```bash
wolframscript -file EFE_variation_xAct.wl
```

Only `Print[...]` output is shown in script mode, so wrap any expression you want to see in `Print[...]`.

## Expected output

- Field equations: `G_ab − κ T_ab = 0`
- Scalar field equation: `□φ − V′(φ) = 0`
- Maxwell equation: `∇_b F^{ba} = 0`
- Every cell marked "expect 0" should return `0`.

## Customizing

Change the matter sector with `LMatter = LScalar`, `LMaxwell`, or their sum. For modified gravity, edit `LGrav` (for example `fR[RicciScalarCD[]]/(2 kap)` for f(R)) and re-run from the Lagrangian section onward.

## Troubleshooting

- `<< xAct`xTensor`` fails: the `xAct` folder is not directly inside `Applications`, or Mathematica was not restarted.
- "Already defined" errors when re-running: run `Quit[]` in a separate cell, then run everything again.
- `VarD` problem with `Detg[]` on an older xAct version: update xAct to the latest release.

## Citation

If you use this code, please also cite xAct: J.-M. Martín-García, *xAct: Efficient tensor computer algebra for the Wolfram Language*, <http://www.xact.es>.

## Author

Ankit C Khunt
