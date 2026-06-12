# Bruggeman Effective Medium Calculator

This Python program calculates the effective refractive index (`n`) and extinction coefficient (`k`) of a mixed optical layer using the Bruggeman effective medium approximation.

The program reads one or two material TXT files, computes the effective complex permittivity, converts it back to `n` and `k`, and saves a new TXT file with the same column structure.

## Features

- Select one or two material TXT files through a file explorer dialog.
- If only one material is selected, the second material is assumed to be air (`n = 1`, `k = 0`).
- Choose the volume fraction of material 2.
- Choose the depolarization factor `q`, with `q = 0.5` as the default value.
- Automatically interpolates the second material if the two wavelength grids are different.
- Saves the mixed-layer optical constants as a TXT file.
- Selects the physically meaningful square-root branch by enforcing passive material behavior (`k >= 0`) and spectral continuity.

## Input File Format

Each material file must be a TXT file with one header row and three numerical columns:

```txt
wavelength (nm)    n    k
400                2.35 0.12
401                2.34 0.12
402                2.33 0.11
```

The script accepts whitespace-separated, tab-separated, comma-separated, or semicolon-separated numerical data.

## Method

The optical constants are first converted to complex permittivity:

```text
epsilon = (n + i k)^2
```

The effective permittivity is obtained using the symmetric Bruggeman equation:

```text
0 = f1 * (epsilon1 - epsilon_eff) / (q * epsilon1 + (1 - q) * epsilon_eff)
  + f2 * (epsilon2 - epsilon_eff) / (q * epsilon2 + (1 - q) * epsilon_eff)
```

where:

- `f2` is the volume fraction of material 2.
- `f1 = 1 - f2`.
- `q` is the depolarization factor.
- `epsilon1` and `epsilon2` are the complex permittivities of the input materials.
- `epsilon_eff` is the effective complex permittivity of the mixed layer.

Finally, the effective optical constants are obtained from:

```text
n + i k = sqrt(epsilon_eff)
```

The square-root branch is selected so that `k >= 0`. If the branch choice is numerically ambiguous, the program uses spectral continuity and chooses the solution closest to the previous wavelength point.

## Installation

Clone the repository and install the required dependency:

```bash
pip install -r requirements.txt
```

## Usage

Run the script:

```bash
python bruggeman.py
```

Then:

1. Select one or two material TXT files.
2. Enter the volume fraction of material 2.
3. Enter the depolarization factor `q`.
4. Choose the output filename and directory.

The output TXT file will contain:

```txt
wavelength (nm)    n    k
```

## Requirements

- Python 3.10 or newer
- NumPy

`tkinter` is used for the file dialogs and is included with most standard Python installations.
