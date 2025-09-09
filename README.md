## Elastic Constants Calculator for Multi-Component Oxides

This Python app computes elastic constants (E\*, K\*, S\*, L\*, μ\*) for mixtures of oxide compounds from user inputs, using periodic table data and ionic radii.

### Features
- Molecular weights from a CSV periodic table (118 elements supported; extendable)
- Formula parsing with Unicode subscripts support (e.g., `B₂O₃`)
- Oxidation state inference for single-cation oxides (oxygen fixed to −2)
- Ionic radii lookup by oxidation state (extendable CSV)
- Mixture density, component/total molar volumes, dissociation energy per unit volume
- Final elastic constants via provided empirical relations

### Project Layout
```
Physics/
  data/
    periodic_table.csv
    ionic_radii.csv
  src/
    data.py
    chem.py
    calcs.py
    cli.py
  README.md
```

### Install
No external dependencies are required (standard library only). Python 3.9+ recommended.

### Usage
Run the CLI with a CSV of inputs or interactively:

1) CSV input mode (recommended)
Create an input CSV with columns: `compound,coefficient,density_g_cm3,Ui_value,Ui_units` where `Ui_units` is `kjmol` or `ev`.

Example `inputs.csv`:
```
compound,coefficient,density_g_cm3,Ui_value,Ui_units
B2O3,1,2.46,806,kjmol
Bi2O3,1,8.9,439,kjmol
SrO,1,4.7,592,kjmol
Li2O,1,2.01,341,kjmol
```

Run:
```
python -m src.cli --inputs inputs.csv
```

2) Interactive mode
```
python -m src.cli
```

### Extending Data
- Add/adjust atomic weights in `data/periodic_table.csv`.
- Add ionic radii in `data/ionic_radii.csv`. Radii are in picometers (pm), and must specify oxidation state. Oxygen radius for O²⁻ is included.

### Notes
- Assumes simple oxides with a single cation species per compound (e.g., `Bi2O3`, `SrO`, `Li2O`, `B2O3`).
- Ionic radii usage: cation radius corresponds to inferred cation oxidation state; oxygen uses O²⁻.
- If a required radius or element is missing, the CLI will report the missing data and exit.


