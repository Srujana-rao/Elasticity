from __future__ import annotations

import argparse
import csv
import sys
from typing import List, Tuple, Dict

from .data import load_atomic_weights, load_ionic_radii
from .calcs import CompoundInput, compute_all


def read_inputs_csv(path: str) -> List[CompoundInput]:
	rows: List[CompoundInput] = []
	with open(path, newline="", encoding="utf-8") as f:
		reader = csv.DictReader(f)
		for row in reader:
			rows.append(
				CompoundInput(
					formula=row["compound"],
					coefficient=float(row["coefficient"]),
					density_g_cm3=float(row["density_g_cm3"]),
					Ui_value=float(row["Ui_value"]),
					Ui_units=row["Ui_units"],
				)
			)
	return rows


def prompt_interactive() -> List[CompoundInput]:
    print("Enter the number of compounds you want to add:")
    try:
        num = int(input("Number of compounds › ").strip())
    except ValueError:
        print("Invalid number.")
        return []
    items: List[CompoundInput] = []
    for i in range(num):
        print(f"\nEnter details for compound #{i+1}:")
        formula = input("  Compound formula (e.g., Bi2O3) › ").strip()
        coef = float(input("  Coefficient (mole units) › ").strip())
        density = float(input("  Density (g/cm^3) › ").strip())
        Ui_value = float(input("  Bond dissociation energy value (kJ/mol) › ").strip())
        # Units are not needed, so we set Ui_units to "kjmol"
        items.append(CompoundInput(formula, coef, density, Ui_value, "kjmol"))
    return items


def format_results(props, M, rho_mix, Vt, Gt, E, K, S, L, mu):
    print("\n=== Results ===")
    print(f"Molecular weight of total composition (M): {M:.3f} g/mol")
    print(f"Total density of mixture: {rho_mix:.6f} g/cm^3")
    print("")
    for p in props:
        print(f"Compound: {p.formula}")
        print(f"  MW: {p.MW:.3f} g/mol")
        print(f"  Coefficient mole fraction (Xi): {p.mole_fraction:.6f}")
        print(f"  Mass: {p.mass:.3f} g")
        print(f"  Volume: {p.volume:.3f} cm^3")
        print(f"  Density: {p.density:.6f} g/cm^3")
        print(f"  Oxidation (cation): +{p.oxidation_cation}")
        print(f"  Ionic radii: RA={p.RA_pm:.2f} pm, RO={p.RO_pm:.2f} pm")
        print(f"  V_i: {p.Vi:.2f} cm^3/mol")  # <-- changed to 2 decimals
        print(f"  G_i: {p.Gi:.3f} kJ/cm^3")   # <-- changed to 3 decimals
        print("")
    print(f"V_t (total packing density): {Vt:.6f}")
    print(f"G_t (total dissociation energy/volume): {Gt:.6f}")
    print("")
    print(f"E*: {E:.6f}")
    print(f"K*: {K:.6f}")
    print(f"S*: {S:.6f}")
    print(f"L*: {L:.6f}")
    print(f"μ*: {mu:.6f}")


def main(argv: List[str] | None = None) -> int:
	parser = argparse.ArgumentParser(description="Elastic constants of multi-component oxide systems")
	parser.add_argument("--inputs", help="CSV file with inputs", default=None)
	args = parser.parse_args(argv)

	atomic = load_atomic_weights()
	ionic = load_ionic_radii()
	ionic_map = {(k[0], k[1]): v.radius_pm for k, v in ionic.items()}

	if args.inputs:
		items = read_inputs_csv(args.inputs)
	else:
		items = prompt_interactive()

	props, M, rho_mix, Vt, Gt, E, K, S, L, mu = compute_all(items, atomic, ionic_map)
	format_results(props, M, rho_mix, Vt, Gt, E, K, S, L, mu)
	return 0


if __name__ == "__main__":
	sys.exit(main())


