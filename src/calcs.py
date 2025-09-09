from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .chem import parse_formula, split_oxide_1cation


AVOGADRO = 6.02214076e23  # 1/mol
PM_TO_CM = 1e-10  # 1 pm = 1e-12 m = 1e-10 cm


@dataclass
class CompoundInput:
	formula: str
	coefficient: float
	density_g_cm3: float
	Ui_value: float  # in kJ/mol or eV depending on Ui_units
	Ui_units: str  # 'kjmol' or 'ev'


@dataclass
class ComponentProps:
	formula: str
	MW: float
	mass: float
	volume: float
	density: float
	oxidation_cation: int
	RA_pm: float
	RO_pm: float
	Vi: float
	Gi: float
	mole_fraction: float


def infer_oxidation_single_cation_oxide(formula: str) -> Tuple[str, int, int, int]:
	"""Return (cation_symbol, x, y, cation_ox_state) for A_x O_y with O fixed at -2.
	Solves x*z + y*(-2) = 0 -> z = 2y/x. Requires integer z.
	"""
	cation, x, y = split_oxide_1cation(formula)
	val = 2 * y / x
	if abs(val - round(val)) > 1e-9:
		raise ValueError(f"Non-integer oxidation state inferred for {formula}: {val}")
	return cation, x, y, int(round(val))


def molecular_weight(formula: str, atomic_weights: Dict[str, float]) -> float:
	counts = parse_formula(formula)
	return sum(atomic_weights[el] * count for el, count in counts.items())


def ui_to_kj_per_mol(value: float, units: str) -> float:
	units_l = units.lower()
	if units_l == "kjmol":
		return value
	if units_l == "ev":
		# 1 eV/molecule = 96.485 kJ/mol
		return value * 96.485
	raise ValueError("Ui_units must be 'kjmol' or 'ev'")


def compute_all(
	inputs: List[CompoundInput],
	atomic_weights: Dict[str, float],
	ionic_radii_pm: Dict[Tuple[str, int], float],
) -> Tuple[List[ComponentProps], float, float, float, float, float, float, float, float, float]:
	# step 1 & 2: MW_i and total M
	MW_list = [molecular_weight(c.formula, atomic_weights) for c in inputs]
	M = sum(MW * c.coefficient for MW, c in zip(MW_list, inputs))

	# step 3: mass_i and volume_i, density_i provided
	component_masses = [c.coefficient * MW for c, MW in zip(inputs, MW_list)]
	component_volumes = [mass / c.density_g_cm3 for mass, c in zip(component_masses, inputs)]

	# step 4: mixture density
	Total_mass = sum(component_masses)
	Total_volume = sum(component_volumes)
	rho_mixture = Total_mass / Total_volume if Total_volume > 0 else float("nan")

	# step 8 needs mole fractions Xi
	moles_each = [c.coefficient for c in inputs]
	M_total_moles = sum(moles_each)
	Xi = [n / M_total_moles for n in moles_each]

	# steps 5,6,7,9 per component
	props: List[ComponentProps] = []
	Vi_list: List[float] = []
	Gi_list: List[float] = []
	for idx, c in enumerate(inputs):
		cat, x, y, z = infer_oxidation_single_cation_oxide(c.formula)
		RA_pm = ionic_radii_pm.get((cat, z))
		if RA_pm is None:
			raise KeyError(f"Missing ionic radius for {cat}^{{{z}}}")
		RO_pm = ionic_radii_pm.get(("O", -2))
		if RO_pm is None:
			raise KeyError("Missing ionic radius for O^{−2}")
		RA_cm = RA_pm * PM_TO_CM
		RO_cm = RO_pm * PM_TO_CM
		Vi = AVOGADRO * (4.0 / 3.0) * math.pi * (x * (RA_cm ** 3) + y * (RO_cm ** 3))
		Ui_kJ_mol = ui_to_kj_per_mol(c.Ui_value, c.Ui_units)
		Gi = (c.density_g_cm3 * Ui_kJ_mol) / MW_list[idx]
		Vi_list.append(Vi)
		Gi_list.append(Gi)
		props.append(
			ComponentProps(
				formula=c.formula,
				MW=MW_list[idx],
				mass=component_masses[idx],
				volume=component_volumes[idx],
				density=c.density_g_cm3,
				oxidation_cation=z,
				RA_pm=RA_pm,
				RO_pm=RO_pm,
				Vi=Vi,
				Gi=Gi,
				mole_fraction=Xi[idx],
			)
		)

	# step 8: total packing density Vt
	Vt = (M / rho_mixture) * sum(Vi * Xi_i for Vi, Xi_i in zip(Vi_list, Xi)) if rho_mixture > 0 else float("nan")

	# step 10: Gt
	Gt = sum(Gi * Xi_i for Gi, Xi_i in zip(Gi_list, Xi))

	# step 11: elastic constants
	E_star = 8.36 * Vt * Gt
	K_star = 10.0 * (Vt ** 2) * Gt
	S_star = (3 * K_star) / (10.2 * Vt - 1.0)
	L_star = K_star + (4.0 / 3.0) * S_star
	mu_star = 0.5 - 1.0 / (7.2 * Vt)

	return props, M, rho_mixture, Vt, Gt, E_star, K_star, S_star, L_star, mu_star


