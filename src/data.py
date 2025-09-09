from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Dict, Tuple


DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


def _csv_to_dict(path: str, key_field: str) -> Dict[str, dict]:
	with open(path, newline="", encoding="utf-8") as f:
		reader = csv.DictReader(f)
		out: Dict[str, dict] = {}
		for row in reader:
			key = row[key_field]
			out[key] = row
		return out


def load_atomic_weights() -> Dict[str, float]:
	path = os.path.join(DATA_DIR, "periodic_table.csv")
	rows = _csv_to_dict(path, "symbol")
	return {k: float(v["atomic_weight"]) for k, v in rows.items()}


@dataclass(frozen=True)
class IonicRadius:
	element: str
	oxidation: int
	radius_pm: float


def load_ionic_radii() -> Dict[Tuple[str, int], IonicRadius]:
	path = os.path.join(DATA_DIR, "ionic_radii.csv")
	with open(path, newline="", encoding="utf-8") as f:
		reader = csv.DictReader(f)
		out: Dict[Tuple[str, int], IonicRadius] = {}
		for row in reader:
			el = row["element"]
			ox = int(row["oxidation"])
			r = float(row["radius_pm"])
			out[(el, ox)] = IonicRadius(el, ox, r)
		return out


