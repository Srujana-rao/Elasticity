from __future__ import annotations

import re
from typing import Dict


SUBSCRIPT_MAP = str.maketrans({
	"₀": "0",
	"₁": "1",
	"₂": "2",
	"₃": "3",
	"₄": "4",
	"₅": "5",
	"₆": "6",
	"₇": "7",
	"₈": "8",
	"₉": "9",
})


ELEMENT_RE = re.compile(r"([A-Z][a-z]?)([0-9]*)")


def normalize_formula(formula: str) -> str:
	"""Convert Unicode subscripts to ASCII digits and remove whitespace.

	Examples: "B₂O₃" -> "B2O3"
	"""
	return formula.translate(SUBSCRIPT_MAP).replace(" ", "")


def parse_formula(formula: str) -> Dict[str, int]:
	"""Parse a simple inorganic formula without parenthesis, e.g., Bi2O3, SrO, Li2O.

	Returns dict of element -> count.
	"""
	text = normalize_formula(formula)
	index = 0
	counts: Dict[str, int] = {}
	for match in ELEMENT_RE.finditer(text):
		el, count_str = match.groups()
		count = int(count_str) if count_str else 1
		counts[el] = counts.get(el, 0) + count
		index = match.end()
	# basic validation
	if index != len(text):
		raise ValueError(f"Unparsed residue in formula '{formula}' at index {index}")
	return counts


def split_oxide_1cation(formula: str) -> tuple[str, int, int]:
	"""Given an oxide with one cation species and oxygen only, return (cation, x, y) where
	formula ~ A_x O_y. Raises if not a single-cation oxide.
	"""
	counts = parse_formula(formula)
	if "O" not in counts or len(counts) != 2:
		raise ValueError("Expected a single-cation oxide containing O only.")
	cation = next(e for e in counts.keys() if e != "O")
	return cation, counts[cation], counts["O"]


