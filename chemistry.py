import re
import pandas as pd
from dataclasses import dataclass
from typing import Dict
from constants import ISOTOPE_DELTAS, ATOMIC_MASSES

@dataclass
class MolecularFormula:
    """Molecular formula with exact mass calculation"""
    elements: Dict[str, int]

    def __post_init__(self):
        self.elements = {k: v for k, v in self.elements.items() if v > 0}

    def mass(self) -> float:
        return sum(ATOMIC_MASSES[element] * count
                  for element, count in self.elements.items())

    def __str__(self) -> str:
        formula = ""
        for element in ['C', 'H', 'N', 'O', 'P', 'S']:
            if element in self.elements and self.elements[element] > 0:
                count = self.elements[element]
                formula += f"{element}{count if count > 1 else ''}"
        return formula

def parse_isotope_label(isotope_str: str) -> Dict[str, int]:
    """
    Parse isotope label notation to extract isotope composition.

    Notation: D(2H), N15(15N), C13(13C), O18(18O)
    Examples:
        "M2DN15" -> {'2H': 2, '15N': 1}
        "M3D" -> {'2H': 3}
        "M4D2N15" -> {'2H': 4, '15N': 2}

    Raises:
        ValueError: If unrecognized characters found in isotope label
    """
    if not isotope_str or not isotope_str.strip():
        return {}

    isotope_dict = {}
    isotope_str = isotope_str.strip().upper()

    # Remove leading 'M' if present
    if isotope_str.startswith('M'):
        isotope_str = isotope_str[1:]

    i = 0
    while i < len(isotope_str):
        # Collect count digits
        count_str = ''
        while i < len(isotope_str) and isotope_str[i].isdigit():
            count_str += isotope_str[i]
            i += 1

        count = int(count_str) if count_str else 1

        if i >= len(isotope_str):
            break

        # Parse isotope type
        if isotope_str[i] == 'D':
            isotope_dict['2H'] = isotope_dict.get('2H', 0) + count
            i += 1
        elif isotope_str[i:i+3] == 'N15':
            isotope_dict['15N'] = isotope_dict.get('15N', 0) + count
            i += 3
        elif isotope_str[i:i+3] == 'C13':
            isotope_dict['13C'] = isotope_dict.get('13C', 0) + count
            i += 3
        elif isotope_str[i:i+3] == 'O18':
            isotope_dict['18O'] = isotope_dict.get('18O', 0) + count
            i += 3
        else:
            # Raise exception with message
            unrecognized = isotope_str[i:i+3] if i+3 <= len(isotope_str) else isotope_str[i:]
            raise ValueError(
                f"Unrecognized isotope label '{unrecognized}' at position {i} "
                f"in '{isotope_str}'. "
                f"Valid labels are: D (deuterium), N15, C13, O18"
            )

    return isotope_dict

def calculate_isotope_mass_shift(isotope_dict: Dict[str, int]) -> float:
    """
    Calculate total mass shift from isotope labeling.

    Example: {'2H': 2, '15N': 1} → 2×1.006277 + 1×0.997035 = 3.00959 Da
    """
    total_shift = 0.0
    for isotope, count in isotope_dict.items():
        if isotope in ISOTOPE_DELTAS:
            total_shift += ISOTOPE_DELTAS[isotope] * count
    return total_shift

def filter_isotope_token(token, formula):
        import re
        if pd.isna(token) or pd.isna(formula) or not str(token).strip():
            return token

        # Parse formula atoms (e.g., C16H32O -> {'C': 16, 'H': 32, 'O': 1})
        matches = re.findall(r'([A-Z][a-z]*)(\d*)', str(formula))
        atoms = {el: int(cnt) if cnt else 1 for el, cnt in matches}

        filtered = str(token).upper()

        # Strip labels if the formula lacks the required element
        if atoms.get('N', 0) == 0: filtered = re.sub(r'N15|15N', '', filtered)
        if atoms.get('C', 0) == 0: filtered = re.sub(r'\d*C13|\d*13C', '', filtered)
        if atoms.get('O', 0) == 0: filtered = re.sub(r'\d*O18|\d*18O', '', filtered)
        if atoms.get('H', 0) == 0: filtered = re.sub(r'\d*D', '', filtered)

        # If everything is stripped and only 'M' is left, return empty string
        return filtered if filtered != 'M' else ''
