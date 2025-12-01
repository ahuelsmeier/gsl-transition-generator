#!/usr/bin/env python3

"""
GSL + Ceramide Transition Generator - Version 1.0

author: Andreas J. Hülsmeier
copyright: Copyright 2025, Andreas J. Hülsmeier / University of Zurich, University Hospital Zurich
license: MIT
version: 1.0.1
maintainer: Andreas J. Hülsmeier
email: andreas.huelsmeier@uzh.ch
status: Prototype
"""

import pandas as pd
import re
import json
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
import argparse
import sys
import os
import logging

# Configure logging
log_level = os.getenv('LOG_LEVEL', 'INFO').upper()
logging.basicConfig(
    level=getattr(logging, log_level),
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Create a logger for this module
logger = logging.getLogger(__name__)

# Atomic masses (IUPAC 2016)
ATOMIC_MASSES = {
    'C': 12.0000000,
    'H': 1.00782503223,
    'N': 14.0030740048,
    'O': 15.9949146223,
    'P': 30.9737619985,
    'S': 31.9720711744
}

PROTON_MASS = 1.007276466812  # Mass of H+


# ============================================================================
# ISOTOPE LABELING WITH M/Z CALCULATION
# ============================================================================

# Exact mass differences for isotope substitutions (IUPAC values)
ISOTOPE_DELTAS = {
    '2H': 1.006276746,      # ²H - ¹H (Deuterium)
    '15N': 0.997034893,     # ¹⁵N - ¹⁴N
    '13C': 1.0033548378,    # ¹³C - ¹²C
    '18O': 2.004244         # ¹⁸O - ¹⁶O (optional)
}

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


class ConfigManager:
    """Manage user configuration for LCB and fatty acid selections"""

    CONFIG_FILE = "gsl_config.json"

    DEFAULT_CONFIG = {
        "lcb_selections": {
            "standard": ["18:0;2", "18:1;2", "18:2;2"],
            "doxCer": ["18:0;1", "18:1;1"]
        },
        "fatty_acid_range": {
            "min_length": 16,
            "max_length": 26,
            "unsaturations": [0, 1],
            "even_chain_only": False
        },
        "selected_fatty_acids": None
    }

    @classmethod
    def load_config(cls):
        """Load configuration from file, or create default"""
        config_path = Path(cls.CONFIG_FILE)
        if config_path.exists():
            try:
                with open(config_path, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, FileNotFoundError, KeyError):
                return cls.DEFAULT_CONFIG.copy()
        return cls.DEFAULT_CONFIG.copy()

    @classmethod
    def save_config(cls, config):
        """Save configuration to file"""
        with open(cls.CONFIG_FILE, 'w') as f:
            json.dump(config, f, indent=2)

    @classmethod
    def get_lcb_list(cls, lipid_class="standard"):
        """Get LCB list for given lipid class"""
        config = cls.load_config()
        key = "doxCer" if lipid_class == "doxCer" else "standard"
        return config["lcb_selections"][key]

    @classmethod
    def get_fatty_acid_list(cls):
        """Get fatty acid list based on configuration"""
        config = cls.load_config()

        if config.get("selected_fatty_acids"):
            return config["selected_fatty_acids"]

        fa_config = config["fatty_acid_range"]
        fatty_acids = []

        # Get chain length range
        min_len = fa_config["min_length"]
        max_len = fa_config["max_length"]
        even_only = fa_config.get("even_chain_only", False)

        # Generate fatty acid list
        for length in range(min_len, max_len + 1):
            # Skip odd-chain lengths if even_only is True
            if even_only and length % 2 != 0:
                continue

            for unsat in fa_config["unsaturations"]:
                fatty_acids.append(f"{length}:{unsat}")

        return fatty_acids

    @classmethod
    def get_default_config(cls):
        """Get complete default configuration including UI state"""
        return {
            "lcb_selections": {
                "standard": ["18:0;2", "18:1;2", "18:2;2"],
                "doxCer": ["18:0;1", "18:1;1"]
            },
            "fatty_acid_range": {
                "min_length": 16,
                "max_length": 26,
                "unsaturations": [0, 1],
                "even_chain_only": False
            },
            "selected_fatty_acids": None,
            "charge_states": [1],
            "selected_adducts": ["M+H", "M-H"],
            "isotope_labeling": {
                "enabled": False,
                "gsl_isotope": "M2DN15",
                "cer_isotope": "M2DN15",
                "doxcer_isotope": "M3D",
                "label_keywords": "LCB,precursor,HG(-Hex",
                "blank_mz": False
            }
        }



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


class LipidDatabase:
    """Lipid database - ALL formulas verified or calculated from verified structures"""
    # (dehydrated residues)
    GSL_HEADGROUP_COMPOSITIONS = {
        'Hex': {'C': 6, 'H': 10, 'O': 5},
        'SM4': {'C': 6, 'H': 10, 'N': 0, 'O': 8, 'S': 1},
        'Lac': {'C': 12, 'H': 20, 'O': 10},
        'LC3': {'C': 20, 'H': 33, 'N': 1, 'O': 15},
        'LC4': {'C': 26, 'H': 43, 'N': 1, 'O': 20},
        'Gb3': {'C': 18, 'H': 30, 'N': 0, 'O': 15},
        'Gb4': {'C': 26, 'H': 43, 'N': 1, 'O': 20},
        'GA2': {'C': 20, 'H': 33, 'N': 1, 'O': 15},
        'GA1': {'C': 26, 'H': 43, 'N': 1, 'O': 20},
        'GM4': {'C': 17, 'H': 27, 'N': 1, 'O': 13},
        'GM3': {'C': 23, 'H': 37, 'N': 1, 'O': 18},
        'GM2': {'C': 31, 'H': 50, 'N': 2, 'O': 23},
        'GM1': {'C': 37, 'H': 60, 'N': 2, 'O': 28},
        'GD3': {'C': 34, 'H': 54, 'N': 2, 'O': 26},
        'GD2': {'C': 42, 'H': 67, 'N': 3, 'O': 31},
        'GD1a': {'C': 48, 'H': 77, 'N': 3, 'O': 36},
        'GD1b': {'C': 48, 'H': 77, 'N': 3, 'O': 36},
        'GT3': {'C': 45, 'H': 71, 'N': 3, 'O': 34},
        'GT2': {'C': 53, 'H': 84, 'N': 4, 'O': 39},
        'GT1a': {'C': 59, 'H': 94, 'N': 4, 'O': 44},
        'GT1b': {'C': 59, 'H': 94, 'N': 4, 'O': 44},
        'GT1c': {'C': 59, 'H': 94, 'N': 4, 'O': 44},
        'GQ1': {'C': 70, 'H': 111, 'N': 5, 'O': 52},
        'GP1': {'C': 81, 'H': 128, 'N': 6, 'O': 60},
        'nLc10': {'C': 68, 'H': 112, 'N': 4, 'O': 50},
        'nLc8': {'C': 54, 'H': 89, 'N': 3, 'O': 40},
        'nLc6': {'C': 40, 'H': 66, 'N': 2, 'O': 30},
        'SHex2': {'C': 12, 'H': 20, 'O': 13, 'S': 1},
    }

    CERAMIDE_COMPOSITIONS = {
        'Cer': {'C': 0, 'H': 0, 'N': 0, 'O': 0},
        'doxCer': {'C': 0, 'H': 0, 'N': 0, 'O': 0},
    }

    SM_COMPOSITIONS = {
        'SM': {'C': 5, 'H': 12, 'N': 1, 'O': 3, 'P': 1},
    }

    @classmethod
    def get_lipid_composition(cls, lipid_class: str) -> Dict[str, int]:
        if lipid_class in cls.GSL_HEADGROUP_COMPOSITIONS:
            return cls.GSL_HEADGROUP_COMPOSITIONS[lipid_class].copy()
        elif lipid_class in cls.CERAMIDE_COMPOSITIONS:
            return cls.CERAMIDE_COMPOSITIONS[lipid_class].copy()
        elif lipid_class in cls.SM_COMPOSITIONS:
            return cls.SM_COMPOSITIONS[lipid_class].copy()
        else:
            raise ValueError(f"Unknown lipid class: {lipid_class}")

    @classmethod
    def get_all_classes(cls) -> List[str]:
        gsl_classes = list(cls.GSL_HEADGROUP_COMPOSITIONS.keys())
        ceramide_classes = list(cls.CERAMIDE_COMPOSITIONS.keys())
        sm_classes = list(cls.SM_COMPOSITIONS.keys())
        return sorted(gsl_classes + ceramide_classes + sm_classes)

    @classmethod
    def is_ceramide_class(cls, lipid_class: str) -> bool:
        return lipid_class in cls.CERAMIDE_COMPOSITIONS

    @classmethod
    def is_gsl_class(cls, lipid_class: str) -> bool:
        return lipid_class in cls.GSL_HEADGROUP_COMPOSITIONS

    @classmethod
    def molecular_weight_range(cls, lipid_class: str) -> str:
        mw_ranges = {
            'Hex': '600-900',
            'Lac': '700-800', 'LC3': '900-1000', 'LC4': '1100-1200',
            'Gb3': '1000-1100', 'Gb4': '1200-1400',
            'GA2': '1000-1100', 'GA1': '1200-1400',
            'GM4': '1000-1200',
            'GM3': '1200-1300', 'GM2': '1400-1500', 'GM1': '1500-1600',
            'GD3': '1500-1600', 'GD2': '1700-1800', 'GD1a': '1800-1900', 'GD1b': '1800-1900',
            'GT3': '1700-1900', 'GT2': '1900-2100', 'GT1a': '2000-2400', 'GT1b': '2000-2400', 'GT1c': '2000-2400',
            'GQ1': '2300-2400', 'GP1': '2600-2700',
            'Cer': '500-750', 'doxCer': '480-730',
            'SM': '650-900', 'SM4': '700-1100',
            'nLc10': '2200-2500',
            'nLc8': '1900-2200',
            'nLc6': '1500-1800',
            'SHex2': '700-1100',
        }
        return mw_ranges.get(lipid_class, '500-2000')

    @classmethod
    def get_sialic_acid_count(cls, lipid_class: str) -> int:
        sialic_counts = {
            'GM4': 1, 'GM3': 1, 'GM2': 1, 'GM1': 1,
            'GD3': 2, 'GD2': 2, 'GD1a': 2, 'GD1b': 2,
            'GT3': 3, 'GT2': 3, 'GT1a': 3, 'GT1b': 3, 'GT1c': 3,
            'GQ1': 4, 'GP1': 5
        }
        return sialic_counts.get(lipid_class, 0)

    @classmethod
    def get_structure_description(cls, lipid_class: str) -> str:
        structures = {
            'doxCer': 'headless, 1-deoxy-Ceramide',
            'Cer': 'HO-',
            'Hex': 'β-D-Glc- or β-D-Gal-linked Ceramide (Hexosylceramide)',
            'SM4': '3-O-sulfated Gal-Cer (Sulfatide)',
            'Lac': 'Galβ1-4Glc-Cer (Lactosyl-Ceramide)',
            'LC3': 'GlcNAcβ1-3Galβ1-4Glc-Cer (Lacto/neoLacto-series), isobaric to GA2',
            'LC4': 'Galβ1-3GlcNAcβ1-3Galβ1-4Glc-Cer (Lacto-series)',
            'Gb3': 'Galα1-4Galβ1-4Glc-Cer (Globotriaosylceramide)',
            'Gb4': 'GalNAcβ1-3Galα1-4Galβ1-4Glc-Cer (isobaric to GA1)',
            'GA2': 'GalNAcβ1-4Galβ1-4Glc-Cer (asialo-GM2), isobaric to Lc3',
            'GA1': 'Galβ1-3GalNAcβ1-4Galβ1-4Glc-Cer (asialo-GM1)',
            'GM4': 'Neu5Acα2-3Galβ-Cer',
            'GM3': 'NeuAcα2-3Galβ1-4Glcβ-Cer',
            'GM2': 'GalNAcβ1-4(NeuAcα2-3)Galβ1-4Glcβ-Cer',
            'GM1': 'Galβ1-3GalNAcβ1-4(NeuAcα2-3)Galβ1-4Glcβ-Cer',
            'GD3': 'NeuAcα2-8NeuAcα2-3Galβ1-4Glcβ-Cer',
            'GD2': 'GalNAcβ1-4(NeuAcα2-8NeuAcα2-3)Galβ1-4Glcβ-Cer',
            'GD1a': 'NeuAcα2-3Galβ1-3GalNAcβ1-4(NeuAcα2-3)Galβ1-4Glcβ-Cer',
            'GD1b': 'Galβ1-3GalNAcβ1-4(NeuAcα2-8NeuAcα2-3)Galβ1-4Glcβ-Cer',
            'GT3': 'Neu5Acα2-8Neu5Acα2-8Neu5Acα2-3Galβ1-4Glcβ-Cer',
            'GT2': 'GalNAcβ1-4(Neu5Acα2-8Neu5Acα2-8Neu5Acα2-3)Galβ1-4Glcβ-Cer',
            'GT1a': 'Neu5Acα2-8Neu5Acα2-3Galβ1-3GalNAcβ1-4(Neu5Acα2-3)Galβ1-4Glcβ-Cer',
            'GT1b': 'Neu5Acα2-3Galβ1-3GalNAcβ1-4(Neu5Acα2-8Neu5Acα2-3)Galβ1-4Glcβ-Cer',
            'GT1c': 'Galβ1-3GalNAcβ1-4(NeuAcα2-8NeuAcα2-8NeuAcα2-3)Galβ1-4Glcβ-Cer',
            'SM': 'Phosphocholine-Cer (Sphingomyelin)',
            'nLc10': 'GlcNAcβ1-3Galβ1-4GlcNAcβ1-3(Galα1-3Galβ1-4GlcNAcβ1-6)Galβ1-4GlcNAcβ1-3Galβ1-4Glcβ-Cer',
            'nLc8': 'Galβ1-4GlcNAcβ1-3(Galβ1-4GlcNAcβ1-6)Galβ1-4GlcNAcβ1-3Galβ1-4Glcβ-Cer',
            'nLc6': 'Galβ1-4GlcNAcβ1-3Galβ1-4GlcNAcβ1-3Galβ1-4Glcβ-Cer',
            'SHex2': 'Sulfated dihexosylceramide',
        }
        return structures.get(lipid_class, 'Structure not specified')



class GSLFragmentRules:
    """Fragment rules for GSL classes, x-, y-, z- fragmets, aglycone fragments"""

    @staticmethod
    def get_ceramide_fragments(lcb_type: str):
        """
        Reuse ceramide LCB fragment rules.
        GSLs use standard dihydroxy LCB fragments (never doxCer).
        """
        return CeramideFragmentRules.get_ceramide_fragments(lcb_type, is_doxcer=False)

    @staticmethod
    def get_hex_hg_loss_fragments(precursor_formula):
        """
        Generate Hex headgroup loss fragments by subtracting hexose units.
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        HEX_HEADGROUP_LOSSES = {
            "HG(-Hex,162)": {'C': 6, 'H': 10, 'O': 5},
            "HG(-Hex,180)": {'C': 6, 'H': 12, 'O': 6},
            "HG(-Hex,198)": {'C': 6, 'H': 14, 'O': 7}
        }
        result = []
        for frag_name, loss_dict in HEX_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result

    @staticmethod
    def get_lac_hg_loss_fragments(precursor_formula):
        """
        Generate Lac headgroup loss fragments by subtracting lactose units.
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        LAC_HEADGROUP_LOSSES = {
            "HG(-Hex2,342)": {'C': 12, 'H': 22, 'O': 11},
            "HG(-Hex2,360)": {'C': 12, 'H': 24, 'O': 12},
            "HG(-Hex2,324)": {'C': 12, 'H': 20, 'O': 10},
        }
        result = []
        for frag_name, loss_dict in LAC_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result

    @staticmethod
    def get_gb3_hg_loss_fragments(precursor_formula):
        """
        Generate Gb3 headgroup loss fragments by subtracting glycan units.
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GB3_HEADGROUP_LOSSES = {
            "HG(-Hex3,504)": {'C': 18, 'H': 32, 'O': 16},   # Hex3
            "HG(-Hex3,522)": {'C': 18, 'H': 34, 'O': 17},   # Hex3 -H2O
            "HG(-Hex3,540)": {'C': 18, 'H': 36, 'O': 18},   # Hex3 -2H2O
            "HG(-Hex2,342)": {'C': 12, 'H': 22, 'O': 11},   # Hex2
            "HG(-Hex,180)": {'C': 6, 'H': 12, 'O': 6},      # Hex
        }
        result = []
        for frag_name, loss_dict in GB3_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result

    @staticmethod
    def get_gb4_hg_loss_fragments(precursor_formula):
        """
        Generate Gb4 headgroup loss fragments by subtracting glycan units.
        Gb4 = GalNAcβ1-3Galα1-4Galβ1-4Glc-Cer
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GB4_HEADGROUP_LOSSES = {
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},           # HexNAc
            "HG(-HexNAcHex,383)": {'C': 14, 'H': 25, 'N': 1, 'O': 11},      # HexNAc + Hex
            "HG(-HexNAcHex2,545)": {'C': 20, 'H': 35, 'N': 1, 'O': 16},     # HexNAc + 2Hex
            "HG(-HexNAcHex3,707)": {'C': 26, 'H': 45, 'N': 1, 'O': 21},     # HexNAc + 3Hex
            "HG(-HexNAcHex3,725)": {'C': 26, 'H': 47, 'N': 1, 'O': 22},     # HexNAc + 3Hex, -H2O
        }
        result = []
        for frag_name, loss_dict in GB4_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result

    @staticmethod
    def get_ga1_hg_loss_fragments(precursor_formula):
        """
        Generate GA1 headgroup loss fragments by subtracting glycan units.
        GA1 (asialo-GM1) = GalNAcβ1-4Galβ1-3GalNAcβ1-4Galβ1-4Glc-Cer
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GA1_HEADGROUP_LOSSES = {
            # Sequential losses from terminal to reducing end
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},           # HexNAc
            "HG(-HexNAcHex,383)": {'C': 14, 'H': 25, 'N': 1, 'O': 11},      # HexNAc + Hex
            "HG(-HexNAc2Hex,586)": {'C': 22, 'H': 38, 'N': 2, 'O': 16},     # 2 HexNAc + Hex
            "HG(-HexNAc2Hex,604)": {'C': 22, 'H': 40, 'N': 2, 'O': 17},     # 2 HexNAc + Hex, -H2O
            "HG(-HexNAc2Hex2,748)": {'C': 28, 'H': 48, 'N': 2, 'O': 21},    # 2 HexNAc + 2 Hex
            "HG(-HexNAc2Hex2,766)": {'C': 28, 'H': 50, 'N': 2, 'O': 22},    # 2 HexNAc + 2 Hex, -H2O
            "HG(-HexNAc2Hex3,910)": {'C': 34, 'H': 58, 'N': 2, 'O': 26},    # Pentasaccharide
            "HG(-HexNAc2Hex3,928)": {'C': 34, 'H': 60, 'N': 2, 'O': 27},    # Pentasaccharide, -H2O
        }
        result = []
        for frag_name, loss_dict in GA1_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result

    @staticmethod
    def get_ga2_hg_loss_fragments(precursor_formula):
        """
        Generate GA2 headgroup loss fragments by subtracting glycan units.
        GA2 (asialo-GM2) = GalNAcβ1-4Galβ1-4Glc-Cer
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GA2_HEADGROUP_LOSSES = {
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},
            "HG(-HexNAcHex,383)": {'C': 14, 'H': 25, 'N': 1, 'O': 11},
            "HG(-HexNAcHex2,545)": {'C': 20, 'H': 35, 'N': 1, 'O': 16},
            "HG(-HexNAcHex2,563)": {'C': 20, 'H': 37, 'N': 1, 'O': 17},
        }
        result = []
        for frag_name, loss_dict in GA2_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result

    @staticmethod
    def get_lc3_hg_loss_fragments(precursor_formula):
        """
        Generate LC3 headgroup loss fragments by subtracting glycan units.
        LC3 (lactotriaosylceramide) = GalNAcβ1-3Galβ1-4Glc-Cer
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        LC3_HEADGROUP_LOSSES = {
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},           # Single HexNAc
            "HG(-HexNAcHex,383)": {'C': 14, 'H': 25, 'N': 1, 'O': 11},      # HexNAc + Hex
            "HG(-HexNAcHex,401)": {'C': 14, 'H': 27, 'N': 1, 'O': 12},      # HexNAc + Hex, -H2O
            "HG(-HexNAcHex2,545)": {'C': 20, 'H': 35, 'N': 1, 'O': 16},     # Trisaccharide
            "HG(-HexNAcHex2,563)": {'C': 20, 'H': 37, 'N': 1, 'O': 17},     # Trisaccharide, -H2O
        }
        result = []
        for frag_name, loss_dict in LC3_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result

    @staticmethod
    def get_lc4_hg_loss_fragments(precursor_formula):
        """
        Generate LC4 headgroup loss fragments by subtracting glycan units.
        LC4 (lactotetraosylceramide) = GalNAcβ1-3GalNAcβ1-3Galβ1-4Glc-Cer
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        LC4_HEADGROUP_LOSSES = {
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},           # Single HexNAc
            "HG(-HexNAc2,442)": {'C': 16, 'H': 28, 'N': 2, 'O': 11},        # 2 HexNAc
            "HG(-HexNAc2Hex,586)": {'C': 22, 'H': 38, 'N': 2, 'O': 16},     # 2 HexNAc + Hex
            "HG(-HexNAc2Hex,604)": {'C': 22, 'H': 40, 'N': 2, 'O': 17},     # 2 HexNAc + Hex, -H2O
            "HG(-HexNAc2Hex2,748)": {'C': 28, 'H': 48, 'N': 2, 'O': 21},    # Tetrasaccharide
            "HG(-HexNAc2Hex2,766)": {'C': 28, 'H': 50, 'N': 2, 'O': 22},    # Tetrasaccharide, -H2O
            "HG(-HexNAc2Hex2,784)": {'C': 28, 'H': 52, 'N': 2, 'O': 23},    # Tetrasaccharide, -2H2O
        }
        result = []
        for frag_name, loss_dict in LC4_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result

    @staticmethod
    def get_sm4_hg_loss_fragments(precursor_formula):
        """
        Generate SM4 headgroup loss fragments by subtracting sulfated hexose units.
        SM4 (sulfatide) = 3-sulfogalactosylceramide (SO3-Galβ1-Cer)
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        SM4_HEADGROUP_LOSSES = {
            "HG(-SO4H2,98)": {'H': 2, 'S': 1, 'O': 4},           # H2SO4 loss
            "HG(-SHex,260)": {'C': 6, 'H': 12, 'S': 1, 'O': 9},  # Sulfated hexose, -H2O
            "HG(-SHex,278)": {'C': 6, 'H': 14, 'S': 1, 'O': 10}, # Sulfated hexose
            "HG(-SHex,242)": {'C': 6, 'H': 10, 'S': 1, 'O': 8},  # Sulfated hexose, -2H2O
        }
        result = []
        for frag_name, loss_dict in SM4_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result


    @staticmethod
    def get_shex2_hg_loss_fragments(precursor_formula):
        """
        Generate SHex2 headgroup loss fragments.
        SHex2 = monosulfated dihexosylceramide (one sulfate on 2 hexoses)
        Common structure: SO3-Gal-Glc-Cer
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        SHEX2_HEADGROUP_LOSSES = {
            "HG(-SO3,80)": {'S': 1, 'O': 3},
            "HG(-HSO3,81)": {'H': 1, 'S': 1, 'O': 3},
            "HG(-H2SO4,98)": {'H': 2, 'S': 1, 'O': 4},
            "HG(-SHex,242)": {'C': 6, 'H': 10, 'S': 1, 'O': 8},
            "HG(-SHex,260)": {'C': 6, 'H': 12, 'S': 1, 'O': 9},
            "HG(-SHex,278)": {'C': 6, 'H': 14, 'S': 1, 'O': 10},
            "HG(-SHexHex,404)": {'C': 12, 'H': 20, 'S': 1, 'O': 13},
            "HG(-SHexHex,422)": {'C': 12, 'H': 22, 'S': 1, 'O': 14},
            "HG(-SHexHex,440)": {'C': 12, 'H': 24, 'S': 1, 'O': 15},
        }
        result = []
        for frag_name, loss_dict in SHEX2_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            frag_formula_str = str(frag_formula)
            frag_mass = frag_formula.mass()
            result.append((frag_name, frag_formula_str, frag_mass))
        return result


    @staticmethod
    def get_gm4_hg_loss_fragments(precursor_formula):
        """
        Generate GM4 headgroup loss fragments (Neu5Ac on GalCer).
        Returns a list of (name, formula_str, monoisotopic_mass).
        """
        GM4_HEADGROUP_LOSSES = {
            # Neu5Ac (residue and free acid)
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},      # Neu5Ac
            # Entire GM4 headgroup
            "HG(-HexNeu5Ac,471)": {'C': 17, 'H': 29, 'N': 1, 'O': 14},
        }
        result = []
        for frag_name, loss_dict in GM4_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_gm3_hg_loss_fragments(precursor_formula):
        """
        Generate GM3 headgroup loss fragments.
        GM3: NeuAc(α2-3)Gal(β1-4)Glc(β)-Cer

        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GM3_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},      # Neu5Ac
            'HG(-HexNeu5Ac,471)': {'C': 17, 'H': 29, 'N': 1, 'O': 14},  # Neu5Ac + Hex
            'HG(-Hex2Neu5Ac,633)': {'C': 23, 'H': 39, 'N': 1, 'O': 19},  # Entire GM3 headgroup
        }

        result = []
        for frag_name, loss_dict in GM3_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))

        return result

    @staticmethod
    def get_gm2_hg_loss_fragments(precursor_formula):
        """
        Generate selected GM2 headgroup loss fragments:
        -Neu5Ac, -Neu5Ac-HexNAc, -Neu5Ac-Hex-HexNAc, -Neu5Ac-HexNAc-Hex-Hex
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GM2_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},      # Neu5Ac
            "HG(-Neu5AcHexNAc,512)": {'C': 19, 'H': 32, 'N': 2, 'O': 14},
            "HG(-Neu5AcHexNAcHex,674)": {'C': 25, 'H': 42, 'N': 2, 'O': 19},
            "HG(-Neu5AcHexNAcHex2,836)": {'C': 31, 'H': 52, 'N': 2, 'O': 24}, # not the dehydrated HG
        }

        result = []
        for frag_name, loss_dict in GM2_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_gm1_hg_loss_fragments(precursor_formula):
        """
        Generate GM1 headgroup loss fragments:
        GM1: Galβ1-3GalNAcβ1-4(NeuAcα2-3)Galβ1-4Glc-Cer
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GM1_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},      # Neu5Ac
            "HG(-Neu5AcHexNAcHex,674)": {'C': 25, 'H': 42, 'N': 2, 'O': 19}, # not the dehydrated HG
            "HG(-Neu5AcHexNAcHex2,836)": {'C': 31, 'H': 52, 'N': 2, 'O': 24}, # not the dehydrated HG
            "HG(-Neu5AcHexNAcHex3,998)": {'C': 37, 'H': 62, 'N': 2, 'O': 29}, # Entire GM1 headgroup, not dehydrated
            "HG(-Neu5AcHexNAcHex3,1016)": {'C': 37, 'H': 64, 'N': 2, 'O': 30}, # Entire GM1 headgroup, not dehydrated, plus H2O
        }

        result = []
        for frag_name, loss_dict in GM1_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_gd3_hg_loss_fragments(precursor_formula):
        """
        Generate GD3 headgroup loss fragments.
        GD3: NeuAcα2-8NeuAcα2-3Galβ1-4Glcβ-Cer (two sialic acids)
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GD3_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},  # One Neu5Ac
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc
            "HG(-HexNeu5Ac2,780)": {'C': 28, 'H': 46, 'N': 2, 'O': 22},  # 2 Neu5Ac + Hex
            "HG(-Hex2Neu5Ac2,942)": {'C': 34, 'H': 56, 'N': 2, 'O': 27},  # Entire GD3 headgroup
        }

        result = []
        for frag_name, loss_dict in GD3_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))

        return result
    @staticmethod
    def get_gd2_hg_loss_fragments(precursor_formula):
        """
        Generate GD2 headgroup loss fragments.
        GD2: GalNAcβ1-4(NeuAcα2-8NeuAcα2-3)Galβ1-4Glcβ-Cer
        Returns a list of (name, formula_str, monoisotopic_mass)
        """
        GD2_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},  # One Neu5Ac
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},  # Single HexNAc
            "HG(-HexNAc,203)": {'C': 8, 'H': 13, 'N': 1, 'O': 5},  # Single HexNAc -H2O, Y2β
            "HG(-Neu5AcHexNAc,512)": {'C': 19, 'H': 32, 'N': 2, 'O': 14},  # Neu5Ac + HexNAc
            "HG(-Neu5Ac2HexNAc,803)": {'C': 30, 'H': 49, 'N': 3, 'O': 22},  # 2 NeuAc + HexNAc
            "HG(-Neu5Ac2HexNAcHex,965)": {'C': 36, 'H': 59, 'N': 3, 'O': 27},  # Hex + HexNAc+ 2 Neu5Ac
            "HG(-HexNAcHex2Neu5Ac2,1127)": {'C': 42, 'H': 69, 'N': 3, 'O': 32},  # Full HG
        }

        result = []
        for frag_name, loss_dict in GD2_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))

        return result

    @staticmethod
    def get_gd1_hg_loss_fragments(precursor_formula, is_a=True):
        """
        Generate GD1a/GD1b headgroup loss fragments.
        Set is_a=True for GD1a (disialo fragments appear); False for GD1b.
        """
        GD1_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},
            "HG(-Neu5Ac2,618)": {'C': 22, 'H': 38, 'N': 2, 'O': 18},   # 2 x Neu5Ac
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # Neu5Acα2-8NeuAc
            "HG(-Neu5Ac2Hex,762)": {'C': 28, 'H': 46, 'N': 2, 'O': 22},  # Neu5Acα2-8NeuAc, Hex
            "HG(-Neu5Ac2HexNAcHex,965)": {'C': 36, 'H': 59, 'N': 3, 'O': 27},  # Neu5Ac2, HexNAc, Hex
            "HG(-Neu5Ac2HexNAcHex2,1127)": {'C': 42, 'H': 69, 'N': 3, 'O': 32},  # Neu5Ac2, HexNAc, Hex2
            "HG(-Neu5Ac2HexNAcHex3,1289)": {'C': 48, 'H': 79, 'N': 3, 'O': 37},  # Neu5Ac2, HexNAc, Hex3, Full HG
        }

        result = []
        for frag_name, loss_dict in GD1_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))

        return result

    @staticmethod
    def get_gt1a_hg_loss_fragments(precursor_formula):
        """
        Generate GT1a headgroup loss fragments.
        GT1a: Neu5Acα2-8Neu5Acα2-3Galβ1-3GalNAcβ1-4(Neu5Acα2-3)Galβ1-4Glcβ-Cer

        GT1a has 2 NeuAc (α2-8 linked) on terminal Gal and 1 NeuAc on internal Gal.
        """
        GT1A_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # Neu5Acα2-8NeuAc
            "HG(-Neu5Ac2,618)": {'C': 22, 'H': 38, 'N': 2, 'O': 18},  # 2 x Neu5Ac
            "HG(-Neu5Ac3,909)": {'C': 33, 'H': 55, 'N': 3, 'O': 26},  # Neu5Acα2-8NeuAc and Neu5Ac
            "HG(-Neu5Ac3,891)": {'C': 33, 'H': 53, 'N': 3, 'O': 25},  # Neu5Acα2-8NeuAc and Neu5Ac -H2O, YY fragment
            "HG(-Neu5Ac2HexNAcHex,965)": {'C': 36, 'H': 59, 'N': 3, 'O': 27},  # Neu5Ac2HexHexNAc
            "HG(-Neu5Ac3HexNAcHex,1274)": {'C': 47, 'H': 78, 'N': 4, 'O': 36},  # Neu5Ac2 and Neu5AcHexNAcHex
            "HG(-Neu5Ac3HexNAcHex,1256)": {'C': 47, 'H': 76, 'N': 4, 'O': 35},  # Neu5Ac2 and Neu5AcHexNAcHex -H2O, YY fragment
        }

        result = []
        for frag_name, loss_dict in GT1A_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result


    @staticmethod
    def get_gt1b_hg_loss_fragments(precursor_formula):
        """
        Generate GT1b headgroup loss fragments.
        GT1b: Neu5Acα2-3Galβ1-3GalNAcβ1-4(Neu5Acα2-8Neu5Acα2-3)Galβ1-4Glcβ-Cer

        GT1b has 1 NeuAc on terminal Gal and 2 NeuAc (α2-8 linked) on internal Gal.
        """
        GT1B_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # Neu5Acα2-8NeuAc
            "HG(-Neu5Ac2,618)": {'C': 22, 'H': 38, 'N': 2, 'O': 18},  # 2 x Neu5Ac
            "HG(-Neu5Ac3,909)": {'C': 33, 'H': 55, 'N': 3, 'O': 26},  # Neu5Acα2-8NeuAc and Neu5Ac
            "HG(-Neu5Ac3,891)": {'C': 33, 'H': 53, 'N': 3, 'O': 25},  # Neu5Acα2-8NeuAc and Neu5Ac -H2O, YY fragment
            "HG(-Neu5AcHexNAcHex,674)": {'C': 25, 'H': 42, 'N': 2, 'O': 19},  # Neu5AcHexNAcHex
            "HG(-Neu5Ac3HexNAcHex,1274)": {'C': 47, 'H': 78, 'N': 4, 'O': 36},  # Neu5Ac2 and Neu5AcHexNAcHex
            "HG(-Neu5Ac3HexNAcHex,1256)": {'C': 47, 'H': 76, 'N': 4, 'O': 35},  # Neu5Ac2 and Neu5AcHexNAcHex -H2O, YY fragment
        }

        result = []
        for frag_name, loss_dict in GT1B_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result


    @staticmethod
    def get_gt1c_hg_loss_fragments(precursor_formula):
        """
        Generate GT1c headgroup loss fragments.
        GT1c: Galβ1-3GalNAcβ1-4(NeuAcα2-8NeuAcα2-8NeuAcα2-3)Galβ1-4Glcβ-Cer

        """
        GT1C_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},
            "HG(-Neu5Ac3,891)": {'C': 33, 'H': 53, 'N': 3, 'O': 25},  # Neu5Acα2-8Neu5Acα2-8NeuAc
            "HG(-Hex,180)": {'C': 6, 'H': 12, 'O': 6},  # terminal Gal
            "HG(-HexHexNAc,383)": {'C': 14, 'H': 25, 'N': 1, 'O': 11},  # terminal Gal-GalNAc
            "HG(-HexNeu5Ac3,1071)": {'C': 39, 'H': 65, 'N': 3, 'O': 31}, # Neu5Ac3 and Gal
            "HG(-Neu5Ac3HexNAcHex,1274)": {'C': 47, 'H': 78, 'N': 4, 'O': 36},  # Neu5Ac3 and Gal-GalNAc
        }

        result = []
        for frag_name, loss_dict in GT1C_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_gt2_hg_loss_fragments(precursor_formula):
        """
        Generate GT2 headgroup loss fragments.
        GT2: GalNAcβ1-4(Neu5Acα2-8Neu5Acα2-8Neu5Acα2-3)Galβ1-4Glcβ-Cer

        """
        GT2_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # Neu5Acα2-8NeuAc
            "HG(-Neu5Ac3,891)": {'C': 33, 'H': 53, 'N': 3, 'O': 25},  # Neu5Acα2-8Neu5Acα2-8NeuAc
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},  # Terminal GalNAc
            "HG(-HexNAcNeu5Ac3,1094)": {'C': 41, 'H': 66, 'N': 4, 'O': 30}, # Neu5Ac3 and GalNAc
            "HG(-Neu5Ac3HexNAc,1112)": {'C': 41, 'H': 68, 'N': 4, 'O': 31},  # Neu5Ac3 and GalNAc and H2O
        }

        result = []
        for frag_name, loss_dict in GT2_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_gt3_hg_loss_fragments(precursor_formula):
        """
        Generate GT3 headgroup loss fragments.
        GT3: Neu5Acα2-8Neu5Acα2-8Neu5Acα2-3Galβ1-4Glcβ-Cer

        """
        GT3_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # Neu5Acα2-8NeuAc
            "HG(-Neu5Ac3,891)": {'C': 33, 'H': 53, 'N': 3, 'O': 25},  # Neu5Acα2-8Neu5Acα2-8NeuAc
            "HG(-Neu5Ac3,909)": {'C': 33, 'H': 55, 'N': 3, 'O': 26},  # # Neu5Acα2-8Neu5Acα2-8NeuAc and H2O
            "HG(-Neu5Ac3Hex,1053)": {'C': 39, 'H': 63, 'N': 3, 'O': 30},  # Neu5Ac3-Gal
            "HG(-Neu5Ac3Hex,1215)": {'C': 45, 'H': 73, 'N': 3, 'O': 35},  # entire HG
        }

        result = []
        for frag_name, loss_dict in GT3_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_gq1_hg_loss_fragments(precursor_formula):
        """
        Generate GQ1 headgroup loss fragments.
        GQ1: Neu5Acα2-8Neu5Acα2-8Neu5Acα2-3Galβ1-3GalNAcβ1-4(Neu5Acα2-3)Galβ1-4Glcβ-Cer

        """
        GQ1_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # Neu5Acα2-8NeuAc
            "HG(-Neu5Ac3,891)": {'C': 33, 'H': 53, 'N': 3, 'O': 25},  # Neu5Acα2-8Neu5Acα2-8NeuAc
            "HG(-Neu5Ac4,1182)": {'C': 44, 'H': 70, 'N': 4, 'O': 33},  # 2 x Neu5Acα2-8NeuAc -H2O, YY ion
            "HG(-Neu5Ac4HexNAcHex,1547)": {'C': 58, 'H': 93, 'N': 5, 'O': 43},  # Neu5Ac4 + HexNAc + Hex
            "HG(-Neu5Ac4HexNAcHex,1727)": {'C': 64, 'H': 105, 'N': 5, 'O': 49},  # Neu5Ac4HexNAcHex2, Z1 ion
            "HG(-Neu5Ac4HexNAcHex,1709)": {'C': 64, 'H': 103, 'N': 5, 'O': 48},  # Neu5Ac4HexNAcHex2, Y1 ion
            "HG(-Neu5Ac4HexNAcHex,1871)": {'C': 70, 'H': 113, 'N': 5, 'O': 53},  # Neu5Ac4HexNAcHex3, Z0 ion
            "HG(-Neu5Ac4HexNAcHex,1887)": {'C': 70, 'H': 113, 'N': 5, 'O': 54},  # Neu5Ac4HexNAcHex3, Y0 ion
        }

        result = []
        for frag_name, loss_dict in GQ1_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_gp1_hg_loss_fragments(precursor_formula):
        """
        Generate GP1 headgroup loss fragments.
        GP1: Neu5Acα2-8Neu5Acα2-8Neu5Acα2-8Neu5Acα2-3Galβ1-3GalNAcβ1-4(Neu5Acα2-3)Galβ1-4Glcβ-Cer
        """
        GP1_HEADGROUP_LOSSES = {
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9},
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # Neu5Acα2-8NeuAc
            "HG(-Neu5Ac3,891)": {'C': 33, 'H': 53, 'N': 3, 'O': 25},  # Neu5Acα2-8Neu5Acα2-8NeuAc
            "HG(-Neu5Ac4,1182)": {'C': 44, 'H': 70, 'N': 4, 'O': 33},  # Neu5Ac4 chain
            "HG(-Neu5Ac4,1200)": {'C': 44, 'H': 72, 'N': 4, 'O': 34},  # Neu5Ac4 + H2O
            "HG(-Neu5Ac5,1473)": {'C': 55, 'H': 87, 'N': 5, 'O': 41},  # Neu5Ac5
            "HG(-Neu5Ac5,1491)": {'C': 55, 'H': 89, 'N': 5, 'O': 42},  # Neu5Ac5 + H2O
            "HG(-Neu5Ac4HexNAcHex,1547)": {'C': 58, 'H': 93, 'N': 5, 'O': 43},  # Neu5Ac4 + HexNAc + Hex
            "HG(-Neu5Ac5HexNAcHex,1838)": {'C': 69, 'H': 110, 'N': 6, 'O': 51},  # Neu5Ac5 + HexNAc + Hex
        }

        result = []
        for frag_name, loss_dict in GP1_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_nlc10_hg_loss_fragments(precursor_formula):
        """
        Generate nLc10 headgroup loss fragments.
        nLc10: GlcNAcβ1-3Galβ1-4GlcNAcβ1-3(Galα1-3Galβ1-4GlcNAcβ1-6)Galβ1-4GlcNAcβ1-3Galβ1-4Glcβ-Cer
        Branched structure with terminal GlcNAc on both antennae
        """
        NLC10_HEADGROUP_LOSSES = {
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},
            "HG(-HexNAcHex,383)": {'C': 14, 'H': 25, 'N': 1, 'O': 11},      # HexNAc + Hex
            "HG(-HexNAc2Hex,586)": {'C': 22, 'H': 38, 'N': 2, 'O': 16},
            "HG(-HexNAc4Hex4,1478)": {'C': 56, 'H': 94, 'N': 4, 'O': 44}  # HexNAc4Hex4
        }

        result = []
        for frag_name, loss_dict in NLC10_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_nlc8_hg_loss_fragments(precursor_formula):
        """
        Generate nLc8 headgroup loss fragments.
        nLc8: Galβ1-4GlcNAcβ1-3(Galβ1-4GlcNAcβ1-6)Galβ1-4GlcNAcβ1-3Galβ1-4Glcβ-Cer
        Branched structure (I antigen precursor)
        """
        NLC8_HEADGROUP_LOSSES = {
            "HG(-Hex,180)": {'C': 6, 'H': 12, 'O': 6},
            "HG(-HexHexNAc,365)": {'C': 14, 'H': 23, 'N': 1, 'O': 10},
            "HG(-HexHexNAc,383)": {'C': 14, 'H': 25, 'N': 1, 'O': 11},
            "HG(-Hex3HexNAc3,1113)": {'C': 42, 'H': 71, 'N': 3, 'O': 31},  # Y2
            "HG(-Hex3HexNAc3,1257)": {'C': 42, 'H': 71, 'N': 3, 'O': 32},  # Z2
        }

        result = []
        for frag_name, loss_dict in NLC8_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result

    @staticmethod
    def get_nlc6_hg_loss_fragments(precursor_formula):
        """
        Generate nLc6 headgroup loss fragments.
        nLc6: Galβ1-4GlcNAcβ1-3Galβ1-4GlcNAcβ1-3Galβ1-4Glcβ-Cer
        Linear structure (i antigen)
        """
        NLC6_HEADGROUP_LOSSES = {
            "HG(-Hex,180)": {'C': 6, 'H': 12, 'O': 6},
            "HG(-HexHexNAc,365)": {'C': 14, 'H': 23, 'N': 1, 'O': 10},
            "HG(-HexHexNAc,383)": {'C': 14, 'H': 25, 'N': 1, 'O': 11},
            "HG(-Hex2HexNAc2,748)": {'C': 28, 'H': 48, 'N': 2, 'O': 21},
            "HG(-Hex2HexNAc2,766)": {'C': 28, 'H': 50, 'N': 2, 'O': 22},
            "HG(-Hex3HexNAc2,910)": {'C': 34, 'H': 58, 'N': 2, 'O': 26},
            "HG(-Hex3HexNAc2,928)": {'C': 34, 'H': 60, 'N': 2, 'O': 27},
        }

        result = []
        for frag_name, loss_dict in NLC6_HEADGROUP_LOSSES.items():
            frag_formula = MolecularFormula(dict(precursor_formula.elements))
            for elem, loss_val in loss_dict.items():
                frag_formula.elements[elem] = frag_formula.elements.get(elem, 0) - loss_val
            frag_formula = MolecularFormula({k: v for k, v in frag_formula.elements.items() if v > 0})
            result.append((frag_name, str(frag_formula), frag_formula.mass()))
        return result


    @staticmethod # positive mode dispatcher, aglycone fragments and head group oxonium ions, a-, b-, c- type fragments
    def get_headgroup_fragments(lipid_class, precursor_formula):
        if lipid_class == 'Hex':
            return (
            GSLFragmentRules.get_hex_hg_loss_fragments(precursor_formula)
            + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'Lac':
            return (
            GSLFragmentRules.get_lac_hg_loss_fragments(precursor_formula)
            + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'Gb3':
            return (
                GSLFragmentRules.get_gb3_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'Gb4':
            return (
                GSLFragmentRules.get_gb4_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GA1':
            return (
                GSLFragmentRules.get_ga1_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GA2':
            return (
                GSLFragmentRules.get_ga2_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'LC3':
            return (
                GSLFragmentRules.get_lc3_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'LC4':
            return (
                GSLFragmentRules.get_lc4_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'SM4':
            return (
                GSLFragmentRules.get_sm4_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'SHex2':
            return (
                GSLFragmentRules.get_shex2_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GM4':
            return (GSLFragmentRules.get_gm4_hg_loss_fragments(precursor_formula) +
                    GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GM3':
            return (GSLFragmentRules.get_gm3_hg_loss_fragments(precursor_formula) +
                    GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GM2':
            return (GSLFragmentRules.get_gm2_hg_loss_fragments(precursor_formula) +
                    GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GM1':
            return (GSLFragmentRules.get_gm1_hg_loss_fragments(precursor_formula) +
                    GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GD3':
            return (GSLFragmentRules.get_gd3_hg_loss_fragments(precursor_formula) +
                    GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GD2':  # ← ADD THIS
            return (GSLFragmentRules.get_gd2_hg_loss_fragments(precursor_formula) +
                    GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class in ['GD1a', 'GD1b']:
            is_a = (lipid_class == 'GD1a')
            return (GSLFragmentRules.get_gd1_hg_loss_fragments(precursor_formula, is_a=is_a)
                    + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GT3':
            return (
                GSLFragmentRules.get_gt3_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
            )
        elif lipid_class == 'GT2':
            return (
                GSLFragmentRules.get_gt2_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
            )

        elif lipid_class == 'GT1a':
            return (
                GSLFragmentRules.get_gt1a_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )

        elif lipid_class == 'GT1b':
            return (
                GSLFragmentRules.get_gt1b_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )

        elif lipid_class == 'GT1c':
            return (
                GSLFragmentRules.get_gt1c_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GQ1':
            return (
                GSLFragmentRules.get_gq1_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'GP1':
            return (
                GSLFragmentRules.get_gp1_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
        )
        elif lipid_class == 'nLc10':
            return (
                GSLFragmentRules.get_nlc10_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
            )
        elif lipid_class == 'nLc8':
            return (
                GSLFragmentRules.get_nlc8_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
            )
        elif lipid_class == 'nLc6':
            return (
                GSLFragmentRules.get_nlc6_hg_loss_fragments(precursor_formula)
                + GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)
            )

        else:
            return GSLFragmentRules.get_headgroup_fragments_positive(lipid_class)

    @staticmethod
    def get_headgroup_fragments_positive(gsl_class: str):
        """
        Return common protonated headgroup a- ,b- ,c-fragments (+ESI mode).

        """
        fragments = []

        if gsl_class in ['Lac', 'LC3', 'LC4']:
            fragments.extend([
                ("HG(Hex,162)", "C6H10O5", 163.0606, True),
                ("HG(Hex,180)", "C6H12O6", 181.0707, True),
                ("HG(Hex2,324)", "C12H20O10", 325.1135, True),
                ("HG(Hex2,342)", "C12H22O11", 343.1235, True),
            ])
            if gsl_class in ['LC3', 'LC4']:
                fragments.extend([
                    ("HG(HexNAc,221)", "C8H15NO6", 222.0972, True),        # N-acetylhexosamine (GalNAc)
                    ("HG(HexNAc,203)", "C8H13NO5", 204.0866, True),        # HexNAc -H2O
                    ("HG(HexNAc,185)", "C8H11NO4", 186.0761, True),        # HexNAc -2H2O
                    ("HG(HexNAc,155)", "C7H9NO3", 156.0655, True),         # HexNAc -H2O -CH4O2
                    ("HG(HexNAc,137)", "C7H7NO2", 138.0550, True),         # HexNAc -3H2O -CH4O2
                ])

        elif gsl_class == 'Hex':
            fragments.extend([
                ("HG(Hex,162)", "C6H10O5", 163.0606, True),
                ("HG(Hex,180)", "C6H12O6", 181.0707, True),
            ])

        elif gsl_class.startswith('GA'):
            fragments.extend([
                ("HG(Hex,180)", "C6H12O6", 181.0707, True),
                ("HG(HexNAc,221)", "C8H15NO6", 222.0972, True),         # N-acetylhexosamine (GalNAc)
                ("HG(HexNAc,203)", "C8H13NO5", 204.0866, True),         # HexNAc -H2O
                ("HG(HexNAc,185)", "C8H11NO4", 186.0761, True),         # HexNAc -2H2O
                ("HG(HexNAc,155)", "C7H9NO3", 156.0655, True),         # HexNAc -H2O -CH4O2
                ("HG(HexNAc,137)", "C7H7NO2", 138.0550, True),         # HexNAc -3H2O -CH4O2
                ("HG(HexHexNAc,383)", "C14H25NO11", 384.1505, True),
            ])

        elif gsl_class in ['Gb3', 'Gb4']:
            fragments.extend([
                ("HG(Hex,162)", "C6H10O5", 163.0606, True),
                ("HG(Hex,180)", "C6H12O6", 181.0707, True),
                ("HG(Hex2,324)", "C12H20O10", 325.1135, True),
                ("HG(Hex3,487)", "C18H30O15", 487.1664, True),
            ])
            if gsl_class == 'Gb4':
                fragments.extend([
                    ("HG(HexNAc,221)", "C8H15NO6", 222.0972, True),         # HexNAc
                    ("HG(HexNAc,203)", "C8H13NO5", 204.0866, True),         # HexNAc -H2O
                    ("HG(HexNAc,185)", "C8H11NO4", 186.0761, True),         # HexNAc -2H2O
                    ("HG(HexNAc,155)", "C7H9NO3", 156.0655, True),          # HexNAc -H2O -CH4O2
                    ("HG(HexNAc,137)", "C7H7NO2", 138.0550, True),          # HexNAc -3H2O -CH4O2
                    ("HG(HexNAcHex,365)", "C14H23NO10", 366.1400, True),
                ])

        elif gsl_class == 'SM':
            fragments.extend([
                ("Phosphocholine", "C5H15NO4P", 184.0739, True),
                ("Phosphocholine-H2O", "C5H13NO3P", 166.0628, True),
            ])

        elif gsl_class == 'SM4':
            fragments.extend([
                ("HG(Hex,180)", "C6H12O6", 181.0707, True),
            ])

        elif gsl_class.startswith('GM'):
            # Common to all sialic acid-containing GSLs
            fragments.extend([
                ('HG(NeuAc,309)', 'C11H19NO9', 310.1133, True),  # Protonated sialic acid
                ('HG(NeuAc,291)', 'C11H17NO8', 292.1027, True),  # Protonated sialic acid -H2O
            ])

        # GM4, NeuAc + Gal
        if gsl_class == 'GM4':
            fragments.extend([
                ('HG(NeuAcGal,471)', 'C17H29NO14', 472.1661, True),  # Full headgroup
                ('HG(NeuAcGal,453)', 'C17H27NO13', 454.1555, True),  # Full headgroup -H2O
            ])
        elif gsl_class == 'GM3':
            fragments.extend([
                ('HG(Hex2,342)', 'C12H22O11', 343.1235, True),
            ])
        elif gsl_class in ['GM2', 'GM1']:
            fragments.extend([
                ("HG(HexNAc,221)", "C8H15NO6", 222.0972, True),         # HexNAc
                ("HG(HexNAc,203)", "C8H13NO5", 204.0866, True),         # HexNAc -H2O
                ("HG(HexNAc,185)", "C8H11NO4", 186.0761, True),         # HexNAc -2H2O
                ("HG(HexNAc,155)", "C7H9NO3", 156.0655, True),         # HexNAc -H2O -CH4O2
                ("HG(HexNAc,137)", "C7H7NO2", 138.0550, True),         # HexNAc -3H2O -CH4O2
                ("HG(HexNAcHex,383)", "C14H25NO11", 384.1500, True),     # HexNAc + Hex
                ("HG(HexNAcHex,365)", "C14H23NO10", 366.1395, True),     # HexNAc + Hex -H2O
                ("HG(HexNAcHex2,545)", "C20H35NO16", 546.2029, True),    # HexNAc + 2x Hex
                ("HG(HexNAcHex2,527)", "C20H33NO15", 528.1923, True),    # HexNAc + 2x Hex -H20
            ])
        elif gsl_class.startswith('GD'):
            fragments.extend([
                ('HG(NeuAc,309)', 'C11H19NO9', 310.1133, True),  # Protonated sialic acid
                ('HG(NeuAc,291)', 'C11H17NO8', 292.1027, True),  # Protonated sialic acid -H2O
                ("HG(NeuAc2,600)", "C22H36N2O17", 601.2087, True),  # Protonated NeuAc2
                ("HG(NeuAc2,582)", "C22H34N2O16", 583.1981, True),  # Protonated NeuAc2 -H2O
            ])
            if gsl_class == 'GD2':
                fragments.extend([
                    ("HG(HexNAc,221)", "C8H15NO6", 222.0972, True),         # HexNAc
                    ("HG(HexNAc,203)", "C8H13NO5", 204.0866, True),         # HexNAc -H2O
                    ("HG(HexNAc,185)", "C8H11NO4", 186.0761, True),         # HexNAc -2H2O
                    ("HG(HexNAc,155)", "C7H9NO3", 156.0655, True),         # HexNAc -H2O -CH4O2
                    ("HG(HexNAc,137)", "C7H7NO2", 138.0550, True),         # HexNAc -3H2O -CH4O2
                    ("HG(HexNAcHex,383)", "C14H25NO11", 384.1500, True),     # HexNAc + Hex
                    ("HG(HexNAcHex,365)", "C14H23NO10", 366.1395, True),     # HexNAc + Hex -H2O
                    ("HG(NeuAc2Hex,744)", "C28H44N2O21", 745.2509, True),     # NeuAc2 + Hex -H2O, double cleaved fragment B3Y2β
                    ("HG(Neu5Ac2HexNAcHex,947)", "C36H57N3O26", 948.3303, True),     # NeuAc2 + HexNAc + Hex -H2O
                    ("HG(Neu5Ac2HexNAcHex2,1109)", "C14H25NO11", 1110.3831, True),     # Full HG, NeuAc2 + HexNAc + Hex2 -H2O
                ])

        elif gsl_class.startswith('GT'):
            fragments.extend([
                ('HG(NeuAc,309)', 'C11H19NO9', 310.1133, True),  # Protonated sialic acid
                ('HG(NeuAc,291)', 'C11H17NO8', 292.1027, True),  # Protonated sialic acid -H2O
                ("HG(NeuAc2,600)", "C22H36N2O17", 601.2087, True),  # Protonated NeuAc2
                ("HG(NeuAc2,582)", "C22H34N2O16", 583.1981, True),  # Protonated NeuAc2 -H2O
                ("HG(HexNAcHex,365)", "C14H23NO10", 366.1395, True),
            ])
            if gsl_class == 'GT3':
                fragments.extend([
                    ("HG(Hex,180)", "C6H12O6", 181.0707, True),  # Single Gal
                ])
            if gsl_class == 'GT2':
                fragments.extend([
                    ("HG(HexNAc,221)", "C8H15NO6", 222.0972, True),         # HexNAc
                    ("HG(HexNAc,203)", "C8H13NO5", 204.0866, True),         # HexNAc -H2O
                    ("HG(HexNAc,185)", "C8H11NO4", 186.0761, True),         # HexNAc -2H2O
                    ("HG(HexNAc,155)", "C7H9NO3", 156.0655, True),         # HexNAc -H2O -CH4O2
                    ("HG(HexNAc,137)", "C7H7NO2", 138.0550, True),         # HexNAc -3H2O -CH4O2
                ])

        elif gsl_class in ['GQ1', 'GP1']:
            fragments.extend([
                ('HG(NeuAc,309)', 'C11H19NO9', 310.1133, True),  # Protonated sialic acid
                ('HG(NeuAc,291)', 'C11H17NO8', 292.1027, True),  # Protonated sialic acid -H2O
                ("HG(NeuAc2,600)", "C22H36N2O17", 601.2087, True),  # Protonated NeuAc2
                ("HG(NeuAc2,582)", "C22H34N2O16", 583.1981, True),  # Protonated NeuAc2 -H2O
                ("HG(HexNAcHex,365)", "C14H23NO10", 366.1395, True),
            ])
            if gsl_class == 'GP1':
                fragments.append(("HG(NeuAc4HexNAcHex,1529)", "C58H91N5O42", 1530.5211, True))

        elif gsl_class in ['nLc10', 'nLc8', 'nLc6']:
            fragments.extend([
                ("HG(Hex,180)", "C6H12O6", 181.0707, True),
                ("HG(HexNAc,221)", "C8H15NO6", 222.0972, True),         # HexNAc
                ("HG(HexNAc,203)", "C8H13NO5", 204.0866, True),         # HexNAc -H2O
                ("HG(HexNAc,185)", "C8H11NO4", 186.0761, True),         # HexNAc -2H2O
                ("HG(HexNAc,155)", "C7H9NO3", 156.0655, True),         # HexNAc -H2O -CH4O2
                ("HG(HexNAc,137)", "C7H7NO2", 138.0550, True),         # HexNAc -3H2O -CH4O2
                ("HG(HexNAcHex,383)", "C14H25NO11", 384.1500, True),
                ("HG(HexNAcHex,365)", "C14H23NO10", 366.1395, True),
            ])

        elif gsl_class.startswith('SHex2'):
            fragments.extend([
                ("HG(Hex,180)", "C6H12O6", 181.0707, True),
                ("HG(Hex2,342)", "C12H22O11", 343.1235, True),
            ])

        return fragments

    @staticmethod
    def get_headgroup_fragments_negative(gsl_class: str, precursor_formula=None):
        """
        Return negative-ion mode headgroup fragments.
        Pre-deprotonated [M-H]- masses.
        """
        fragments = []

        if gsl_class == 'SM4':
            fragments.extend([
                ("HG(HSO4,97)", "H2SO4", 96.9596, True),          # Bisulfate [M-H]-
                ("HG(SHexCer,242)", "C6H10O8S", 241.0017, True),  # Sulfohexose [M-H]-
                ("HG(SHexCer)+(C2H5NO)", "C8H15NO9S", 300.0389, True),  # +C2H5NO
            ])
            # Add headgroup loss fragments if a precursor formula is provided
            if precursor_formula:
                fragments.extend(GSLFragmentRules.get_sm4_hg_loss_fragments(precursor_formula))


        elif gsl_class == 'SHex2':
            fragments.extend([
                ("HG(HSO4,97)", "H2SO4", 96.9596, True),          # Bisulfate [M-H]-
                ("HG(SO3,80)", "HSO3", 79.9568, True),            # Sulfate radical anion
                ("HG(SHex,260)", "C6H12SO9", 259.0124, True),    # Sulfated hexose [M-H]-
                ("HG(SHex,242)", "C6H10SO8", 241.0018, True),    # Dehydrated
                ("HG(SHexHex,404)", "C12H20SO13", 403.0546, True), # Sulfated lactose [M-H]-
                ("HG(SHexHex,386)", "C12H18SO12", 385.0441, True), # -H2O
                ("HG(Hex,180)", "C6H12O6", 179.0556, True),      # Hexose [M-H]-
                ("HG(Hex2,342)", "C12H22O11", 341.1084, True),   # Lactose [M-H]-
            ])
            if precursor_formula:
                fragments.extend(GSLFragmentRules.get_shex2_hg_loss_fragments(precursor_formula))

        # Handle all sialic acid-containing GSLs
        elif gsl_class in ['GM4', 'GM3', 'GM2', 'GM1', 'GD3', 'GD2', 'GD1a', 'GD1b', 'GD1c', 'GT3', 'GT2', 'GT1a', 'GT1b', 'GT1c', 'GQ1', 'GP1']:
            # 1. Common diagnostic fragment for ALL sialic GSLs
            fragments.extend([
                ('HG(NeuAc,291)', 'C11H17NO8', 290.0876, True),  # Sialic acid -H2O
            ])

            # 2. Motif-specific diagnostic fragments (B-ions)
            # Disialo B-ions for GSLs with the Neu5Ac-Neu5Ac motif
            if gsl_class in ['GD3', 'GD2', 'GD1b']:
                fragments.extend([
                    ('HG(NeuAc2,582)', 'C22H34N2O16', 581.1836, True),  # Neu5Ac2 -H2O
                    ('HG(NeuAc2-CO2,538)', 'C21H34N2O14', 537.1937, True),  # Neu5Ac2 -H2O -CO2
                ])

            # Isomer-specific fragments for GD1a and GD1b
            if gsl_class in ['GD1a', 'GD1b']:
                # Shared fragment for both isomers
                fragments.extend([
                    ('HG(HexNAcHex,365)', 'C14H25NO10', 366.1395, True),
                ])
                # Unique fragment for GD1a
                if gsl_class == "GD1a":
                    fragments.extend([
                        ('HG(NeuAcHexNAcHex,656)', 'C25H40N2O18', 655.2203, True),
                    ])

            # Trisialo B-ions for GT1 series (3 sialic acids)
            if gsl_class in ['GT1a', 'GT1b', 'GT1c']:
                fragments.extend([
                    ('HG(NeuAc2,582)', 'C22H34N2O16', 581.1836, True),  # Neu5Ac2 -H2O
                    ('HG(NeuAc2-CO2,538)', 'C21H34N2O14', 537.1937, True),  # Neu5Ac2 -H2O -CO2
                ])
                # GT1b-specific: unsialylated terminal Gal
                if gsl_class == 'GT1b':
                    fragments.extend([
                        ('HG(NeuAcHexHexNAc,656)', 'C25H40N2O18', 655.2203, True),  # Neu5Ac-Gal-GalNAc B3α fragment, diagnostic
                    ])
                # GT1c-specific: unsialylated terminal Gal
                if gsl_class == 'GT1c':
                    fragments.extend([
                        ('HG(HexHexNAc,365)', 'C14H23NO10', 364.1249, True),
                        ('HG(NeuAc3,873)', 'C33H51N3O24', 872.2790, True),  # Neu5Ac3 -H2O
                        ('HG(NeuAc3-CO2,829)', 'C32H51N3O22', 828.2891, True),  # Neu5Ac3 -H2O -CO2
                    ])
            if gsl_class == 'GT2':
                fragments.extend([
                    ('HG(NeuAc2,582)', 'C22H34N2O16', 581.1836, True),
                    ('HG(NeuAc2-CO2,538)', 'C21H34N2O14', 537.1937, True),
                    ('HG(NeuAc3,873)', 'C33H51N3O24', 872.2790, True),      # Trisialo
                    ('HG(NeuAc3-CO2,829)', 'C32H51N3O22', 828.2891, True),  # Trisialo -CO2
                    ('HG(HexNAc,203)', 'C8H13NO5', 203.0721, True),         # Terminal GalNAc -H2O
                ])
            if gsl_class == 'GT3':
                fragments.extend([
                    ('HG(NeuAc2,582)', 'C22H34N2O16', 581.1836, True),
                    ('HG(NeuAc2-CO2,538)', 'C21H34N2O14', 537.1937, True),
                    ('HG(NeuAc3,873)', 'C33H51N3O24', 872.2790, True),
                    ('HG(NeuAc3-CO2,829)', 'C32H51N3O22', 828.2891, True),
                ])
            if gsl_class == 'GQ1':
                fragments.extend([
                    ('HG(NeuAc2,582)', 'C22H34N2O16', 581.1836, True),
                    ('HG(NeuAc2-CO2,538)', 'C21H34N2O14', 537.1937, True),
                    ('HG(NeuAc2Hex,762)', 'C28H46N2O22', 761.2469, True),
                    ('HG(NeuAc3,873)', 'C33H51N3O24', 872.2790, True),
                    ('HG(NeuAc3-CO2,829)', 'C32H51N3O22', 828.2891, True),
                ])
            if gsl_class == 'GP1':
                fragments.extend([
                    ('HG(NeuAc2,582)', 'C22H34N2O16', 581.1836, True),
                    ('HG(NeuAc2-CO2,538)', 'C21H34N2O14', 537.1937, True),
                    ('HG(NeuAc3,873)', 'C33H51N3O24', 872.2790, True),
                    ('HG(NeuAc3-CO2,829)', 'C32H51N3O22', 828.2891, True),
                    ('HG(NeuAc4,1164)', 'C44H68N4O32', 1163.3744, True),
                    ('HG(NeuAc4-CO2,1120)', 'C43H68N4O30', 1119.3846, True),
                ])

            # 3. Headgroup loss fragments (Y-ions) for each class
            if precursor_formula is not None:
                if gsl_class == "GM4":
                    fragments.extend(GSLFragmentRules.get_gm4_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GM3":
                    fragments.extend(GSLFragmentRules.get_gm3_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GM2":
                    fragments.extend(GSLFragmentRules.get_gm2_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GM1":
                    fragments.extend(GSLFragmentRules.get_gm1_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GD3":
                    fragments.extend(GSLFragmentRules.get_gd3_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GD2":
                    fragments.extend(GSLFragmentRules.get_gd2_hg_loss_fragments(precursor_formula))
                elif gsl_class in ["GD1a", "GD1b"]:
                    fragments.extend(GSLFragmentRules.get_gd1_hg_loss_fragments(precursor_formula, is_a=(gsl_class == "GD1a")))
                elif gsl_class == "GT1a":
                    fragments.extend(GSLFragmentRules.get_gt1a_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GT1b":
                    fragments.extend(GSLFragmentRules.get_gt1b_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GT1c":
                    fragments.extend(GSLFragmentRules.get_gt1c_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GT2":
                    fragments.extend(GSLFragmentRules.get_gt2_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GT3":
                    fragments.extend(GSLFragmentRules.get_gt3_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GQ1":
                    fragments.extend(GSLFragmentRules.get_gq1_hg_loss_fragments(precursor_formula))
                elif gsl_class == "GP1":
                    fragments.extend(GSLFragmentRules.get_gp1_hg_loss_fragments(precursor_formula))

        elif gsl_class in ['nLc10', 'nLc8', 'nLc6']:
            fragments.extend([
                ('HG(Hex,180)', 'C6H12O6', 179.0561, True),
                ('HG(Hex,162)', 'C6H10O5', 161.0455, True),
                ('HG(HexNAc,221)', 'C8H15NO6', 220.0827, True),
                ('HG(HexNAc,203)', 'C8H13NO5', 202.0721, True),
                ('HG(HexNAcHex,383)', 'C14H25NO11', 382.1355, True),
                ('HG(HexNAcHex,365)', 'C14H23NO10', 364.1249, True),
            ])
            if precursor_formula:
                if gsl_class == 'nLc10':
                    fragments.extend(GSLFragmentRules.get_nlc10_hg_loss_fragments(precursor_formula))
                elif gsl_class == 'nLc8':
                    fragments.extend(GSLFragmentRules.get_nlc8_hg_loss_fragments(precursor_formula))
                elif gsl_class == 'nLc6':
                    fragments.extend(GSLFragmentRules.get_nlc6_hg_loss_fragments(precursor_formula))

        return fragments


# ============================================================================
# NEGATIVE ION MODE FRAGMENT RULES
# ============================================================================

class NegativeFragmentRules:
    """
    Fragment rules for negative ion mode (universal for all GSLs).
    All masses are pre-deprotonated [M-H]- values.
    """
    @staticmethod
    def get_lcb_fragments_negative(lcb_type: str):
        """
        Generate LCB fragments for negative ion mode.
        Covers all LCB types used in positive mode with negative-mode specific fragments.
        Returns list of (name, formula_str, [M-H]- m/z, is_pre_deprotonated)

        Common negative-mode LCB fragments:
        - LCB(-CH3O): Loss of methanol group
        - LCB(-C2H8NO): Loss of ethanolamine-like fragment
        """
        fragments = []

        # ========== Standard dihydroxy LCBs (;2) ==========

        # d16:0;2
        if lcb_type == "16:0;2":
            fragments.extend([
                ("LCB(-CH3O)", "C15H32NO", 241.2406, True),
                ("LCB(-C2H8NO)", "C14H28O", 211.2062, True),
            ])

        # d16:1;2
        elif lcb_type == "16:1;2":
            fragments.extend([
                ("LCB(-CH3O)", "C15H30NO", 239.2249, True),
                ("LCB(-C2H8NO)", "C14H26O", 209.1905, True),
            ])

        # d17:0;2
        elif lcb_type == "17:0;2":
            fragments.extend([
                ("LCB(-CH3O)", "C16H34NO", 255.2562, True),
                ("LCB(-C2H8NO)", "C15H30O", 225.2218, True),
            ])

        # d17:1;2
        elif lcb_type == "17:1;2":
            fragments.extend([
                ("LCB(-CH3O)", "C16H32NO", 253.2406, True),
                ("LCB(-C2H8NO)", "C15H28O", 223.2062, True),
            ])

        # d18:0;2 (sphinganine)
        elif lcb_type == "18:0;2":
            fragments.extend([
                ("LCB(-CH3O)", "C17H37NO", 270.2797, True),
                ("LCB(-C2H8NO)", "C16H32O", 239.2375, True),
            ])

        # d18:1;2 (sphingosine) - Most common
        elif lcb_type == "18:1;2":
            fragments.extend([
                ("LCB(-CH3O)", "C17H35NO", 268.2640, True),
                ("LCB(-C2H8NO)", "C16H30O", 237.2218, True),
            ])

        # d18:2;2 (sphingadiene)
        elif lcb_type == "18:2;2":
            fragments.extend([
                ("LCB(-CH3O)", "C17H33NO", 266.2484, True),
                ("LCB(-C2H8NO)", "C16H28O", 235.2062, True),
            ])

        # d18:0;3 (phytosphingosine - trihydroxy)
        elif lcb_type == "18:0;3":
            fragments.extend([
                ("LCB(-CH3O)", "C17H36NO2", 285.2668, True),
                ("LCB(-C2H8NO)", "C16H32O2", 255.2324, True),
            ])

        # d19:0;2
        elif lcb_type == "19:0;2":
            fragments.extend([
                ("LCB(-CH3O)", "C18H38NO", 283.2875, True),
                ("LCB(-C2H8NO)", "C17H34O", 253.2531, True),
            ])

        # d19:1;2
        elif lcb_type == "19:1;2":
            fragments.extend([
                ("LCB(-CH3O)", "C18H36NO", 281.2719, True),
                ("LCB(-C2H8NO)", "C17H32O", 251.2375, True),
            ])

        # d20:0;2
        elif lcb_type == "20:0;2":
            fragments.extend([
                ("LCB(-CH3O)", "C19H40NO", 297.3032, True),
                ("LCB(-C2H8NO)", "C18H36O", 267.2688, True),
            ])

        # d20:1;2
        elif lcb_type == "20:1;2":
            fragments.extend([
                ("LCB(-CH3O)", "C19H38NO", 295.2875, True),
                ("LCB(-C2H8NO)", "C18H34O", 265.2531, True),
            ])

        # ========== 1-Deoxy LCBs (;1) - if needed for negative mode ==========
        # Note: 1-deoxy sphingolipids are less common in negative mode
        # but included for completeness

        # d18:0;1 (1-deoxysphingosine)
        elif lcb_type == "18:0;1":
            fragments.extend([
                ("doxLCB(-CH3O)", "C17H36N", 253.2769, True),
                ("doxLCB(-C2H8N)", "C16H32", 223.2426, True),
            ])

        # d18:1;1
        elif lcb_type == "18:1;1":
            fragments.extend([
                ("doxLCB(-CH3O)", "C17H34N", 251.2613, True),
                ("doxLCB(-C2H8N)", "C16H30", 221.2269, True),
            ])

        # d19:0;1
        elif lcb_type == "19:0;1":
            fragments.extend([
                ("doxLCB(-CH3O)", "C18H38N", 267.2926, True),
                ("doxLCB(-C2H8N)", "C17H34", 237.2582, True),
            ])

        # d20:0;1
        elif lcb_type == "20:0;1":
            fragments.extend([
                ("doxLCB(-CH3O)", "C19H40N", 281.3082, True),
                ("doxLCB(-C2H8N)", "C18H36", 251.2739, True),
            ])

        return fragments

    @staticmethod
    def get_fa_fragments_negative(fa_type: str):
        """
        Generate fatty acid fragments for negative ion mode from ceramide-based lipids.

        In sphingolipids, FA is amide-linked: LCB-NH-CO-FA

        """
        fragments = []

        import re
        match = re.match(r'(\d+):(\d+)', fa_type)
        if not match:
            return fragments

        carbon = int(match.group(1))
        double_bonds = int(match.group(2))
        h_count = 2 * carbon - 2 * double_bonds

        # =====================================================================
        # 1. FA+(HN) - Retains amide fragment with carbonyl
        # =====================================================================
        # Fragment: R-CO-NH with only the carbonyl oxygen
        # Formula: CnH(2n-2db+1)NO (ONE oxygen from C=O)
        fa_hn_formula = f"C{carbon}H{h_count+1}NO"
        fa_hn_mass = (
            carbon * ATOMIC_MASSES['C'] +
            (h_count + 1) * ATOMIC_MASSES['H'] +
            1 * ATOMIC_MASSES['N'] +
            1 * ATOMIC_MASSES['O'] -  # Only 1 O (carbonyl)
            ATOMIC_MASSES['H']  # [M-H]-
        )

        # =====================================================================
        # 2. FA+(C2H3N) - Acetonitrile adduct with carbonyl
        # =====================================================================
        # Adduct: R-CO-NH + CH3CN
        # Formula: C(n+2)H(2n-2db+3)NO (ONE oxygen from C=O)
        fa_acn_formula = f"C{carbon+2}H{h_count+3}NO"
        fa_acn_mass = (
            (carbon + 2) * ATOMIC_MASSES['C'] +
            (h_count + 3) * ATOMIC_MASSES['H'] +
            1 * ATOMIC_MASSES['N'] +
            1 * ATOMIC_MASSES['O'] -  # Only 1 O (carbonyl)
            ATOMIC_MASSES['H']  # [M-H]-
        )

        # =====================================================================
        # 3. FA+(C2H3NO) - Acetamide-like adduct
        # =====================================================================
        # Adduct with formyl/acetyl group: adds C2H3NO
        # Formula: C(n+2)H(2n-2db+3)NO2 (TWO oxygens: original + adduct)
        fa_c2h3no_formula = f"C{carbon+2}H{h_count+3}NO2"
        fa_c2h3no_mass = (
            (carbon + 2) * ATOMIC_MASSES['C'] +
            (h_count + 3) * ATOMIC_MASSES['H'] +
            1 * ATOMIC_MASSES['N'] +
            2 * ATOMIC_MASSES['O'] -  # 2 O (carbonyl + adduct)
            ATOMIC_MASSES['H']  # [M-H]-
        )

        fragments.extend([
            (f"FA {fa_type}+(HN)", fa_hn_formula, fa_hn_mass, True),
            (f"FA {fa_type}+(C2H3N)", fa_acn_formula, fa_acn_mass, True),
            (f"FA {fa_type}+(C2H3NO)", fa_c2h3no_formula, fa_c2h3no_mass, True),
        ])

        return fragments




class CeramideFragmentRules:
    """Fragment rules for ceramides, positive ion mode"""

    @staticmethod
    def get_ceramide_fragments(lcb_type: str, is_doxcer: bool = False):
        fragments = []
        if is_doxcer:
            if lcb_type == "18:0;1":
                fragments.extend([
                    ("doxLCB 18:0;1", "C18H39NO", 286.3104),
                    ("doxLCB 18:0;1(-H2O)", "C18H37N", 268.2999),
                ])
            elif lcb_type == "18:1;1":
                fragments.extend([
                    ("doxLCB 18:1;1", "C18H37NO", 284.2948),
                    ("doxLCB 18:1;1(-H2O)", "C18H35N", 266.2842),
                ])
            elif lcb_type == "19:0;1":
                fragments.extend([
                    ("doxLCB 19:0;1", "C19H41NO", 300.3261),
                    ("doxLCB 19:0;1(-H2O)", "C19H39N", 282.3156),
                ])
            elif lcb_type == "20:0;1":
                fragments.extend([
                    ("doxLCB 20:0;1", "C20H43NO", 314.3417),
                    ("doxLCB 20:0;1(-H2O)", "C20H41N", 296.3312),
                ])
         # Standard LCB fragments (for Cer, SM, and all GSLs)
        else:
            if lcb_type == "16:0;2":
                fragments.extend([
                    ("LCB 16:0;2(-HO)", "C16H33NO", 256.2635),
                    ("LCB 16:0;2(-H3O2)", "C16H31N", 238.2529),
                    ("LCB 16:0;2(-CH3O2)", "C15H31N", 226.2529),
                ])
            elif lcb_type == "16:1;2":
                fragments.extend([
                    ("LCB 16:1;2(-HO)", "C16H31NO", 254.2478),
                    ("LCB 16:1;2(-H3O2)", "C16H29N", 236.2373),
                    ("LCB 16:1;2(-CH3O2)", "C15H29N", 224.2373),
                ])
            elif lcb_type == "17:0;2":
                fragments.extend([
                    ("LCB 17:0;2(-HO)", "C17H35NO", 270.2791),
                    ("LCB 17:0;2(-H3O2)", "C17H33N", 252.2686),
                    ("LCB 17:0;2(-CH3O2)", "C16H33N", 240.2686),
                ])
            elif lcb_type == "17:1;2":
                fragments.extend([
                    ("LCB 17:1;2(-HO)", "C17H33NO", 268.2635),
                    ("LCB 17:1;2(-H3O2)", "C17H31N", 250.2529),
                    ("LCB 17:1;2(-CH3O2)", "C16H31N", 238.2529),
                ])
            elif lcb_type == "18:0;2":
                fragments.extend([
                    ("LCB 18:0;2(-HO)", "C18H37NO", 284.2948),
                    ("LCB 18:0;2(-H3O2)", "C18H35N", 266.2842),
                    ("LCB 18:0;2(-CH3O2)", "C17H35N", 254.2842),
                ])
            elif lcb_type == "18:1;2":
                fragments.extend([
                    ("LCB 18:1;2(-HO)", "C18H35NO", 282.2791),
                    ("LCB 18:1;2(-H3O2)", "C18H33N", 264.2686),
                    ("LCB 18:1;2(-CH3O2)", "C17H33N", 252.2686),
                ])
            elif lcb_type == "18:2;2":
                fragments.extend([
                    ("LCB 18:2;2(-HO)", "C18H33NO", 280.2635),
                    ("LCB 18:2;2(-H3O2)", "C18H31N", 262.2529),
                    ("LCB 18:2;2(-CH3O2)", "C17H31N", 250.2529),
                ])
            elif lcb_type == "18:0;3":
                fragments.extend([
                    ("LCB 18:0;3(-HO)", "C18H37NO2", 300.2897),
                    ("LCB 18:0;3(-H3O2)", "C18H35NO", 282.2791),
                    ("LCB 18:0;3(-CH3O2)", "C17H35NO", 270.2791),
                ])
            elif lcb_type == "19:0;2":
                fragments.extend([
                    ("LCB 19:0;2(-HO)", "C19H39NO", 298.3104),
                    ("LCB 19:0;2(-H3O2)", "C19H37N", 280.2999),
                    ("LCB 19:0;2(-CH3O2)", "C18H37N", 268.2999),
                ])
            elif lcb_type == "19:1;2":
                fragments.extend([
                    ("LCB 19:1;2(-HO)", "C19H37NO", 296.2948),
                    ("LCB 19:1;2(-H3O2)", "C19H35N", 278.2842),
                    ("LCB 19:1;2(-CH3O2)", "C18H35N", 266.2842),
                ])
            elif lcb_type == "20:0;2":
                fragments.extend([
                    ("LCB 20:0;2(-HO)", "C20H41NO", 312.3261),
                    ("LCB 20:0;2(-H3O2)", "C20H39N", 294.3156),
                    ("LCB 20:0;2(-CH3O2)", "C19H39N", 282.3156),
                ])
            elif lcb_type == "20:1;2":
                fragments.extend([
                    ("LCB 20:1;2(-HO)", "C20H39NO", 310.3104),
                    ("LCB 20:1;2(-H3O2)", "C20H37N", 292.2999),
                    ("LCB 20:1;2(-CH3O2)", "C19H37N", 280.2999),
                ])

        return fragments

@dataclass
class AdductInfo:
    name: str
    mass_delta: float
    charge: int
    polarity: str

def get_adduct_definitions(charge_states: List[int] = [1],
                          selected_adducts: Optional[List[str]] = None) -> List[AdductInfo]:
    """Get adduct definitions for specified charge states and adduct types"""
    all_adducts = []

    adduct_map = {
        1: {
            '[M+H]+': AdductInfo("[M+H]1+", PROTON_MASS, 1, "positive"),
            '[M-H]-': AdductInfo("[M-H]1-", -PROTON_MASS, -1, "negative"),
            '[M+Na]+': AdductInfo("[M+Na]1+", 22.98976928, 1, "positive"),
            '[M+NH4]+': AdductInfo("[M+NH4]1+", 18.03383, 1, "positive"),
            '[M+CH3COO]-': AdductInfo("[M+CH3COO]1-", 59.013851, -1, "negative"),
            '[M+HCOO]-': AdductInfo("[M+HCOO]1-", 44.998201, -1, "negative"),
        },
        2: {
            '[M+2H]2+': AdductInfo("[M+2H]2+", 2 * PROTON_MASS, 2, "positive"),
            '[M-2H]2-': AdductInfo("[M-2H]2-", -2 * PROTON_MASS, -2, "negative"),
            '[M+H+Na]2+': AdductInfo("[M+H+Na]2+", PROTON_MASS + 22.98976928, 2, "positive"),
            '[M+2Na]2+': AdductInfo("[M+2Na]2+", 2 * 22.98976928, 2, "positive"),
        },
        3: {
            '[M+3H]3+': AdductInfo("[M+3H]3+", 3 * PROTON_MASS, 3, "positive"),
            '[M-3H]3-': AdductInfo("[M-3H]3-", -3 * PROTON_MASS, -3, "negative"),
            '[M+2H+Na]3+': AdductInfo("[M+2H+Na]3+", 2 * PROTON_MASS + 22.98976928, 3, "positive"),
            '[M+H+2Na]3+': AdductInfo("[M+H+2Na]3+", PROTON_MASS + 2 * 22.98976928, 3, "positive"),
            '[M+3Na]3+': AdductInfo("[M+3Na]3+", 3 * 22.98976928, 3, "positive"),
        },
        4: {
            '[M+4H]4+': AdductInfo('[M+4H]4+', 4*PROTON_MASS, 4, 'positive'),
            '[M-4H]4-': AdductInfo('[M-4H]4-', -4*PROTON_MASS, -4, 'negative'),
        },
        5: {
            '[M+5H]5+': AdductInfo('[M+5H]5+', 5*PROTON_MASS, 5, 'positive'),
            '[M-5H]5-': AdductInfo('[M-5H]5-', -5*PROTON_MASS, -5, 'negative'),
        }
    }

    if selected_adducts is None:
        for charge in charge_states:
            if charge == 1:
                all_adducts.extend([
                    adduct_map[1]['[M+H]+'],
                    adduct_map[1]['[M-H]-'],
                ])
            elif charge == 2:
                all_adducts.extend([
                    adduct_map[2]['[M+2H]2+'],
                    adduct_map[2]['[M-2H]2-'],
                ])
            elif charge == 3:
                all_adducts.extend([
                    adduct_map[3]['[M+3H]3+'],
                    adduct_map[3]['[M-3H]3-'],
                ])
            elif charge == 4:
                all_adducts.extend([
                    adduct_map[4]['[M+4H]4+'],
                    adduct_map[4]['[M-4H]4-'],
                ])
            elif charge == 5:
                all_adducts.extend([
                    adduct_map[5]['[M+5H]5+'],
                    adduct_map[5]['[M-5H]5-'],
                ])

    else:
        for adduct_name in selected_adducts:
            for charge in charge_states:
                if charge in adduct_map and adduct_name in adduct_map[charge]:
                    all_adducts.append(adduct_map[charge][adduct_name])

    return all_adducts

def get_recommended_charges_for_lipid(lipid_class: str) -> List[int]:
    """Get recommended charge states based on lipid class"""
    if LipidDatabase.is_ceramide_class(lipid_class):
        return [1]
    else:
        small_gsl = ['Hex', 'SM4', 'Lac', 'SHex2', 'LC3', 'LC4', 'GA1', 'GA2']
        medium_gsl = ['GM3', 'GM2', 'GM1', 'GD3']
        large_gsl = ['GD2', 'GD1a', 'GD1b', 'GD1c', 'GT3', 'GT2']
        very_large_gsl = ['GT1a', 'GT1b', 'GQ1', 'GP1']

        if lipid_class in small_gsl:
            return [1]
        elif lipid_class in medium_gsl:
            return [1, 2]
        elif lipid_class in large_gsl:
            return [2, 3]
        elif lipid_class in very_large_gsl:
            if lipid_class == 'GP1':
                return [3, 4, 5]
            elif lipid_class == 'GQ1':
                return [3, 4]
            else:  # GT1a, GT1b
                return [2, 3]
        else:
            return [1, 2]

def generate_lipid_formulas(lipid_class: str,
                           selected_lcbs: Optional[List[str]] = None,
                           selected_fatty_acids: Optional[List[str]] = None):
    """Generate molecular formulas for any lipid class with optional selections"""
    formulas = {}

    if selected_lcbs:
        lcb_list = selected_lcbs
    else:
        lcb_list = ConfigManager.get_lcb_list(lipid_class)

    lcb_variations = {}
    for lcb in lcb_list:
        parts = lcb.split(':')
        if len(parts) == 2:
            carbons = int(parts[0])
            unsat_hydox = parts[1].split(';')
            unsaturations = int(unsat_hydox[0])
            hydroxyls = int(unsat_hydox[1]) if len(unsat_hydox) > 1 else (1 if lipid_class == 'doxCer' else 2)

            hydrogens = 2 * carbons + 3 - 2 * unsaturations
            lcb_variations[lcb] = {
                'C': carbons,
                'H': hydrogens,
                'N': 1,
                'O': hydroxyls
            }

    if selected_fatty_acids:
        fa_list = selected_fatty_acids
    else:
        fa_list = ConfigManager.get_fatty_acid_list()

    fatty_acids = {}
    for fa in fa_list:
        parts = fa.split(':')
        if len(parts) == 2:
            carbons = int(parts[0])
            unsaturations = int(parts[1])
            hydrogens = 2 * carbons - 2 * unsaturations
            fatty_acids[fa] = {'C': carbons, 'H': hydrogens, 'O': 2}

    headgroup_comp = LipidDatabase.get_lipid_composition(lipid_class)

    for lcb_name, lcb_comp in lcb_variations.items():
        for fa_name, fa_comp in fatty_acids.items():
            ceramide_backbone = {
                'C': lcb_comp['C'] + fa_comp['C'],
                'H': lcb_comp['H'] + fa_comp['H'] - 2,
                'N': lcb_comp['N'],
                'O': lcb_comp['O'] + fa_comp['O'] - 1
            }

            final_formula = {
                'C': ceramide_backbone['C'] + headgroup_comp['C'],
                'H': ceramide_backbone['H'] + headgroup_comp['H'],
                'N': ceramide_backbone['N'] + headgroup_comp.get('N', 0),
                'O': ceramide_backbone['O'] + headgroup_comp['O'],
                'P': headgroup_comp.get('P', 0),
                'S': headgroup_comp.get('S', 0)
            }

            species_name = f"{lcb_name}/{fa_name}"
            formulas[species_name] = MolecularFormula(final_formula)

    return formulas

def generate_transitions(lipid_class: str, charge_states: List[int] = [1],
                         selected_adducts: Optional[List[str]] = None,
                         selected_lcbs: Optional[List[str]] = None,
                         selected_fatty_acids: Optional[List[str]] = None):
    """
    Generate transitions with m/z calculations for both positive and negative ion modes.

    Parameters:
        lipid_class: Lipid class identifier (e.g., 'Cer', 'SM', 'HexCer', 'SM4')
        charge_states: List of charge states to generate
        selected_adducts: Optional list of specific adducts to use
        selected_lcbs: Optional list of specific LCB types
        selected_fatty_acids: Optional list of specific fatty acid types

    Returns:
        List of transition dictionaries
    """

    logger.info(f"Starting transition generation for {lipid_class} with charge states {charge_states}")

    formulas = generate_lipid_formulas(lipid_class, selected_lcbs, selected_fatty_acids)

    transitions = []
    adducts = get_adduct_definitions(charge_states, selected_adducts)

    is_ceramide = LipidDatabase.is_ceramide_class(lipid_class)
    is_doxcer = (lipid_class == 'doxCer')
    is_gsl = LipidDatabase.is_gsl_class(lipid_class)

    for species_name, formula in formulas.items():
        lcb_type = species_name.split('/')[0]

        # Format precursor name
        if is_ceramide:
            precursor_name = f"{lipid_class}({species_name})"
        else:
            precursor_name = f"{lipid_class} {species_name}"

        formula_str = str(formula)
        precursor_mass = formula.mass()

        for adduct in adducts:
            precursor_mz = (precursor_mass + adduct.mass_delta) / abs(adduct.charge)

            # ============================================================================
            # PRECURSOR TRANSITION
            # ============================================================================
            transitions.append({
                'Molecule List Name': lipid_class,
                'Molecule': precursor_name,
                'Molecule Formula': formula_str,
                'Precursor Adduct': adduct.name,
                'Precursor m/z': round(precursor_mz, 4),
                'Precursor Charge': adduct.charge,
                'Product Name': 'precursor',
                'Product Formula': formula_str,
                'Product Adduct': adduct.name,
                'Product m/z': round(precursor_mz, 4),
                'Product Charge': adduct.charge
            })

            # ============================================================================
            # WATER LOSS (POSITIVE MODE ONLY)
            # ============================================================================
            if adduct.polarity == "positive":
                h2o_loss_formula = MolecularFormula(dict(formula.elements))
                h2o_loss_formula.elements['H'] = h2o_loss_formula.elements.get('H', 0) - 2
                h2o_loss_formula.elements['O'] = h2o_loss_formula.elements.get('O', 0) - 1
                h2o_loss_formula = MolecularFormula({k: v for k, v in h2o_loss_formula.elements.items() if v > 0})
                h2o_loss_formula_str = str(h2o_loss_formula)
                h2o_loss_mass = h2o_loss_formula.mass()
                h2o_loss_mz = (h2o_loss_mass + adduct.mass_delta) / abs(adduct.charge)

                transitions.append({
                    'Molecule List Name': lipid_class,
                    'Molecule': precursor_name,
                    'Molecule Formula': formula_str,
                    'Precursor Adduct': adduct.name,
                    'Precursor m/z': round(precursor_mz, 4),
                    'Precursor Charge': adduct.charge,
                    'Product Name': 'precursor-(H2O,18)',
                    'Product Formula': h2o_loss_formula_str,
                    'Product Adduct': adduct.name,
                    'Product m/z': round(h2o_loss_mz, 4),
                    'Product Charge': adduct.charge
                })

            # Set fragment charge and adduct based on polarity
            fragment_charge = 1 if adduct.charge > 0 else -1
            fragment_adduct = "[M+H]1+" if adduct.charge > 0 else "[M-H]1-"

            # ============================================================================
            # LCB FRAGMENTS (POSITIVE MODE ONLY)
            # ============================================================================
            if adduct.polarity == 'positive':
                if is_ceramide:
                    # Ceramides (including doxCer)
                    lcb_fragments = CeramideFragmentRules.get_ceramide_fragments(lcb_type, is_doxcer)
                elif lipid_class == 'SM':
                    # SM uses ceramide fragments, never doxCer variants
                    lcb_fragments = CeramideFragmentRules.get_ceramide_fragments(lcb_type, is_doxcer=False)
                elif lipid_class == 'SM4':
                    # SM4 (sulfatide) uses ceramide fragments
                    lcb_fragments = CeramideFragmentRules.get_ceramide_fragments(lcb_type, is_doxcer=False)
                elif is_gsl:
                    lcb_fragments = GSLFragmentRules.get_ceramide_fragments(lcb_type)
                else:
                    lcb_fragments = []

                for frag_name, frag_formula, frag_mass in lcb_fragments:
                    frag_mz = frag_mass  # Already [M+H]+
                    transitions.append({
                        'Molecule List Name': lipid_class,
                        'Molecule': precursor_name,
                        'Molecule Formula': formula_str,
                        'Precursor Adduct': adduct.name,
                        'Precursor m/z': round(precursor_mz, 4),
                        'Precursor Charge': adduct.charge,
                        'Product Name': frag_name,
                        'Product Formula': frag_formula,
                        'Product Adduct': fragment_adduct,
                        'Product m/z': round(frag_mz, 4),
                        'Product Charge': fragment_charge
                    })

            # ============================================================================
            # HEADGROUP AND POLARITY-SPECIFIC FRAGMENTS
            # ============================================================================
            if is_gsl or lipid_class in ['SM', 'SM4']:
                if adduct.polarity == 'positive':
                    # === POSITIVE ION MODE ===
                    hg_fragments = GSLFragmentRules.get_headgroup_fragments(lipid_class, formula)

                    for frag in hg_fragments:
                        if len(frag) == 4 and frag[3] is True:
                            frag_name, frag_formula, frag_mass, _ = frag
                            frag_mz = frag_mass  # Already [M+H]+
                        else:
                            frag_name, frag_formula, frag_mass = frag
                            frag_mz = frag_mass + PROTON_MASS  # Single proton only for fragment ions

                        transitions.append({
                            'Molecule List Name': lipid_class,
                            'Molecule': precursor_name,
                            'Molecule Formula': formula_str,
                            'Precursor Adduct': adduct.name,
                            'Precursor m/z': round(precursor_mz, 4),
                            'Precursor Charge': adduct.charge,
                            'Product Name': frag_name,
                            'Product Formula': frag_formula,
                            'Product Adduct': fragment_adduct,
                            'Product m/z': round(frag_mz, 4),
                            'Product Charge': 1  # Always z=1 for product ions
                        })

                elif adduct.polarity == 'negative':
                    # === NEGATIVE ION MODE ===
                    # Negative-mode headgroup fragments (lipid-class specific)
                    hg_fragments_neg = GSLFragmentRules.get_headgroup_fragments_negative(lipid_class, formula)
                    for frag in hg_fragments_neg:
                        if len(frag) == 4:  # 4 elements (with diagnostic flag)
                            frag_name, frag_formula, frag_mass, is_anion = frag
                        else:  # 3 elements (from loss generators)
                            frag_name, frag_formula, frag_mass = frag
                            is_anion = False  # Loss fragments are neutral
                        # Calculate m/z based on fragment type, always singly charged
                        if is_anion:
                            # Already deprotonated (e.g., sialic acid anion)
                            frag_mz = frag_mass
                        elif frag_name.startswith("HG(-"):
                            # Headgroup loss: add H2O for glycosidic cleavage, then deprotonate
                            H2O_MASS = 18.0105647
                            frag_mz = frag_mass + H2O_MASS - PROTON_MASS
                        else:
                            # Other fragments: just deprotonate
                            frag_mz = frag_mass - PROTON_MASS

                        transitions.append({
                            'Molecule List Name': lipid_class,
                            'Molecule': precursor_name,
                            'Molecule Formula': formula_str,
                            'Precursor Adduct': adduct.name,
                            'Precursor m/z': round(precursor_mz, 4),
                            'Precursor Charge': adduct.charge,
                            'Product Name': frag_name,
                            'Product Formula': frag_formula,
                            'Product Adduct': fragment_adduct,
                            'Product m/z': round(frag_mz, 4),
                            'Product Charge': -1  # Always z=-1 for product ions
                        })

                    # === DOUBLY-CHARGED FRAGMENTS FOR GT1 (NEGATIVE MODE) ===
                    if is_gsl and lipid_class in ['GT1a', 'GT1b', 'GT1c'] and abs(adduct.charge) >= 2:
                        logger.debug(f"Processing doubly-charged fragments for {lipid_class} in {adduct.polarity} mode")
                        matched_count = 0

                        for frag in hg_fragments_neg:  # ← Use hg_fragments_neg, not hg_fragments!
                            if len(frag) == 4:
                                frag_name, frag_formula, frag_mass, is_anion = frag
                            else:
                                frag_name, frag_formula, frag_mass = frag

                            # Match single Neu5Ac loss
                            if '-Neu5Ac,309' in frag_name or '-Neu5Ac,291' in frag_name:
                                matched_count += 1
                                logger.debug(f"Matched fragment: {frag_name}")

                                # Negative mode: (M + H2O - 2H) / 2
                                if frag_name.startswith("HG(-"):
                                    H2O_MASS = 18.0105647
                                    doubly_charged_mz = (frag_mass + H2O_MASS - 2 * PROTON_MASS) / 2
                                else:
                                    doubly_charged_mz = (frag_mass - 2 * PROTON_MASS) / 2

                                transitions.append({
                                    'Molecule List Name': lipid_class,
                                    'Molecule': precursor_name,
                                    'Molecule Formula': formula_str,
                                    'Precursor Adduct': adduct.name,
                                    'Precursor m/z': round(precursor_mz, 4),
                                    'Precursor Charge': adduct.charge,
                                    'Product Name': f"{frag_name} [Z=2]",
                                    'Product Formula': frag_formula,
                                    'Product Adduct': '[M-2H]2-',
                                    'Product m/z': round(doubly_charged_mz, 4),
                                    'Product Charge': -2
                                })
                        if matched_count > 0:
                            logger.info(f"Generated {matched_count} doubly-charged fragments for {lipid_class}")

                    if is_gsl and lipid_class == 'GP1' and abs(adduct.charge) >= 2:
                        logger.debug(f"Processing doubly-charged fragments for {lipid_class} in {adduct.polarity} mode")
                        matched_count = 0

                        for frag in hg_fragments_neg:
                            if len(frag) == 4:
                                frag_name, frag_formula, frag_mass, is_anion = frag
                            else:
                                frag_name, frag_formula, frag_mass = frag

                            # Only generate doubly-charged for these losses
                            if frag_name in ["HG(-Neu5Ac,309)", "HG(-Neu5Ac2,600)"]:
                                matched_count += 1
                                logger.debug(f"Matched fragment: {frag_name}")

                                doubly_charged_mz = (frag_mass - 2 * PROTON_MASS) / 2
                                transitions.append({
                                    'Molecule List Name': lipid_class,
                                    'Molecule': precursor_name,
                                    'Molecule Formula': formula_str,
                                    'Precursor Adduct': adduct.name,
                                    'Precursor m/z': round(precursor_mz, 4),
                                    'Precursor Charge': adduct.charge,
                                    'Product Name': f"{frag_name} [Z=2]",
                                    'Product Formula': frag_formula,
                                    'Product Adduct': '[M-2H]2-',
                                    'Product m/z': round(doubly_charged_mz, 4),
                                    'Product Charge': -2
                                })
                        if matched_count > 0:
                            logger.info(f"Generated {matched_count} doubly-charged fragments for {lipid_class}")

                    # === DOUBLY-CHARGED FRAGMENTS FOR nLc SERIES (NEGATIVE MODE) ===
                    if is_gsl and lipid_class in ['nLc10', 'nLc8'] and abs(adduct.charge) >= 2:
                        logger.debug(f"Processing doubly-charged fragments for {lipid_class} in {adduct.polarity} mode")
                        matched_count = 0

                        for frag in hg_fragments_neg:
                            if len(frag) == 4:
                                frag_name, frag_formula, frag_mass, is_anion = frag
                            else:
                                frag_name, frag_formula, frag_mass = frag

                            # nLc10: Match specific fragments for doubly charged
                            if lipid_class == 'nLc10':
                                if frag_name in ["HG(-HexNAc,221)", "HG(-HexNAcHex,383)", "HG(-HexNAc2Hex,586)"]:
                                    matched_count += 1
                                    logger.debug(f"Matched fragment: {frag_name}")

                                    # Negative mode: (M - 2H) / 2
                                    doubly_charged_mz = (frag_mass - 2 * PROTON_MASS) / 2

                                    transitions.append({
                                        'Molecule List Name': lipid_class,
                                        'Molecule': precursor_name,
                                        'Molecule Formula': formula_str,
                                        'Precursor Adduct': adduct.name,
                                        'Precursor m/z': round(precursor_mz, 4),
                                        'Precursor Charge': adduct.charge,
                                        'Product Name': f"{frag_name} [Z=2]",
                                        'Product Formula': frag_formula,
                                        'Product Adduct': '[M-2H]2-',
                                        'Product m/z': round(doubly_charged_mz, 4),
                                        'Product Charge': -2
                                    })

                            # nLc8: Match HG(-Hex,180) for doubly charged
                            elif lipid_class == 'nLc8':
                                if frag_name == "HG(-Hex,180)":
                                    matched_count += 1
                                    logger.debug(f"Matched fragment: {frag_name}")

                                    # Negative mode: (M - 2H) / 2
                                    doubly_charged_mz = (frag_mass - 2 * PROTON_MASS) / 2

                                    transitions.append({
                                        'Molecule List Name': lipid_class,
                                        'Molecule': precursor_name,
                                        'Molecule Formula': formula_str,
                                        'Precursor Adduct': adduct.name,
                                        'Precursor m/z': round(precursor_mz, 4),
                                        'Precursor Charge': adduct.charge,
                                        'Product Name': f"{frag_name} [Z=2]",
                                        'Product Formula': frag_formula,
                                        'Product Adduct': '[M-2H]2-',
                                        'Product m/z': round(doubly_charged_mz, 4),
                                        'Product Charge': -2
                                    })

                        if matched_count > 0:
                            logger.info(f"Generated {matched_count} doubly-charged fragments for {lipid_class}")


                    # Negative-mode LCB fragments (universal for all GSLs)
                    lcb_fragments_neg = NegativeFragmentRules.get_lcb_fragments_negative(lcb_type)

                    for frag_name, frag_formula, frag_mass, _ in lcb_fragments_neg:
                        transitions.append({
                            'Molecule List Name': lipid_class,
                            'Molecule': precursor_name,
                            'Molecule Formula': formula_str,
                            'Precursor Adduct': adduct.name,
                            'Precursor m/z': round(precursor_mz, 4),
                            'Precursor Charge': adduct.charge,
                            'Product Name': frag_name,
                            'Product Formula': frag_formula,
                            'Product Adduct': fragment_adduct,
                            'Product m/z': round(frag_mass, 4),
                            'Product Charge': fragment_charge
                        })

                    # Negative-mode FA fragments (universal for all GSLs)
                    if '/' in species_name:
                        fa_type = species_name.split('/')[1]
                        fa_fragments_neg = NegativeFragmentRules.get_fa_fragments_negative(fa_type)

                        for frag_name, frag_formula, frag_mass, _ in fa_fragments_neg:
                            transitions.append({
                                'Molecule List Name': lipid_class,
                                'Molecule': precursor_name,
                                'Molecule Formula': formula_str,
                                'Precursor Adduct': adduct.name,
                                'Precursor m/z': round(precursor_mz, 4),
                                'Precursor Charge': adduct.charge,
                                'Product Name': frag_name,
                                'Product Formula': frag_formula,
                                'Product Adduct': fragment_adduct,
                                'Product m/z': round(frag_mass, 4),
                                'Product Charge': fragment_charge
                            })

    return transitions


ADDUCT_INSERT_RE = re.compile(r'^\[(?P<prefix>\d*)M')

def blank_mz_values(df: pd.DataFrame) -> pd.DataFrame:
    """
    Blank all m/z values in the dataframe.

    """
    df = df.copy()
    # Target columns that contain m/z values
    target_cols = ['Precursor m/z', 'Product m/z']

    for col in target_cols:
        if col in df.columns:
            df[col] = ''
        else:
            # Check for stripped versions if exact match fails
            stripped_col = col.strip()
            if stripped_col in df.columns:
                df[stripped_col] = ''

    return df

def add_isotope_labels(df: pd.DataFrame, isotope: str = 'M2DN15',
                       doxcer_isotope: str = 'M3D', cer_isotope: str = 'M2DN15',
                       lcb: str = "LCB,precursor,HG(-Hex") -> pd.DataFrame:
    """Add isotope labels WITH calculated m/z values."""

    df.rename(columns=lambda c: c.lstrip('\ufeff').strip(), inplace=True)

    light = df.copy()
    light['Isotope Label Type'] = 'light'

    heavy = df.copy()
    heavy['Isotope Label Type'] = 'heavy'

    # Get isotope token for each lipid class
    def get_isotope_token(list_name):
        if list_name == 'doxCer':
            return doxcer_isotope
        elif list_name == 'Cer':
            return cer_isotope
        else:
            return isotope

    heavy.loc[:, 'Isotope Token'] = heavy['Molecule List Name'].apply(get_isotope_token)

    # Calculate mass shift for each row
    def get_mass_shift(iso_token):
        isotope_composition = parse_isotope_label(iso_token)
        return calculate_isotope_mass_shift(isotope_composition)

    heavy.loc[:, 'Mass Shift'] = heavy['Isotope Token'].apply(get_mass_shift)

    def to_heavy_adduct(adduct, iso_token):
        """
        Insert isotope label into adduct notation with uppercase normalization.

        Examples:
            'M+H' + 'm2dn15' -> 'M[M2DN15+H' (normalized to uppercase)
            'M+H' + 'M2DN15' -> 'M[M2DN15+H' (already uppercase)
        """
        if not isinstance(adduct, str) or not adduct.strip():
            return adduct

        # Normalize isotope token to uppercase for consistent display
        iso_token_normalized = iso_token.strip().upper()

        return ADDUCT_INSERT_RE.sub(
            lambda m: f'[{m.group("prefix")}{iso_token_normalized}',
            adduct,
            count=1
        )

    heavy.loc[:, 'Precursor Adduct'] = [
        to_heavy_adduct(a, iso) for a, iso in
        zip(heavy['Precursor Adduct'], heavy['Isotope Token'])
    ]

    # Calculate PRECURSOR m/z with mass shift
    if 'Precursor m/z' in heavy.columns:
        heavy.loc[:, 'Precursor m/z'] = heavy.apply(
            lambda row: round(row['Precursor m/z'] + row['Mass Shift'] / abs(row['Precursor Charge']), 4)
            if pd.notna(row['Precursor m/z']) and row['Precursor m/z'] != '' else row['Precursor m/z'],
            axis=1
        )
    # Determine which products to label
    keywords = [k.strip() for k in lcb.split(',') if k.strip()]

    def should_label_product(product_name, product_formula):
        product_lower = str(product_name).lower()
        return any(keyword.lower() in product_lower for keyword in keywords)

    # Update product adduct text
    heavy.loc[:, 'Product Adduct'] = [
        to_heavy_adduct(a, iso) if should_label_product(n, f) else a
        for a, n, f, iso in zip(heavy['Product Adduct'], heavy['Product Name'],
                                heavy['Product Formula'], heavy['Isotope Token'])
    ]

    # Calculate PRODUCT m/z for labeled products
    if 'Product m/z' in heavy.columns:
        heavy.loc[:, 'Product m/z'] = heavy.apply(
            lambda row: round(row['Product m/z'] + row['Mass Shift'], 4)
            if pd.notna(row['Product m/z']) and row['Product m/z'] != '' and
               should_label_product(row['Product Name'], row['Product Formula'])
            else row['Product m/z'],
            axis=1
        )

    heavy = heavy.drop(['Isotope Token', 'Mass Shift'], axis=1)

    return pd.concat([light, heavy], ignore_index=True)


def main():
    """CLI interface"""
    parser = argparse.ArgumentParser(
        description='GSL + Ceramide Transition Generator v1.0 - WITH CONFIGURABLE CHAINS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gslgen_v10.py --lipid-class GM1 --charge-states 1 2 --adducts '[M+H]+' '[M+Na]+' --add-labels -v
  python gslgen_v10.py --lipid-class Cer --adducts '[M+H]+' '[M+NH4]+' -v
  python gslgen_v10.py --lipid-class GD3 --charge-states 2 3 -v
        """
    )

    all_classes = LipidDatabase.get_all_classes()
    parser.add_argument('--lipid-class', required=True, choices=all_classes,
                       help='Lipid class to generate')
    parser.add_argument('--output', '-o', help='Output CSV file')
    parser.add_argument('--charge-states', nargs='*', type=int, choices=[1, 2, 3, 4, 5],
                       help='Charge states to include')
    parser.add_argument('--auto-charges', action='store_true',
                       help='Use recommended charge states')

    parser.add_argument('--adducts', nargs='*',
                       choices=['[M+H]+', '[M-H]-', '[M+Na]+', '[M+NH4]+',
                               '[M+CH3COO]-', '[M+HCOO]-',
                               '[M+2H]2+', '[M-2H]2-', '[M+H+Na]2+', '[M+2Na]2+',
                               '[M+3H]3+', '[M-3H]3-', '[M+2H+Na]3+', '[M+H+2Na]3+', '[M+3Na]3+',
                               '[M+4H]4+', '[M-4H]4-',
                               '[M+5H]5+', '[M-5H]5-'],
                       help='Specific adducts to generate')

    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('--add-labels', action='store_true', help='Add isotope labels')
    parser.add_argument('--isotope', default='M2DN15', help='Isotope token for GSLs')
    parser.add_argument('--cer-isotope', default='M2DN15', help='Isotope token for Cer')
    parser.add_argument('--doxcer-isotope', default='M3D', help='Isotope token for doxCer')
    parser.add_argument('--lcb', default='LCB,precursor,HG(-Hex', help='Label keywords')
    parser.add_argument('--blank-mz', action='store_true', help='Blank ALL m/z values in light and heavy rows')

    args = parser.parse_args()

    if not args.output:
        args.output = f"{args.lipid_class.lower()}_transitions_final.csv"

    if not args.charge_states:
        if args.auto_charges:
            args.charge_states = get_recommended_charges_for_lipid(args.lipid_class)
        else:
            args.charge_states = [1]

    print(f"🧬 GSL + Ceramide Transition Generator v1.0")
    print(f"📋 Generating transitions for {args.lipid_class}")

    structure = LipidDatabase.get_structure_description(args.lipid_class)
    if structure != 'Structure not specified':
        print(f"📐 Structure: {structure}")

    print(f"📏 MW range: {LipidDatabase.molecular_weight_range(args.lipid_class)} Da")
    print(f"⚡ Charge states: {args.charge_states}")

    if args.adducts:
        print(f"🔬 Selected adducts: {', '.join(args.adducts)}")

    transitions = generate_transitions(args.lipid_class, args.charge_states, args.adducts)
    df = pd.DataFrame(transitions)

    print(f"✅ Generated {len(transitions)} base transitions")

    if args.add_labels:
        print(f"🏷️ Adding isotope labels...")
        df = add_isotope_labels(df,
                                isotope=args.isotope,
                                doxcer_isotope=args.doxcer_isotope,
                                cer_isotope=args.cer_isotope,
                                lcb=args.lcb)
        print(f"🏷️ Final: {len(df)} transitions (light/heavy)")

    # Apply blanking
    if args.blank_mz:
        print(f"🧹 Blanking all m/z values...")
        df = blank_mz_values(df)

    df.to_csv(args.output, index=False)
    print(f"💾 Saved to: {os.path.abspath(args.output)}")

    if args.verbose:
        print(f"\n📈 Sample transitions:")
        sample_cols = ['Molecule', 'Precursor Adduct', 'Precursor m/z', 'Product Name', 'Product m/z']
        if args.add_labels:
            sample_cols.append('Isotope Label Type')
        print(df[sample_cols].head(10).to_string(index=False))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main()
    else:
        print("🧬 GSL + Ceramide Transition Generator v1.0")
        print("\nUsage:")
        print("  python gslgen_v10.py --lipid-class GM1 --charge-states 1 2 --adducts '[M+H]+' '[M+Na]+' -v")
        print("  python gslgen_v10.py --lipid-class Cer --adducts '[M+H]+' '[M+NH4]+' -v")
        print("\nRun with --help for all options")

