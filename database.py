from typing import Dict, List, Tuple, Optional
import json
from pathlib import Path

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
            'Cer': 'Ceramide',
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

class ConfigManager:
    """Manage user configuration for LCB and fatty acid selections"""

    # Define the Universe of supported LCBs (Source of Truth)
    SUPPORTED_STANDARD_LCBS = [
        "16:0;2", "16:1;2", "17:0;2", "17:1;2",
        "18:0;2", "18:0;3", "18:1;2", "18:2;2",
        "19:0;2", "19:1;2", "20:0;2", "20:1;2"
    ]

    SUPPORTED_DOX_LCBS = [
        "16:0;1", "16:1;1", "17:0;1", "17:1;1",
        "18:0;1", "18:1;1", "18:2;1",
        "19:0;1", "19:1;1", "20:0;1", "20:1;1"
    ]

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
                "label_keywords": "LCB,precursor,HG(-",
                "blank_mz": False
            }
        }
