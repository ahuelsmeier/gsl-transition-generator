from chemistry import MolecularFormula
from constants import ATOMIC_MASSES, PROTON_MASS

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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
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
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8},      # dehydrated Neu5Ac
            "HG(-Neu5AcHexNAcHex,674)": {'C': 25, 'H': 42, 'N': 2, 'O': 19}, # not dehydrated
            "HG(-Neu5AcHexNAcHex,656)": {'C': 25, 'H': 40, 'N': 2, 'O': 18}, # dehydrated
            "HG(-Neu5AcHexNAcHex2,836)": {'C': 31, 'H': 52, 'N': 2, 'O': 24}, # not dehydrated
            "HG(-Neu5AcHexNAcHex2,818)": {'C': 31, 'H': 50, 'N': 2, 'O': 23}, # dehydrated
            "HG(-Neu5AcHexNAcHex3,998)": {'C': 37, 'H': 62, 'N': 2, 'O': 29}, # Entire GM1 headgroup, not dehydrated
            "HG(-Neu5AcHexNAcHex3,980)": {'C': 37, 'H': 60, 'N': 2, 'O': 28}, # Entire GM1 headgroup, dehydrated
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
            "HG(-HexNAc,221)": {'C': 8, 'H': 15, 'N': 1, 'O': 6},  # Single HexNAc
            "HG(-HexNAc,203)": {'C': 8, 'H': 13, 'N': 1, 'O': 5},  # Single HexNAc -H2O, Y2β
            "HG(-Neu5AcHexNAc,512)": {'C': 19, 'H': 32, 'N': 2, 'O': 14},  # Neu5Ac + HexNAc
            "HG(-Neu5Ac2HexNAc,803)": {'C': 30, 'H': 49, 'N': 3, 'O': 22},  # 2 NeuAc + HexNAc, Z-ion
            "HG(-Neu5Ac2HexNAc,785)": {'C': 30, 'H': 47, 'N': 3, 'O': 21},  # 2 NeuAc + HexNAc, Y-ion
            "HG(-Neu5Ac2HexNAcHex,965)": {'C': 36, 'H': 59, 'N': 3, 'O': 27},  # Hex + HexNAc+ 2 Neu5Ac, Z-ion
            "HG(-Neu5Ac2HexNAcHex,947)": {'C': 36, 'H': 57, 'N': 3, 'O': 26},  # Hex + HexNAc+ 2 Neu5Ac, Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,618)": {'C': 22, 'H': 38, 'N': 2, 'O': 18},   # 2 x Neu5Ac
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
            "HG(-Neu5Ac2Hex,762)": {'C': 28, 'H': 46, 'N': 2, 'O': 22},  # Neu5Acα2-8NeuAc, Hex
            "HG(-Neu5Ac2HexNAcHex,965)": {'C': 36, 'H': 59, 'N': 3, 'O': 27},  # Hex + HexNAc+ 2 Neu5Ac, Z-ion
            "HG(-Neu5Ac2HexNAcHex,947)": {'C': 36, 'H': 57, 'N': 3, 'O': 26},  # Hex + HexNAc+ 2 Neu5Ac, Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
            "HG(-Neu5Ac2,618)": {'C': 22, 'H': 38, 'N': 2, 'O': 18},  # 2 x Neu5Ac
            "HG(-Neu5Ac3,909)": {'C': 33, 'H': 55, 'N': 3, 'O': 26},  # Neu5Acα2-8NeuAc and Neu5Ac
            "HG(-Neu5Ac3,891)": {'C': 33, 'H': 53, 'N': 3, 'O': 25},  # Neu5Acα2-8NeuAc and Neu5Ac -H2O, YY fragment
            "HG(-Neu5Ac2HexNAcHex,965)": {'C': 36, 'H': 59, 'N': 3, 'O': 27},  # Hex + HexNAc+ 2 Neu5Ac, Z-ion
            "HG(-Neu5Ac2HexNAcHex,947)": {'C': 36, 'H': 57, 'N': 3, 'O': 26},  # Hex + HexNAc+ 2 Neu5Ac, Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
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
            "HG(-Neu5Ac,309)": {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # Neu5Ac Z-ion
            "HG(-Neu5Ac,291)": {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # Neu5Ac Y-ion
            "HG(-Neu5Ac2,600)": {'C': 22, 'H': 36, 'N': 2, 'O': 17},  # NeuAcα2-8NeuAc, Z-ion
            "HG(-Neu5Ac2,582)": {'C': 22, 'H': 34, 'N': 2, 'O': 16},  # NeuAcα2-8NeuAc, Y-ion
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
            "HG(-HexNAc4Hex4,1460)": {'C': 56, 'H': 92, 'N': 4, 'O': 40}  # HexNAc4Hex4 Y-ion
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
                ("Phosphocholine", "C5H14NO4P", 184.0739, True),
                ("Phosphocholine-H2O", "C5H12NO3P", 166.0628, True),
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
                    ("HG(Neu5Ac2HexNAcHex2,1109)", "C42H67N3O31", 1110.3831, True),     # Full HG, NeuAc2 + HexNAc + Hex2 -H2O
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
        Generate fatty acid fragments for negative ion mode.
        Dominant fragment: Carboxylate anion [RCOO]-
        """
        fragments = []

        import re
        match = re.match(r'(\d+):(\d+)', fa_type)
        if not match:
            return fragments

        carbon = int(match.group(1))
        double_bonds = int(match.group(2))
        h_count = 2 * carbon - 2 * double_bonds  # H count for Neutral FA (CnH2nO2)

        # =====================================================================
        # 1. Carboxylate Anion [R-COO]-
        # =====================================================================
        # We define the NEUTRAL fatty acid (R-COOH).
        # The main loop sees the flag is False, so it calculates: (Mass - Proton)
        # resulting in the correct [R-COO]- m/z.

        neutral_formula = f"C{carbon}H{h_count}O2"

        neutral_fa_mass = (
            carbon * ATOMIC_MASSES['C'] +
            h_count * ATOMIC_MASSES['H'] +
            2 * ATOMIC_MASSES['O']
        )

        fragments.append((f"FA {fa_type} [RCOO]-", neutral_formula, neutral_fa_mass, False))

        # =====================================================================
        # 2. Amide/Adduct Fragments (Minor/Diagnostic)
        # =====================================================================
        # FA+(HN) - Retains amide fragment with carbonyl
        fa_hn_formula = f"C{carbon}H{h_count+1}NO"
        fa_hn_mass = (
            carbon * ATOMIC_MASSES['C'] +
            (h_count + 1) * ATOMIC_MASSES['H'] +
            1 * ATOMIC_MASSES['N'] +
            1 * ATOMIC_MASSES['O'] -
            ATOMIC_MASSES['H']  # Pre-calculated [M-H]-
        )

        # FA+(C2H3N) - Acetonitrile adduct
        fa_acn_formula = f"C{carbon+2}H{h_count+3}NO"
        fa_acn_mass = (
            (carbon + 2) * ATOMIC_MASSES['C'] +
            (h_count + 3) * ATOMIC_MASSES['H'] +
            1 * ATOMIC_MASSES['N'] +
            1 * ATOMIC_MASSES['O'] -
            ATOMIC_MASSES['H']  # Pre-calculated [M-H]-
        )

        # FA+(C2H3NO) - Acetamide-like adduct
        fa_c2h3no_formula = f"C{carbon+2}H{h_count+3}NO2"
        fa_c2h3no_mass = (
            (carbon + 2) * ATOMIC_MASSES['C'] +
            (h_count + 3) * ATOMIC_MASSES['H'] +
            1 * ATOMIC_MASSES['N'] +
            2 * ATOMIC_MASSES['O'] -
            ATOMIC_MASSES['H']  # Pre-calculated [M-H]-
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
