import pandas as pd
import logging
from typing import List, Optional
from dataclasses import dataclass

# Import from your new modules
from constants import PROTON_MASS
from chemistry import MolecularFormula
from database import LipidDatabase, ConfigManager
from fragment_rules import CeramideFragmentRules, GSLFragmentRules, NegativeFragmentRules

logger = logging.getLogger(__name__)

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
                    # Negative-mode headgroup fragments
                    hg_fragments_neg = GSLFragmentRules.get_headgroup_fragments_negative(lipid_class, formula)
                    for frag in hg_fragments_neg:
                        if len(frag) == 4:  # 4 elements (with diagnostic flag)
                            frag_name, frag_formula, frag_mass, is_anion = frag
                        else:  # 3 elements (from loss generators)
                            frag_name, frag_formula, frag_mass = frag
                            is_anion = False  # Loss fragments are neutral
                        # Calculate m/z based on fragment type
                        if is_anion:
                            # Already deprotonated (e.g., sialic acid anion)
                            frag_mz = frag_mass
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

                                # Negative mode: [M-2H]2- calculation
                                # (Neutral_Mass - 2 * Proton_Mass) / 2
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

                    # Check is_anion flag
                    for frag_name, frag_formula, frag_mass, is_anion in lcb_fragments_neg:

                        # Calculate m/z based on flag
                        if is_anion:
                            frag_mz = frag_mass
                        else:
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
                            'Product Charge': fragment_charge
                        })

                    # Negative-mode FA fragments (universal for all GSLs)
                    if '/' in species_name:
                        fa_type = species_name.split('/')[1]
                        fa_fragments_neg = NegativeFragmentRules.get_fa_fragments_negative(fa_type)

                        # Check is_anion flag
                        for frag_name, frag_formula, frag_mass, is_anion in fa_fragments_neg:

                            # Calculate m/z based on flag
                            if is_anion:
                                frag_mz = frag_mass
                            else:
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
                                'Product Charge': fragment_charge
                            })

    return transitions

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
        small_gsl = ['Hex', 'SM4', 'Lac', 'SHex2', 'LC3', 'LC4', 'GA1', 'GA2', 'GM3', 'GM2', 'GM1']
        medium_gsl = ['GD3']
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
