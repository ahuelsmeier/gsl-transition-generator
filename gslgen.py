#!/usr/bin/env python3

"""
GSL + Ceramide Transition Generator - Version 1.1.0

author: Andreas J. Hülsmeier
copyright: Copyright 2025, Andreas J. Hülsmeier / University of Zurich, University Hospital Zurich
license: MIT
version: 1.1.0
maintainer: Andreas J. Hülsmeier
email: andreas.huelsmeier@uzh.ch
status: Prototype
"""

import pandas as pd
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
logger = logging.getLogger(__name__)

# ============================================================================
# MODULE ROUTER (Pass-throughs for the GUI)
# ============================================================================
from constants import ATOMIC_MASSES, ISOTOPE_DELTAS, ADDUCT_INSERT_RE
from chemistry import parse_isotope_label, calculate_isotope_mass_shift, filter_isotope_token, MolecularFormula
from fragment_rules import GSLFragmentRules, NegativeFragmentRules, CeramideFragmentRules
from database import LipidDatabase, ConfigManager
from core_transitions import generate_transitions, get_recommended_charges_for_lipid
from isotope_labeling import add_isotope_labels, blank_mz_values

# Expose these specifically for the GUI's configuration dialog
if hasattr(ConfigManager, 'SUPPORTED_STANDARD_LCBS'):
    SUPPORTED_STANDARD_LCBS = ConfigManager.SUPPORTED_STANDARD_LCBS
if hasattr(ConfigManager, 'SUPPORTED_DOX_LCBS'):
    SUPPORTED_DOX_LCBS = ConfigManager.SUPPORTED_DOX_LCBS

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================
def main():
    """CLI interface"""
    parser = argparse.ArgumentParser(
        description='GSL + Ceramide Transition Generator v1.1.0 - WITH CONFIGURABLE CHAINS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gslgen.py --lipid-class GM1 --charge-states 1 2 --adducts '[M+H]+' '[M+Na]+' --add-labels -v
  python gslgen.py --lipid-class Cer --adducts '[M+H]+' '[M+NH4]+' -v
  python gslgen.py --lipid-class GD3 --charge-states 2 3 -v
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
    parser.add_argument('--lcb', default='LCB,precursor,HG(-', help='Label keywords')
    parser.add_argument('--blank-mz', action='store_true',
                       help='Blank ALL m/z values in light and heavy rows')

    args = parser.parse_args()

    if not args.output:
        args.output = f"{args.lipid_class.lower()}_transitions_final.csv"

    if not args.charge_states:
        args.charge_states = (
            get_recommended_charges_for_lipid(args.lipid_class)
            if args.auto_charges
            else [1]
        )

    print(f"🧬 GSL + Ceramide Transition Generator v1.1.0")
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
        print(f"🏷️  Adding isotope labels...")
        df = add_isotope_labels(
            df,
            isotope=args.isotope,
            doxcer_isotope=args.doxcer_isotope,
            cer_isotope=args.cer_isotope,
            lcb=args.lcb
        )
        print(f"🏷️  Final: {len(df)} transitions (light/heavy)")

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
        print("🧬 GSL + Ceramide Transition Generator v1.1.0")
        print("\nUsage:")
        print("  python gslgen.py --lipid-class GM1 --charge-states 1 2 --adducts '[M+H]+' '[M+Na]+' -v")
        print("  python gslgen.py --lipid-class Cer --adducts '[M+H]+' '[M+NH4]+' -v")
        print("\nRun with --help for all options")
