GSL & Ceramide Transition Generator (v1.0)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17702284.svg)](https://doi.org/10.5281/zenodo.17702284)

User Manual & Documentation
Author: Andreas J. Hülsmeier
License: MIT

1. Overview

The GSL & Ceramide Transition Generator is a Python-based tool designed for the targeted mass spectrometry analysis of sphingolipids.
Developed as an independent refinement of the [LipidCreator](https://lifs-tools.org/lipidcreator) concept ([Peng et al., 2020, *Nat Commun*, 11, 2057](https://doi.org/10.1038/s41467-020-15960-z)), this application moves beyond general lipidomics to focus exclusively on Glycosphingolipids (GSLs) and Ceramides. It addresses the specific analytical challenges of high-molecular-weight gangliosides (up to GP1) by offering high-precision mass calculation, fragmentation rules, and extensive multi-charge state support.

Key Features

-Specialized Coverage: Verified molecular formulas for 30+ classes, from simple Ceramides to complex Gangliosides (GM, GD, GT, GQ, GP series).

-Math calculations based on IUPAC 2016 atomic masses.

-User-defined ranges for Fatty Acids (FA) and Long Chain Bases (LCB), including hydroxylation and unsaturation specificities.

-Isotope Labeling: Support for metabolic labeling studies (Deuterium, ¹³C, ¹⁵N, ¹⁸O).

-Skyline Compatibility: Generates native .csv transition lists ready for import into [Skyline](https://skyline.ms) MS.


2. Installation & Requirements

System Requirements

-OS: Windows, macOS, or Linux (Cross-platform Python application)

-Python: Version 3.8 or higher

-Screen Resolution: Minimum 1280x720 recommended for GUI usage

Dependencies

The application relies on the following Python libraries:
-pandas (Data handling)

-PySide6 (Graphical User Interface)

Installation Steps

I. Ensure Python 3 is installed on your system.
II. Install the required dependencies via terminal/command prompt:
    
bash
pip install pandas PySide6

III. Place gslgen.py (Backend), gui_gslgen.py (Frontend) and icon_256.png in the same directory.
IV. Launch the application:

    
GUI user interface

bash
python gui_gslgen.py

Basic examples for cli usage

python gslgen.py --lipid-class GM1 --charge-states 1 2 --output gm1_transitions.csv

With isotope labels

python gslgen.py --lipid-class Cer --add-labels --output cer_labeled.csv

See all options

python gslgen.py --help

3. User Interface Guide

The GUI is divided into three logical sections: Selection, Configuration, and Output.

A. Lipid Class Selection

Located at the top-left. The dropdown menu allows selection of the specific lipid headgroup chemistry.

-Scope: Includes neutral GSLs (GlcCer, LacCer, Gb3/4), acidic Gangliosides (GM1-3, GD1-3, GT1-3, GQ1, GP1), Ceramides, 1-deoxy-Ceramides and Sphingomyelins.

-Note: Selecting a class automatically loads the related fragmentation rules (e.g., sialic acid loss for gangliosides, phosphocholine for SM, LCB fragments for ceramides).

    
B. MS Parameters

-Charge States: Checkboxes for 1 to 5.

Recommendation: Use lower states (1) for neutral GSLs and high states (2, 3, 4, 5) for polysialylated gangliosides (e.g., GD1b, GT1c or GQ1) to bring precursors into measurable m/z ranges.

-Adducts: Select the ionization mode.

Positive: [M+H]+, [M+Na]+, [M+NH4]+

Negative: [M-H]-, [M+HCOO]- (Formate), [M+CH3COO]- (Acetate)

C. Isotope Labeling (Optional)

Enable "Isotope Labeling" to generate transitions for metabolically labeled samples.

-Elements: Support for Deuterium (D), ¹⁵N, ¹³C.

-Configuration: The specific number of labeled atoms per molecule is defined in the backend configuration (see Section 4).

-"Label Keywords" determine which products get isotope labeling. Enter a comma-separated list of substrings (case-insensitive) that the script will search for in each "Product Name". When a "Product Name" contains any of these keywords, the script marks it as label-eligible and creates an isotope-shifted "heavy" version.
The default keywords are LCB, precursor, HG(-Hex, but you can edit this list to control exactly which fragment groups receive labeling in your exported transition list.

4. Configuration

Unlike static lists, the GSL & Ceramide Transition Generator builds molecules dynamically. You must configure the building blocks via the Configuration button (or Edit > Settings).

Long Chain Base (LCB) Settings

Select the sphingoid backbone.

-Carbon Range: e.g., 18 to 20 (Generates d18:1, d20:1).

-Hydroxylation e.g.:

18:0;1 = 1-deoxy-sphinganine

18:1;2 = Sphingosine

18:0;3 = Phytosphingosine


Fatty Acid (FA) Settings

Define the N-acyl chain.

-Range: Min/Max carbon length (e.g., 16 to 24).

-Unsaturation: Max double bonds (e.g., 0 for saturated only, 1 for mono-unsaturated).

-Parity: Choose Even only (biological standard), Odd, or Both.

    
Tip: Keep these ranges narrow. Selecting "All Carbon Lengths" combined with "All Classes" will generate hundreds of thousands of transitions and may freeze the application. Choose specific ranges that match your biological hypothesis.

5. Workflow: Generating a Method

    1. Launch the GUI.
    2. Configure your LCB and FA ranges (e.g., d18:1, FA 16:0-24:0).
    3. Select Lipid Class (e.g., GD1a).
    4. Select Charge States: For GD1a (2 sialic acids), select -2 and -3.
    5. Select Adducts: Choose [M-H]-.
    6. Click "Generate Transitions".
        ◦ The software calculates exact masses using IUPAC values.
        ◦ It applies fragmentation rules (e.g., Precursor - 291Da for sialic acid loss).
    7. Review the results in the table view.
    8. Export to CSV.

6. Integration with Skyline MS. This tool uses a File-Based Integration workflow.

    1. Export: Save your generated list as your_filename.csv from the GSL Generator.
    2. Open Skyline: Create a new Small Molecule document.
    3. Import:
        ◦ Go to File > Import > Transition List...
        ◦ Select your CSV file.
    4. Verify: Skyline will recognize the headers (Precursor m/z, Product m/z, Molecule Name) automatically.

CSV Output Format

The application exports the following columns:

-Molecule List Name: The lipid class selected (e.g., GD1a).

-Molecule: The specific species name (e.g., GD1a 18:1;2/18:0).

-Molecule Formula: The molecular sum formula of the precursor.

-Precursor Adduct: The selected precursor adduct and charge state.

-Precursor m/z: Exact mass of the precursor ion.

-Precursor Charge: The selected charge state for the precursor ion.

-Product Name: The name of the product ion.

-Product Formula: The molecular sum formula of the fragment.

-Product m/z: Exact mass of the fragment.

-Product Charge: The charge state of the fragment.

-Label: Light/Heavy status, when isotope labels are selected.



7. Scientific Rationale & Calculations

Exact Mass Calculation

Masses are calculated using IUPAC 2016 standards:

-Carbon: 12.0000000

-Hydrogen: 1.00782503223

-Nitrogen: 14.0030740048

-Oxygen: 15.9949146223

-Proton: 1.007276466812


Fragmentation Logic

The tool implements verified fragmentation pathways:

-Ceramides: LCB fragments (e.g., m/z 264.27 for d18:1).

-GSLs: Glycosidic bond cleavages (B-, C-, Y-, Z-ions) and headgroup losses.

-Gangliosides: Specific logic for sialic acid loss (Neu5Ac, -291 Da) and characteristic backbone fragments.


8. Troubleshooting

Issue: The application is slow/unresponsive during generation.

-Cause: You may have selected too wide a range of Fatty Acids or LCBs combined with too many charge states.

-Fix: Narrow the configuration, be specific to your biological hypothesis.

Issue: Masses in Skyline slightly differ.

-Cause: Rounding differences. Skyline uses its own atomic weight table.

-Fix: If importing into Skyline, check "Blank all m/z values" to let Skyline recalculate m/z from the explicit formula. This approach leverages Skyline's internal mass calculation engine, which calculates precise m/z values and isotope distributions from molecular formulas. You can then choose between monoisotopic and average mass calculations for both precursor and product ions using Skyline's Mass Type Choice feature. Alternatively, you can manually delete the m/z values from the exported CSV file before importing, or simply trust the generator's m/z calculations, which are based on IUPAC-standard molecular formulas.

© 2025 Andreas J. Hülsmeier
University of Zurich / University Hospital Zurich
