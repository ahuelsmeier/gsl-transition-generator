# GSL & Ceramide Transition Generator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17702284.svg)](https://doi.org/10.5281/zenodo.17702284)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A Python-based tool for generating high-precision mass spectrometry transition lists for glycosphingolipids (GSLs) and ceramides, with specialized support for complex gangliosides up to GP1.

##  Key Features

- **Specialized Coverage**: 30+ lipid classes from simple ceramides to complex gangliosides (GM, GD, GT, GQ, GP series)
- **High-Precision Calculations**: IUPAC 2016 atomic masses with exact mass computation
- **Multi-Charge Support**: Charge states 1+ to 5+ for high-molecular-weight species
- **Isotope Labeling**: Support for ²H, ¹³C, ¹⁵N, ¹⁸O metabolic labeling studies
- **[Skyline](https://skyline.ms) Compatible**: Native CSV export for direct import into Skyline MS
- **Dual Interface**: Command-line (CLI) and graphical (GUI) options
- **Dynamic Configuration**: User-defined ranges for fatty acids and long-chain bases

##  Installation

### Requirements

- Python 3.8 or higher
- Dependencies: `pandas`, `PySide6`

### Quick Start

Install dependencies

pip install pandas PySide6

Clone repository

git clone https://github.com/ahuelsmeier/gsl-transition-generator.git
cd gsl-transition-generator

Launch GUI

python gui_gslgen.py


##  Usage

### Graphical Interface (Recommended)

python gui_gslgen.py


### Command-Line Interface

Basic example

python gslgen.py --lipid-class GM1 --charge-states 1 2 --output gm1_transitions.csv

With isotope labels

python gslgen.py --lipid-class Cer --add-labels --output cer_labeled.csv

View all options

python gslgen.py --help


##  Documentation

For detailed usage instructions, configuration options, and scientific rationale, see the **[User Manual](USER_MANUAL.md)**.

### Quick Links
- [Supported Lipid Classes](USER_MANUAL.md#lipid-class-selection)
- [Configuration Guide](USER_MANUAL.md#configuration)
- [Skyline Integration](USER_MANUAL.md#integration-with-skyline-ms)
- [Troubleshooting](USER_MANUAL.md#troubleshooting)

##  Supported Lipid Classes

**Neutral GSLs**: GlcCer, LacCer, Gb3, Gb4  
**Gangliosides**: GM3, GM2, GM1, GD3, GD2, GD1a/b, GT3, GT1a/b/c, GQ1a/b, GP1  
**Ceramides**: Cer, doxCer  
**Others**: Sphingomyelin (SM)

##  Citation

If you use this software in your research, please cite:

> Hülsmeier, A.J. (2025). GSL & Ceramide Transition Generator (v1.0.5). Zenodo. https://doi.org/10.5281/zenodo.17702284

<details>
<summary>BibTeX format</summary>
@software{huelsmeier_2025_gsl,
author = {Hülsmeier, Andreas J.},
title = {GSL & Ceramide Transition Generator},
year = 2025,
publisher = {Zenodo},
version = {v1.0.5},
doi = {10.5281/zenodo.17702284},
url = {https://doi.org/10.5281/zenodo.17702284}
}

</details>

##  Author

**Andreas J. Hülsmeier**  
University of Zurich, University Hospital Zurich  
 andreas.huelsmeier@uzh.ch  
 [ORCID: 0000-0001-6987-8979](https://orcid.org/0000-0001-6987-8979)

##  License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

##  Acknowledgments

Developed as an independent refinement of the [LipidCreator](https://lifs-tools.org/lipidcreator) concept ([Peng et al., 2020, *Nat Commun*, 11, 2057](https://doi.org/10.1038/s41467-020-15960-z)), focusing exclusively on glycosphingolipids and ceramides with enhanced support for high-molecular-weight gangliosides.

---

**⭐ If you find this tool useful, please consider citing it and starring the repository!**
