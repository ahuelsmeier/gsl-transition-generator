import re

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

# Regex pattern for inserting isotope labels into adducts safely
ADDUCT_INSERT_RE = re.compile(r'(?P<prefix>[M\]\+\-])')
