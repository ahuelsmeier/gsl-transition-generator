import pandas as pd
import re
from constants import ADDUCT_INSERT_RE
from chemistry import parse_isotope_label, calculate_isotope_mass_shift, filter_isotope_token

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
                       lcb: str = "LCB,precursor,HG(-") -> pd.DataFrame:
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

    # ====================================================================
    # 1. TOKEN FILTERING ENGINE
    # Filters the requested Isotope Token based on the atoms actually present in the formula.
    # ====================================================================


    # ====================================================================
    # 2. APPLY FILTER TO PRECURSORS
    # ====================================================================
    heavy['Filtered_Precursor_Token'] = heavy.apply(
        lambda r: filter_isotope_token(r['Isotope Token'], r.get('Molecule Formula', '')), axis=1
    )

    heavy.loc[:, 'Precursor Adduct'] = [
        to_heavy_adduct(a, tok) if tok else a
        for a, tok in zip(heavy['Precursor Adduct'], heavy['Filtered_Precursor_Token'])
    ]

    if 'Precursor m/z' in heavy.columns:
        heavy.loc[:, 'Precursor m/z'] = heavy.apply(
            lambda row: round(row['Precursor m/z'] + (calculate_isotope_mass_shift(parse_isotope_label(row['Filtered_Precursor_Token'])) / abs(row['Precursor Charge'])), 4)
            if pd.notna(row['Precursor m/z']) and row['Precursor m/z'] != '' and row['Filtered_Precursor_Token'] else row['Precursor m/z'],
            axis=1
        )

    # ====================================================================
    # 3. APPLY FILTER TO PRODUCTS
    # ====================================================================
    keywords = [k.strip() for k in lcb.split(',') if k.strip()]

    def should_label_product(product_name, product_formula):
        product_lower = str(product_name).lower()
        return any(keyword.lower() in product_lower for keyword in keywords)

    heavy['Filtered_Product_Token'] = heavy.apply(
        lambda r: filter_isotope_token(r['Isotope Token'], r.get('Product Formula', '')), axis=1
    )

    heavy.loc[:, 'Product Adduct'] = [
        to_heavy_adduct(a, tok) if should_label_product(n, f) and tok else a
        for a, n, f, tok in zip(heavy['Product Adduct'], heavy['Product Name'], heavy['Product Formula'], heavy['Filtered_Product_Token'])
    ]

    if 'Product m/z' in heavy.columns:
        heavy.loc[:, 'Product m/z'] = heavy.apply(
            lambda row: round(row['Product m/z'] + (calculate_isotope_mass_shift(parse_isotope_label(row['Filtered_Product_Token'])) / abs(row['Product Charge'])), 4)
            if pd.notna(row['Product m/z']) and row['Product m/z'] != '' and row['Filtered_Product_Token'] and should_label_product(row['Product Name'], row['Product Formula']) else row['Product m/z'],
            axis=1
        )

    heavy = heavy.drop(['Isotope Token', 'Mass Shift', 'Filtered_Precursor_Token', 'Filtered_Product_Token'], axis=1, errors='ignore')

    return pd.concat([light, heavy], ignore_index=True)
