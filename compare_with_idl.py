"""
compare_with_idl.py

Load the IDL save file 'N2264.save' and compare the IDL identifications
(N2264_id array) with our CSV master catalog ('ngc2264_combined.csv').
Reports match statistics for each catalog.
"""

import os
os.environ['IDL_STARTUP'] = ''
import numpy as np
import pandas as pd
import sys
sys.path.append('/Applications/NV5/idl/lib/bridges')
from idlpy import *


def load_idl_data():
    """Load IDL save file via IDL bridge and return N2264_id array and lists struct."""

    # idl = IDL()  # start IDL session
    # Use READSAV to load only the needed variables and avoid side-effects
    # idl.execute(f'READSAV, "{SAVE_PATH}", N2264_id=N2264_id, N2264_lists=N2264_lists')
    # Retrieve variables from IDL
    # N2264_id = idl.get_value('N2264_id')
    # N2264_lists = idl.get_value('N2264_lists')
    # idl.exit()  # close IDL session
    
    IDL.run('restore,"/Users/ettoref/ASTRONOMY/DATA/N2264_XMM_alt/N2264.save",/v')
    N2264_id = IDL.N2264_id
    N2264_lists = IDL.N2264_lists

    return N2264_id, N2264_lists

def load_csv_catalog():
    """Load our combined CSV catalog into a DataFrame."""
    df = pd.read_csv('/Users/ettoref/ASTRONOMY/DATA/ngc2264-catalog/ngc2264_combined.csv', dtype=str)
    # Convert numeric ID columns to integers, filling missing with -1
    for col in ['gaia_id', '2MASS', 'wise_id', 'chandra_id', 'xmm_id']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(-1).astype(int)
    return df

def compare_catalogs(idl_ids, df, field_name, colname, index):
    """
    Compare IDL identification numbers vs DataFrame IDs.

    idl_ids: numpy array of shape (N,) from IDL N2264_id[:, index]
    df: pandas DataFrame with our IDs in df[colname]
    field_name: name of the catalog, e.g. 'gaia'
    colname: column in df, e.g. 'gaia_id'
    index: zero-based column index into idl_ids
    """
    idl_col = idl_ids[:, index]
    script_col = df[colname].values
    matches = idl_col == script_col
    total = len(matches)
    nmatch = np.count_nonzero(matches)
    nmismatch = total - nmatch
    print(f"{field_name.upper():10s}: {nmatch}/{total} matches, {nmismatch} mismatches")
    if nmismatch > 0:
        bad = np.where(~matches)[0]
        print("  Example mismatches at row indices:", bad[:10].tolist())

def main():
    # Load data
    print("Loading IDL data...")
    idl_ids, lists = load_idl_data()
    print("Loading combined CSV catalog...")
    df = load_csv_catalog()

    # Define mapping from IDL struct fields to DataFrame columns and index
    mapping = {
        'gaia':    ('gaia_id',    lists.gaia.i - 1),
        'tmass':   ('2MASS',      lists.tmass.i - 1),
        'wise':    ('wise_id',    lists.wise.i - 1),
        'chandra': ('chandra_id', lists.chandra.i - 1),
        'xmm':     ('xmm_id',     lists.xmm.i - 1),
    }

    # Compare for each catalog
    print("Comparing identification fields:")
    for field, (colname, idx) in mapping.items():
        if colname in df.columns:
            compare_catalogs(idl_ids, df, field, colname, idx)
        else:
            print(f"{field.upper():10s}: column '{colname}' not in CSV.")

if __name__ == '__main__':
    main()