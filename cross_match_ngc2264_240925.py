"""
cross_match_ngc2264.py
======================

This script builds a multi‑wavelength catalogue for the NGC 2264 region by
querying several astronomical archives and cross‑matching the sources on
the sky.  It uses Gaia as the master list and links each Gaia object to
any 2MASS, WISE, Chandra and XMM‑Newton counterparts within a given
matching radius.  When more than one counterpart is found within the
tolerance the identifiers are concatenated into a comma‑separated string.

The resulting catalogue is saved as `ngc2264_combined.csv` in the
current working directory.

Requirements:
  - astroquery
  - astropy
  - pandas

Usage:
    python cross_match_ngc2264.py

This file is auto‑generated in the sandbox environment for testing.
It mirrors the content provided by the user and is used to validate the
script behaviour and propose improvements.
"""

from __future__ import annotations

import sys
from typing import List, Optional
import numpy as np
import numpy.linalg as la

import os as _os
import matplotlib as _mpl
# Force a non-interactive backend when used in servers/threads unless explicitly overridden
try:
    _bk = str(_os.environ.get('MPLBACKEND', '')).lower()
    if _bk not in ('agg', 'module://matplotlib_inline.backend_inline'):
        _mpl.use('Agg', force=True)
except Exception:
    pass
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.transforms import offset_copy

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, AsinhStretch, PercentileInterval
from astroquery.skyview import SkyView
from astropy.utils.data import download_file
import io
import requests
import time

try:
    from astroquery.gaia import Gaia
    # Remove default row limit to retrieve all sources
    Gaia.ROW_LIMIT = -1
    from astroquery.vizier import Vizier
    # Note: HEASARC queries were used in earlier versions of this script to
    # obtain Chandra and XMM‑Newton sources. In order to incorporate the
    # updated X‑ray catalogues from Flaccomio et al. (2023), we now query
    # the relevant tables in VizieR (J/A+A/670/A37/tablea1 and table2) via
    # astroquery.vizier. The HEASARC import is therefore no longer needed.
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import pandas as pd
    # Optional SAMP support (Aladin Desktop)
    try:
        from astropy.samp import SAMPIntegratedClient  # type: ignore
    except Exception:
        SAMPIntegratedClient = None
    # No IDL .sav reading: we rely on a pre-exported CSV for Chandra
except ImportError as exc:
    print("Missing required packages: {}".format(exc))
    print("Install dependencies with: pip install astropy astroquery pandas")
    sys.exit(1)


import os
# Directory to cache downloaded catalogs

CACHE_DIR = "catalogs"
os.makedirs(CACHE_DIR, exist_ok=True)

# Subdirectory for cached sky images (e.g., 2MASS J cutouts)
IMAGES_DIR = os.path.join(CACHE_DIR, "2MASS_J")
os.makedirs(IMAGES_DIR, exist_ok=True)

import argparse


def query_gaia(center: SkyCoord, radius: u.Quantity) -> pd.DataFrame:
    """Query Gaia DR3 around the given sky position.

    Parameters
    ----------
    center : SkyCoord
        Centre of the cone search.
    radius : Quantity
        Search radius.

    Returns
    -------
    DataFrame
        Table of Gaia sources with RA, Dec and source_id.
    """
    cache_file = os.path.join(
        CACHE_DIR,
        f"gaia_RA{center.ra.deg:.4f}_DEC{center.dec.deg:.4f}_R{radius.value:.4f}.csv"
    )
    # Resolve the global REFRESH flag at runtime.  If REFRESH has not been
    # defined yet, default to ``False`` so that the cache is used if present.
    refresh = globals().get("REFRESH", False)
    if not refresh and os.path.exists(cache_file):
        # Read Gaia cache, ensuring gaia_id is int64
        df = pd.read_csv(cache_file, dtype={'gaia_id': 'int64'})
        # Backward compatibility: older caches may lack PM columns; add safe defaults
        for col, default in (
            ('pmra', 0.0), ('pmdec', 0.0),
            ('pmra_error', 0.0), ('pmdec_error', 0.0),
            ('ref_epoch', 2016.0),
        ):
            if col not in df.columns:
                df[col] = default
        return df
    # Use ADQL query to explicitly fetch ra_dec_corr
    ra_center = center.ra.deg
    dec_center = center.dec.deg
    rad_deg = radius.to(u.deg).value
    adql = f"""
    SELECT ra, dec, ra_error, dec_error, ra_dec_corr,
           pmra, pmdec, pmra_error, pmdec_error, ref_epoch,
           source_id, phot_g_mean_mag, ruwe
      FROM gaiadr3.gaia_source
     WHERE CONTAINS(
             POINT('ICRS', ra, dec),
             CIRCLE('ICRS', {ra_center}, {dec_center}, {rad_deg})
           )=1
    """
    job = Gaia.launch_job_async(adql)
    results = job.get_results()
    # Convert to DataFrame and keep relevant columns, including ra_dec_corr
    df = results.to_pandas()[[
        'ra', 'dec',
        'ra_error', 'dec_error', 'ra_dec_corr',
        'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'ref_epoch',
        'source_id', 'phot_g_mean_mag', 'ruwe'
    ]]
    # Rename and cast Gaia ID to integer
    df = df.rename(columns={
        'ra': 'ra_deg',
        'dec': 'dec_deg',
        'ra_error': 'ra_err_mas',
        'dec_error': 'dec_err_mas',
        'ra_dec_corr': 'ra_dec_corr',
        'pmra': 'pmra',
        'pmdec': 'pmdec',
        'pmra_error': 'pmra_error',
        'pmdec_error': 'pmdec_error',
        'ref_epoch': 'ref_epoch',
        'source_id': 'gaia_id',
        'phot_g_mean_mag': 'Gmag',
        'ruwe': 'ruwe'
    })
    df['gaia_id'] = df['gaia_id'].astype('int64')
    # Convert milliarcsecond errors to degrees
    df['ra_err_deg'] = df['ra_err_mas'] * u.mas.to(u.deg)
    df['dec_err_deg'] = df['dec_err_mas'] * u.mas.to(u.deg)
    # Compute full error ellipse (major, minor, PA) from RA/Dec errors and covariance
    # Correct RA error for sky projection
    ra_err_sky = df['ra_err_deg'] * np.cos(np.deg2rad(df['dec_deg']))
    dec_err    = df['dec_err_deg']
    cov_rd     = df['ra_dec_corr'] * ra_err_sky * dec_err

    a = ra_err_sky**2
    b = dec_err**2
    mean = 0.5 * (a + b)
    diff = 0.5 * (a - b)
    sqrt_term = np.sqrt(diff**2 + cov_rd**2)

    lam1 = mean + sqrt_term   # major-axis variance
    lam2 = mean - sqrt_term   # minor-axis variance

    # Convert variances to 1σ radii in arcsec
    df['errMaj'] = np.sqrt(lam1) * 3600.0
    df['errMin'] = np.sqrt(lam2) * 3600.0

    # Position angle east of north (0°=north)
    df['errPA'] = (np.degrees(np.arctan2(cov_rd, lam1 - a)) % 360.0)
    # Drop per-axis uncertainty columns (_mas converted to _deg), raw mas columns, and ra_dec_corr
    drop_cols = [
        'ra_err_mas', 'dec_err_mas',
        'ra_err_deg', 'dec_err_deg',
        'ra_dec_corr'
    ]
    for _c in list(drop_cols):
        if _c in df.columns:
            df = df.drop(columns=[_c])
    df.to_csv(cache_file, index=False)
    return df


def query_2mass(center: SkyCoord, radius: u.Quantity) -> pd.DataFrame:
    """Query 2MASS Point Source Catalogue via Vizier.

    Returns a DataFrame with positions and photometry.
    """
    cache_file = os.path.join(
        CACHE_DIR,
        f"2MASS_RA{center.ra.deg:.4f}_DEC{center.dec.deg:.4f}_R{radius.value:.4f}.csv"
    )
    # Resolve REFRESH at runtime; use False if undefined
    refresh = globals().get("REFRESH", False)
    if not refresh and os.path.exists(cache_file):
        return pd.read_csv(cache_file)
    # Request a broad set of columns including quality flags so that
    # we can cache the complete table and filter later.
    v = Vizier(
        columns=[
            'RAJ2000', 'DEJ2000',
            'errMaj', 'errMin', 'errPA',
            'Jmag', 'Hmag', 'Kmag', '2MASS',
            # Common 2MASS quality/contamination flags (if available)
            'Qflg', 'Cflg', 'Xflg', 'prox', 'bl_flg'
        ],
        catalog='II/246', row_limit=-1
    )
    tables = v.query_region(center, radius=radius)
    if len(tables) == 0:
        return pd.DataFrame(columns=[
            'ra_deg', 'dec_deg',
            'pos_err_maj_arcsec', 'pos_err_min_arcsec', 'pos_err_pa_deg',
            'Jmag', 'Hmag', 'Kmag', '2MASS'
        ])
    tab = tables[0]
    df = tab.to_pandas()
    df = df.rename(columns={
        'RAJ2000': 'ra_deg',
        'DEJ2000': 'dec_deg',
        'errMaj': 'pos_err_maj_arcsec',
        'errMin': 'pos_err_min_arcsec',
        'errPA': 'pos_err_pa_deg'
    })
    # Define uncertainty ellipse parameters
    df['errMaj'] = df['pos_err_maj_arcsec']
    df['errMin'] = df['pos_err_min_arcsec']
    df['errPA'] = df['pos_err_pa_deg']
    # Derive RA/Dec uncertainties from ellipse axes (arcsec → deg)
    df['ra_err_deg']  = df['pos_err_min_arcsec'] / 3600.0
    df['dec_err_deg'] = df['pos_err_maj_arcsec'] / 3600.0
    # Drop per-axis uncertainty columns and raw positional errors
    for _c in ['ra_err_deg', 'dec_err_deg', 'pos_err_maj_arcsec', 'pos_err_min_arcsec', 'pos_err_pa_deg']:
        if _c in df.columns:
            df = df.drop(columns=[_c])
    df.to_csv(cache_file, index=False)
    # Return the full dataframe (with added errMaj/errMin/errPA) so filters can use flags
    return df


def query_wise(center: SkyCoord, radius: u.Quantity) -> pd.DataFrame:
    """Query the AllWISE catalogue via VizieR.

    Returns a DataFrame with positions, photometry, and designation.
    """
    cache_file = os.path.join(
        CACHE_DIR,
        f"WISE_RA{center.ra.deg:.4f}_DEC{center.dec.deg:.4f}_R{radius.value:.4f}.csv"
    )
    # Resolve REFRESH at runtime; use False if undefined
    refresh = globals().get("REFRESH", False)
    if not refresh and os.path.exists(cache_file):
        return pd.read_csv(cache_file)
    # Use Vizier to query the AllWISE PSC (II/328/allwise)
    v = Vizier(
        columns=[
            'RAJ2000', 'DEJ2000',
            'W1mpro', 'W2mpro', 'W3mpro', 'W4mpro',
            'AllWISE', 'eeMaj', 'eeMin', 'eePA',
            # Quality flags commonly used in AllWISE
            'ph_qual', 'ccf', 'ext_flg', 'w1snr', 'w2snr'
        ],
        catalog='II/328/allwise', row_limit=-1
    )
    tables = v.query_region(center, radius=radius)
    if len(tables) == 0:
        return pd.DataFrame(columns=[
            'ra_deg', 'dec_deg', 'W1mpro', 'W2mpro', 'W3mpro', 'W4mpro',
            'wise_id', 'pos_err_maj_arcsec', 'pos_err_min_arcsec', 'pos_err_pa_deg'
        ])
    tab = tables[0]
    df = tab.to_pandas()
    df.rename(columns={
        'RAJ2000': 'ra_deg',
        'DEJ2000': 'dec_deg',
        'AllWISE': 'wise_id',
        'eeMaj':   'pos_err_maj_arcsec',
        'eeMin':   'pos_err_min_arcsec',
        'eePA':    'pos_err_pa_deg'
    }, inplace=True)
    # Define uncertainty ellipse parameters
    df['errMaj'] = df['pos_err_maj_arcsec']
    df['errMin'] = df['pos_err_min_arcsec']
    df['errPA']  = df['pos_err_pa_deg']
    # Derive RA/Dec uncertainties from ellipse axes (arcsec → deg)
    df['ra_err_deg']  = df['pos_err_min_arcsec'] / 3600.0
    df['dec_err_deg'] = df['pos_err_maj_arcsec'] / 3600.0
    # Identify WISE photometry columns (case-insensitive)
    phot_fields = ['W1MPRO', 'W2MPRO', 'W3MPRO', 'W4MPRO']
    phot_cols = [c for c in df.columns if c.upper() in phot_fields]
    # Sort photometry columns in order
    phot_cols_sorted = [
        next((c for c in phot_cols if c.upper() == f), None)
        for f in phot_fields
    ]
    final_cols = (
        ['ra_deg', 'dec_deg']
        + [c for c in phot_cols_sorted if c]
        + ['wise_id', 'errMaj', 'errMin', 'errPA']
    )
    # Drop per-axis uncertainty and raw positional error columns (if present)
    for _c in ['ra_err_deg', 'dec_err_deg', 'pos_err_maj_arcsec', 'pos_err_min_arcsec', 'pos_err_pa_deg']:
        if _c in df.columns:
            df = df.drop(columns=[_c])
    df.to_csv(cache_file, index=False)
    # Return full df with added errMaj/errMin/errPA and quality columns
    return df


def query_chandra(center: SkyCoord, radius: u.Quantity) -> pd.DataFrame:
    """Query the Flaccomio et al. (2023) Chandra source list via VizieR.

    The catalogue `J/A+A/670/A37/tablea1` contains 1043 Chandra/ACIS
    detections in the NGC 2264 field with positions in columns RAJ2000 and
    DEJ2000 and a source identifier in the `ACIS` column【355715869545438†L106-L112】.  This function
    performs a cone search on that table around ``center`` with search
    ``radius`` and returns a DataFrame of decimal‑degree coordinates along
    with the ACIS identifier.  If no sources are found, an empty DataFrame
    with the expected columns is returned.
    """
    cache_file = os.path.join(
        CACHE_DIR,
        f"Chandra_RA{center.ra.deg:.4f}_DEC{center.dec.deg:.4f}_R{radius.value:.4f}.csv"
    )
    # Resolve REFRESH at runtime; use False if undefined
    refresh = globals().get("REFRESH", False)
    if not refresh and os.path.exists(cache_file):
        return pd.read_csv(cache_file)
    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'ACIS', 'epos'], catalog='J/A+A/670/A37/tablea1', row_limit=-1)
    tables = v.query_region(center, radius=radius)
    if len(tables) == 0:
        return pd.DataFrame(columns=['ra_deg', 'dec_deg', 'errMaj', 'errMin', 'errPA', 'chandra_id'])
    tab = tables[0].to_pandas()
    # Parse RAJ2000/DEJ2000 strings into decimal degrees
    coords = SkyCoord(tab['RAJ2000'], tab['DEJ2000'], unit=(u.hourangle, u.deg), frame='icrs')
    tab['ra_deg'] = coords.ra.deg
    tab['dec_deg'] = coords.dec.deg
    tab = tab.rename(columns={
        'ACIS': 'chandra_id',
        'epos': 'pos_err_arcsec'
    })
    # Fallback: if rename didn't work, find any 'epos' column ignoring case
    if 'pos_err_arcsec' not in tab.columns:
        for col in tab.columns:
            if col.lower() == 'epos':
                tab['pos_err_arcsec'] = tab[col]
                break
    # Use symmetric positional error for ellipse
    tab['errMaj'] = tab['pos_err_arcsec']
    tab['errMin'] = tab['pos_err_arcsec']
    tab['errPA'] = 0.0
    # Derive RA/Dec uncertainties (symmetric error, arcsec → deg)
    tab['ra_err_deg']  = tab['pos_err_arcsec'] / 3600.0
    tab['dec_err_deg'] = tab['pos_err_arcsec'] / 3600.0
    # Drop per-axis uncertainty and raw positional error columns (if present)
    for _c in ['ra_err_deg', 'dec_err_deg', 'pos_err_arcsec']:
        if _c in tab.columns:
            tab = tab.drop(columns=[_c])
    tab.to_csv(cache_file, index=False)
    return tab[[
        'ra_deg', 'dec_deg',
        'errMaj', 'errMin', 'errPA',
        'chandra_id'
    ]]




def query_xmm(center: SkyCoord, radius: u.Quantity) -> pd.DataFrame:
    """Query the Flaccomio et al. (2023) XMM‑Newton source list via VizieR.

    The catalogue `J/A+A/670/A37/table2` provides 944 XMM‑Newton detections
    with positions in columns RAJ2000 and DEJ2000 and a source
    identification number in the `XMM` column【583771511037106†L107-L120】.  This function
    uses astroquery.vizier to perform a cone search and returns the
    results with uniform column names.  If the query yields no sources,
    an empty DataFrame with the expected columns is returned.
    """
    cache_file = os.path.join(
        CACHE_DIR,
        f"XMM_RA{center.ra.deg:.4f}_DEC{center.dec.deg:.4f}_R{radius.value:.4f}.csv"
    )
    # Resolve REFRESH at runtime; use False if undefined
    refresh = globals().get("REFRESH", False)
    if not refresh and os.path.exists(cache_file):
        return pd.read_csv(cache_file)
    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'XMM', 'Id.rad'], catalog='J/A+A/670/A37/table2', row_limit=-1)
    tables = v.query_region(center, radius=radius)
    if len(tables) == 0:
        return pd.DataFrame(columns=['ra_deg', 'dec_deg', 'errMaj', 'errMin', 'errPA', 'xmm_id'])
    tab = tables[0].to_pandas()
    # Parse RAJ2000/DEJ2000 strings into decimal degrees
    coords = SkyCoord(tab['RAJ2000'], tab['DEJ2000'], unit=(u.hourangle, u.deg), frame='icrs')
    tab['ra_deg'] = coords.ra.deg
    tab['dec_deg'] = coords.dec.deg
    tab = tab.rename(columns={
        'XMM': 'xmm_id',
        'Id.rad': 'pos_err_arcsec'
    })
    # Use symmetric positional error for ellipse
    tab['errMaj'] = tab['pos_err_arcsec']
    tab['errMin'] = tab['pos_err_arcsec']
    tab['errPA'] = 0.0
    # Derive RA/Dec uncertainties (symmetric error, arcsec → deg)
    tab['ra_err_deg']  = tab['pos_err_arcsec'] / 3600.0
    tab['dec_err_deg'] = tab['pos_err_arcsec'] / 3600.0
    # Drop per-axis uncertainty and raw positional error columns (if present)
    for _c in ['ra_err_deg', 'dec_err_deg', 'pos_err_arcsec']:
        if _c in tab.columns:
            tab = tab.drop(columns=[_c])
    tab.to_csv(cache_file, index=False)
    return tab[[
        'ra_deg', 'dec_deg',
        'errMaj', 'errMin', 'errPA',
        'xmm_id'
    ]]


# -------------------- Quality filters per catalog --------------------
def _print_filter_counts(name: str, total: int, kept: int) -> None:
    dropped = total - kept
    frac = (kept / total * 100.0) if total > 0 else 0.0
    print(f"{name}: total={total}, after quality filters={kept} ({frac:.1f}% kept, {dropped} dropped)")


def filter_gaia_quality(df: pd.DataFrame) -> pd.DataFrame:
    total = len(df)
    if total == 0:
        _print_filter_counts('Gaia', 0, 0)
        return df
    mask = pd.Series(True, index=df.index)
    if 'ruwe' in df.columns:
        # Common astrometric quality cut
        with np.errstate(invalid='ignore'):
            mask &= (df['ruwe'] <= 1.4)
    kept = int(mask.sum())
    _print_filter_counts('Gaia', total, kept)
    return df.loc[mask].reset_index(drop=True)


def filter_2mass_quality(df: pd.DataFrame) -> pd.DataFrame:
    total = len(df)
    if total == 0:
        _print_filter_counts('2MASS', 0, 0)
        return df
    mask = pd.Series(True, index=df.index)
    # Photometric quality: require at least one of J/H/K with quality in {A,B,C}
    qual_series = None
    for col in ['ph_qual', 'Qflg', 'qflg']:
        if col in df.columns:
            qual_series = df[col].astype(str)
            break
    if qual_series is not None:
        def good_any_abc(s: str) -> bool:
            s = s.strip()
            if len(s) == 0:
                return False
            # consider first three chars as J/H/K if present
            candidates = list(s[:3])
            return any(ch in ('A', 'B', 'C') for ch in candidates)
        mask &= qual_series.apply(good_any_abc)
    # Contamination flags: require zeros if available
    for col in ['Cflg', 'cc_flg', 'ccflag']:
        if col in df.columns:
            s = df[col].astype(str)
            mask &= s.apply(lambda x: all((c == '0') for c in list(x[:3])))
            break
    # Extended flag: prefer point sources
    for col in ['Xflg', 'ext_flg']:
        if col in df.columns:
            with np.errstate(invalid='ignore'):
                mask &= (pd.to_numeric(df[col], errors='coerce').fillna(0) == 0)
            break
    kept = int(mask.sum())
    _print_filter_counts('2MASS', total, kept)
    return df.loc[mask].reset_index(drop=True)


def filter_wise_quality(df: pd.DataFrame) -> pd.DataFrame:
    total = len(df)
    if total == 0:
        _print_filter_counts('WISE', 0, 0)
        return df
    mask = pd.Series(True, index=df.index)
    # Photometric quality: W1/W2 not upper limits
    if 'ph_qual' in df.columns:
        pq = df['ph_qual'].astype(str)
        def ok_ph(s: str) -> bool:
            s = s.strip()
            if len(s) == 0:
                return True
            w1 = s[0] if len(s) >= 1 else 'U'
            w2 = s[1] if len(s) >= 2 else 'U'
            return (w1 in ('A','B','C')) and (w2 in ('A','B','C'))
        mask &= pq.apply(ok_ph)
    # Artifact flags: require '0' in W1/W2 if ccf is present
    for col in ['ccf', 'cc_flags', 'ccflag']:
        if col in df.columns:
            cf = df[col].astype(str)
            def ok_cc(s: str) -> bool:
                s = s.strip()
                if len(s) == 0:
                    return True
                w1 = s[0] if len(s) >= 1 else '0'
                w2 = s[1] if len(s) >= 2 else '0'
                return (w1 == '0') and (w2 == '0')
            mask &= cf.apply(ok_cc)
            break
    # Extended flag: prefer point sources
    if 'ext_flg' in df.columns:
        with np.errstate(invalid='ignore'):
            mask &= (pd.to_numeric(df['ext_flg'], errors='coerce').fillna(0) == 0)
    kept = int(mask.sum())
    _print_filter_counts('WISE', total, kept)
    return df.loc[mask].reset_index(drop=True)


def filter_chandra_quality(df: pd.DataFrame) -> pd.DataFrame:
    total = len(df)
    if total == 0:
        _print_filter_counts('Chandra', 0, 0)
        return df
    mask = pd.Series(True, index=df.index)
    # Minimal sanity filters: finite coordinates and positive positional error
    with np.errstate(invalid='ignore'):
        mask &= df['ra_deg'].notna() & df['dec_deg'].notna()
        if 'errMaj' in df.columns:
            mask &= (pd.to_numeric(df['errMaj'], errors='coerce') > 0)
    kept = int(mask.sum())
    _print_filter_counts('Chandra', total, kept)
    return df.loc[mask].reset_index(drop=True)


def filter_xmm_quality(df: pd.DataFrame) -> pd.DataFrame:
    total = len(df)
    if total == 0:
        _print_filter_counts('XMM', 0, 0)
        return df
    mask = pd.Series(True, index=df.index)
    # Minimal sanity filters: finite coordinates and positive positional error
    with np.errstate(invalid='ignore'):
        mask &= df['ra_deg'].notna() & df['dec_deg'].notna()
        if 'errMaj' in df.columns:
            mask &= (pd.to_numeric(df['errMaj'], errors='coerce') > 0)
    kept = int(mask.sum())
    _print_filter_counts('XMM', total, kept)
    return df.loc[mask].reset_index(drop=True)


def query_chandra_csv(center: SkyCoord, radius: u.Quantity, csv_path: str) -> pd.DataFrame:
    """Load Chandra sources from a pre-exported CSV and filter to cone.

    The CSV is expected to provide at least: ra_deg, dec_deg, and an ID
    column (pub_n or chandra_id). Positional error column (epos or
    errMaj/errMin) is optional; defaults to 1.0 arcsec if missing.
    """
    if not csv_path:
        raise ValueError("CSV path is empty; use --chandra-csv-path to provide a file.")
    df_raw = pd.read_csv(csv_path)
    # Normalize columns
    def pick(*names):
        for n in names:
            if n in df_raw.columns:
                return n
        return None
    ra_col = pick('ra_deg','RA_deg','RA')
    dec_col = pick('dec_deg','DEC_deg','Dec')
    if ra_col is None or dec_col is None:
        raise KeyError(f"CSV must contain 'ra_deg' and 'dec_deg' columns; has {list(df_raw.columns)}")
    id_col = pick('chandra_id','pub_n','id','srcid','num')
    if id_col is None:
        raise KeyError("CSV must contain an identifier column (e.g., 'pub_n' or 'chandra_id').")
    err_col = pick('epos','errMaj','err_arcsec')

    ra = pd.to_numeric(df_raw[ra_col], errors='coerce')
    dec = pd.to_numeric(df_raw[dec_col], errors='coerce')
    cid = df_raw[id_col]
    if err_col is not None:
        err = pd.to_numeric(df_raw[err_col], errors='coerce').fillna(1.0)
    else:
        err = pd.Series(1.0, index=df_raw.index)
    df_full = pd.DataFrame({
        'ra_deg': ra,
        'dec_deg': dec,
        'errMaj': err,
        'errMin': err,
        'errPA':  0.0,
        'chandra_id': cid
    }).dropna(subset=['ra_deg','dec_deg'])

    # Filter spatially
    coords = SkyCoord(ra=df_full['ra_deg'].values * u.deg,
                      dec=df_full['dec_deg'].values * u.deg, frame='icrs')
    sep = coords.separation(center)
    mask = sep <= radius
    df = df_full.loc[mask].reset_index(drop=True)
    # Write cache consistent with other queries
    cache_file = os.path.join(
        CACHE_DIR,
        f"Chandra_RA{center.ra.deg:.4f}_DEC{center.dec.deg:.4f}_R{radius.value:.4f}.csv"
    )
    df.to_csv(cache_file, index=False)
    return df


# -------------------- Blend-aware association for 2MASS --------------------
def attach_2mass_blends(combined_df: pd.DataFrame,
                        tmass_df: pd.DataFrame,
                        *,
                        r_group_arcsec: float = 1.5,
                        max_pair_sep_arcsec: float = 2.0,
                        r_centroid_arcsec: float = 0.5,
                        max_dG_mag: float = 1.0,
                        prefer_blflag: bool = True) -> pd.DataFrame:
    """Attach a single 2MASS source to a compact group of Gaia stars (blend).

    The normal D²≤1 matching remains untouched. This pass only handles
    2MASS sources that remain unmatched after the standard merge, by
    identifying nearby compact groups of Gaia stars consistent with an
    unresolved blend in 2MASS.

    Criteria (low‑spurious):
      - at least two Gaia within ``r_group_arcsec`` of the 2MASS pos
      - pairwise Gaia separations ≤ ``max_pair_sep_arcsec``
      - 2MASS vs Gaia centroid distance ≤ ``r_centroid_arcsec``
      - optional brightness similarity: ΔG ≤ ``max_dG_mag``
      - if ``prefer_blflag`` and 2MASS has ``bl_flg>0``, treat as strong evidence
    """
    if combined_df.empty or tmass_df.empty or 'gaia_id' not in combined_df.columns:
        return combined_df

    # Only consider Gaia master rows (exclude synthetic rows)
    is_gaia = (~combined_df['gaia_id'].isna()) & (combined_df['gaia_id'] != -1)
    if not np.any(is_gaia):
        return combined_df
    gaia = combined_df.loc[is_gaia, ['ra_deg', 'dec_deg', 'Gmag', 'gaia_id']].copy()
    gaia.rename(columns={'ra_deg': 'ra', 'dec_deg': 'dec'}, inplace=True)

    # Consider as already-attached only 2MASS IDs that are linked to at least
    # one VALID Gaia master row (gaia_id != -1). Synthetic rows should not block.
    if '2MASS' in combined_df.columns:
        attached_ids = set(
            str(x) for x in combined_df.loc[is_gaia, '2MASS'].dropna().astype(str).values
        )
    else:
        attached_ids = set()

    def _deg_to_arcsec(dx_deg, dy_deg, dec0_deg):
        return (dx_deg * np.cos(np.deg2rad(dec0_deg)) * 3600.0,
                dy_deg * 3600.0)

    for _, tm in tmass_df.iterrows():
        tm_id = str(tm.get('2MASS'))
        if tm_id in attached_ids:
            continue
        ra0 = float(tm['ra_deg']); dec0 = float(tm['dec_deg'])

        # Geometry around this 2MASS source
        dx_as, dy_as = _deg_to_arcsec(gaia['ra'].values - ra0, gaia['dec'].values - dec0, dec0)
        sep_as = np.hypot(dx_as, dy_as)
        cand = np.where(sep_as <= r_group_arcsec)[0]
        if len(cand) < 2:
            continue

        # Pairwise compactness
        ok_pairs = True
        idxs = cand.tolist()
        for i in range(len(idxs)):
            for j in range(i + 1, len(idxs)):
                a = gaia.iloc[idxs[i]]; b = gaia.iloc[idxs[j]]
                dxa, dya = _deg_to_arcsec(float(a['ra'] - b['ra']), float(a['dec'] - b['dec']), dec0)
                if np.hypot(dxa, dya) > max_pair_sep_arcsec:
                    ok_pairs = False
                    break
            if not ok_pairs:
                break
        if not ok_pairs:
            continue

        # Centroid proximity
        ra_c = float(np.mean(gaia.iloc[idxs]['ra'].values))
        dec_c = float(np.mean(gaia.iloc[idxs]['dec'].values))
        dx_c_as, dy_c_as = _deg_to_arcsec(ra_c - ra0, dec_c - dec0, dec0)
        if np.hypot(dx_c_as, dy_c_as) > r_centroid_arcsec:
            continue

        # Brightness similarity (optional)
        gvals = gaia.iloc[idxs]['Gmag'].values
        if np.nanmax(gvals) - np.nanmin(gvals) > max_dG_mag:
            continue

        # Prefer 2MASS blend flag when available
        if prefer_blflag and ('bl_flg' in tmass_df.columns):
            try:
                blv = int(tm.get('bl_flg')) if not pd.isna(tm.get('bl_flg')) else 0
            except Exception:
                blv = 0
            # If explicitly non-blended and we already have strict geometry, we can still attach; else proceed
            # No hard rejection here; the geometry acts as the main guardrail.
            pass

        # Attach to all Gaia members in the group
        gaia_ids = [int(g) for g in gaia.iloc[idxs]['gaia_id'].values]
        for gid in gaia_ids:
            mask = combined_df['gaia_id'] == gid
            combined_df.loc[mask, '2MASS'] = tm.get('2MASS')
            combined_df.loc[mask, 'blend_2mass'] = True
            combined_df.loc[mask, 'blend_n_2mass'] = len(gaia_ids)
            combined_df.loc[mask, 'blend_members_2mass'] = ';'.join(str(x) for x in gaia_ids)
        attached_ids.add(tm_id)

    return combined_df


def _dedupe_unmatched_for_id(combined_df: pd.DataFrame, id_col: str) -> pd.DataFrame:
    """Remove synthetic rows (gaia_id == -1) for a given catalog id when
    the same id is already attached to at least one valid Gaia row.

    This prevents duplicate entries like one row with GAIA='—' and another
    with a GAIA match for the same catalog identifier (e.g., 2MASS).
    """
    if id_col not in combined_df.columns or 'gaia_id' not in combined_df.columns:
        return combined_df
    try:
        mask_valid = (~combined_df['gaia_id'].isna()) & (combined_df['gaia_id'] != -1)
        ids_with_valid = set(combined_df.loc[mask_valid, id_col].dropna().astype(str).values)
        drop_mask = (~mask_valid) & (combined_df[id_col].astype(str).isin(ids_with_valid))
        n_before = len(combined_df)
        if drop_mask.any():
            combined_df = combined_df.loc[~drop_mask].reset_index(drop=True)
            n_removed = n_before - len(combined_df)
            if n_removed > 0:
                print(f"[dedupe] removed {n_removed} synthetic rows for {id_col}")
    except Exception:
        pass
    return combined_df


def cross_match(source_df: pd.DataFrame,
                other_df: pd.DataFrame,
                other_id_col: str,
                min_radius: u.Quantity,
                pdf_path: Optional[str] = None,
                scale_factor: float = 2.0,
                all_catalogs: Optional[dict] = None,
                plot_mode: str = 'match') -> List[List]:
    """Cross‑match two tables based on positions.

    Parameters
    ----------
    source_df : DataFrame
        DataFrame containing 'ra_deg' and 'dec_deg' columns; used as the
        reference table (Gaia).
    other_df : DataFrame
        Secondary DataFrame also containing 'ra_deg' and 'dec_deg'.
    other_id_col : str
        Column name in `other_df` whose value(s) should be attached to the
        reference table when a match occurs.
    min_radius : u.Quantity
        Minimum matching radius (in arcsec) enforced on the semi‑axes of
        the positional uncertainty ellipses.  Ellipses smaller than this
        radius (converted to degrees) will be clipped to `min_radius`.
        Matching uses a Mahalanobis distance threshold of 1.0 (i.e.,
        a 1‑σ ellipse); if a more permissive match is desired, increase
        the ``scale_factor`` accordingly.
    pdf_path : str, optional
        If provided, write diagnostic plots of matching ellipses to this
        PDF file.
    scale_factor : float, optional
        Factor to scale the positional uncertainty ellipses by before
        performing the Mahalanobis distance test.
    all_catalogs : dict, optional
        Dictionary of catalog metadata and dataframes for plotting.
    plot_mode : {'native','match'}
        Controls which ellipses are plotted in the PDF. If 'native', the
        plots show per-catalog positional uncertainties as provided by the
        archives. If 'match', the plots show the ellipses actually used for
        identification (i.e., errors scaled by `scale_factor` and clipped by
        `min_radius`). Default is 'match'.

    Returns
    -------
    list of lists
        Each element is a list of identifiers (strings) matched to the
        corresponding row in `source_df`.  Empty lists indicate no match.
    """
    if other_df.empty:
        return [[] for _ in range(len(source_df))]

    pdf = PdfPages(pdf_path) if pdf_path else None

    # The SkyCoord constructions below were unused; matching is done purely
    # via small–angle approximations in Cartesian space.  They have been
    # removed to avoid unnecessary overhead.

    # Helper to extract per-row target epoch (Julian year) for the secondary catalog.
    def _epoch_of_row(j: int) -> Optional[float]:
        meta = None
        if isinstance(all_catalogs, dict):
            meta = all_catalogs.get(other_id_col.replace('_id', ''), {})
        # Candidate numeric columns
        candidates = [
            'ref_epoch', 'epoch', 'EPOCH', 'Epoch',
            'MJD', 'mjd', 'JD', 'jd',
            'w1mjdmean', 'w2mjdmean', 'w3mjdmean', 'w4mjdmean'
        ]
        row = other_df.iloc[j]
        # Try direct year first
        for c in ['ref_epoch', 'epoch', 'EPOCH', 'Epoch']:
            if c in other_df.columns:
                try:
                    val = float(row[c])
                    if not np.isnan(val):
                        return val
                except Exception:
                    pass
        # Try MJD/averaged WISE MJDs
        mjds = []
        for c in ['MJD', 'mjd', 'w1mjdmean', 'w2mjdmean', 'w3mjdmean', 'w4mjdmean']:
            if c in other_df.columns:
                try:
                    v = float(row[c])
                    if not np.isnan(v) and v > 0:
                        mjds.append(v)
                except Exception:
                    pass
        if mjds:
            mjd = float(np.nanmean(mjds))
            # Convert MJD to Julian year (approx)
            return 2000.0 + (mjd - 51544.5) / 365.25
        # Try JD
        for c in ['JD', 'jd']:
            if c in other_df.columns:
                try:
                    jd = float(row[c])
                    if not np.isnan(jd) and jd > 0:
                        return 2000.0 + (jd - 2451545.0) / 365.25
                except Exception:
                    pass
        # Try ISO date strings
        for c in ['date', 'dateobs', 'obs_date']:
            if c in other_df.columns:
                try:
                    from datetime import datetime
                    s = str(row[c])
                    dt = datetime.fromisoformat(s)
                    year_start = datetime(dt.year, 1, 1)
                    doy = (dt - year_start).days + 0.5
                    return dt.year + doy / 365.25
                except Exception:
                    pass
        # Fallback to catalog-level epoch if provided
        try:
            return float(meta.get('epoch')) if meta and (meta.get('epoch') is not None) else None
        except Exception:
            return None

    # Precompute 2x2 positional covariance matrices for source objects (no PM inflation here)
    # Apply scale factor to ellipse axes before converting to degrees
    min_deg = min_radius.to(u.deg).value
    src_maj = np.maximum((source_df['errMaj'] * scale_factor).values / 3600.0, min_deg)  # deg
    src_min = np.maximum((source_df['errMin'] * scale_factor).values / 3600.0, min_deg)  # deg
    src_pa  = np.deg2rad(source_df['errPA'].values)  # radians
    n_src = len(source_df)
    src_cov = np.zeros((n_src, 2, 2))
    # Gather per-row proper motion and uncertainties if available
    has_pm = ('pmra' in source_df.columns) and ('pmdec' in source_df.columns)
    has_pm_err = ('pmra_error' in source_df.columns) and ('pmdec_error' in source_df.columns)
    has_epoch = ('ref_epoch' in source_df.columns)
    for i in range(n_src):
        a = src_maj[i]
        b = src_min[i]
        phi = src_pa[i]
        cphi = np.cos(phi)
        sphi = np.sin(phi)
        # Major axis unit vector degrees: (sin(PA), cos(PA))
        vmaj = np.array([sphi, cphi])
        # Minor axis unit vector: (-cos(PA), sin(PA))
        vmin = np.array([-cphi, sphi])
        cov = a*a * np.outer(vmaj, vmaj) + b*b * np.outer(vmin, vmin)
        src_cov[i] = cov

    # Precompute 2x2 covariance matrices for other objects
    oth_maj = np.maximum((other_df['errMaj'] * scale_factor).values / 3600.0, min_deg)  # deg
    oth_min = np.maximum((other_df['errMin'] * scale_factor).values / 3600.0, min_deg)  # deg
    oth_pa  = np.deg2rad(other_df['errPA'].values)  # radians
    n_oth = len(other_df)
    oth_cov = np.zeros((n_oth, 2, 2))
    for j in range(n_oth):
        a = oth_maj[j]
        b = oth_min[j]
        phi = oth_pa[j]
        cphi = np.cos(phi)
        sphi = np.sin(phi)
        vmaj = np.array([sphi, cphi])
        vmin = np.array([-cphi, sphi])
        oth_cov[j] = a*a * np.outer(vmaj, vmaj) + b*b * np.outer(vmin, vmin)

    # Precompute fast bounds and arrays for candidate pruning
    other_ra = other_df['ra_deg'].values
    other_dec = other_df['dec_deg'].values
    r_src_bound = np.maximum(src_maj, src_min)  # deg
    r_oth_bound = np.maximum(oth_maj, oth_min)  # deg
    r_oth_max = float(r_oth_bound.max()) if len(r_oth_bound) else 0.0
    # Precompute per-row epochs for other_df once
    epoch_list = []
    for _j in range(len(other_df)):
        try:
            epoch_list.append(_epoch_of_row(_j))
        except Exception:
            epoch_list.append(None)
    # Epoch span across the secondary catalog
    if any(e is not None for e in epoch_list):
        _e_arr = np.array([e if e is not None else np.nan for e in epoch_list], dtype=float)
        try:
            min_epoch = float(np.nanmin(_e_arr))
            max_epoch = float(np.nanmax(_e_arr))
        except Exception:
            min_epoch = max_epoch = None
    else:
        min_epoch = max_epoch = None
    # PM arrays
    pmra_arr = source_df['pmra'].fillna(0.0).values if 'pmra' in source_df.columns else np.zeros(len(source_df))
    pmdec_arr = source_df['pmdec'].fillna(0.0).values if 'pmdec' in source_df.columns else np.zeros(len(source_df))
    pmra_err_arr = source_df['pmra_error'].fillna(0.0).values if 'pmra_error' in source_df.columns else np.zeros(len(source_df))
    pmdec_err_arr = source_df['pmdec_error'].fillna(0.0).values if 'pmdec_error' in source_df.columns else np.zeros(len(source_df))
    ref_epoch_arr = source_df['ref_epoch'].values if 'ref_epoch' in source_df.columns else np.full(len(source_df), np.nan)
    # Precompute PM speed (deg/yr) for each source
    pm_speed_deg_yr = np.hypot(pmra_arr, pmdec_arr) / 1000.0 / 3600.0

    matches = []
    for i in range(len(source_df)):
        # Offsets in degrees: RA scaled by cos(dec) for sky projection
        dec0_deg = float(source_df.iloc[i]['dec_deg'])
        dec0 = np.deg2rad(dec0_deg)
        ra0  = float(source_df.iloc[i]['ra_deg'])
        cosd = np.cos(dec0)
        dra  = (other_ra - ra0) * cosd
        ddec = (other_dec - dec0_deg)

        # Compute a conservative search radius (deg) for candidates
        pm_shift_max = 0.0
        if has_pm and has_epoch and (min_epoch is not None) and (max_epoch is not None):
            try:
                ref_ep = float(ref_epoch_arr[i])
            except Exception:
                ref_ep = np.nan
            if not (ref_ep is None or np.isnan(ref_ep)):
                dt_max = max(abs(min_epoch - ref_ep), abs(max_epoch - ref_ep))
                pm_shift_max = pm_speed_deg_yr[i] * float(dt_max)
        R_bound = float(r_src_bound[i]) + r_oth_max + pm_shift_max
        if R_bound <= 0:
            R_bound = float(r_src_bound[i]) + r_oth_max
        R2 = R_bound * R_bound
        # Vectorized prefilter by Euclidean distance in projected plane
        mask = (dra * dra + ddec * ddec) <= R2
        cand_idx = np.nonzero(mask)[0]

        ids = []
        # Evaluate only candidates with full Mahalanobis test
        for j in cand_idx:
            # Per-row epoch for PM shift
            epoch_j = epoch_list[j]
            dx_pm = 0.0
            dy_pm = 0.0
            if (epoch_j is not None) and has_pm and has_epoch:
                try:
                    ref_ep = float(ref_epoch_arr[i])
                except Exception:
                    ref_ep = None
                if ref_ep is not None and not np.isnan(ref_ep):
                    dt = float(epoch_j - ref_ep)
                    dx_pm = (pmra_arr[i] / 1000.0 / 3600.0) * dt
                    dy_pm = (pmdec_arr[i] / 1000.0 / 3600.0) * dt

            # Combined covariance matrix with PM error inflation if available
            C = src_cov[i] + oth_cov[j]
            if (epoch_j is not None) and has_pm_err and has_epoch:
                try:
                    ref_ep = float(ref_epoch_arr[i])
                except Exception:
                    ref_ep = None
                if ref_ep is not None and not np.isnan(ref_ep):
                    dt = float(epoch_j - ref_ep)
                    sig_ra_deg = abs(pmra_err_arr[i]) / 1000.0 / 3600.0 * abs(dt)
                    sig_de_deg = abs(pmdec_err_arr[i]) / 1000.0 / 3600.0 * abs(dt)
                    C = C + np.diag([sig_ra_deg**2, sig_de_deg**2])

            # Invert covariance, with fallback regularization if singular
            try:
                invC = la.inv(C)
            except la.LinAlgError:
                jitter = np.eye(2) * 1e-12
                try:
                    invC = la.inv(C + jitter)
                except la.LinAlgError:
                    continue
            vec = np.array([dra[j] - dx_pm, ddec[j] - dy_pm])
            D2 = vec.dot(invC).dot(vec)
            if D2 <= 1.0:
                if other_id_col in other_df.columns:
                    ids.append(other_df.iloc[j][other_id_col])
                else:
                    ids.append(other_df.index[j])
        matches.append(ids)

    # --- Precompute reverse assignment (other -> master) only if plotting PDF ---
    if pdf:
        assign_of_other = [None] * n_oth
        D2_best_arr = np.full(n_oth, np.inf, dtype=float)
        sep_best_arr = np.full(n_oth, np.nan, dtype=float)
        dec_src_rad_all = np.deg2rad(source_df['dec_deg'].values)
        for j in range(n_oth):
            dra_vec = (other_df['ra_deg'].iloc[j] - source_df['ra_deg'].values) * np.cos(dec_src_rad_all)
            ddec_vec = (other_df['dec_deg'].iloc[j] - source_df['dec_deg'].values)
            D2_vals = np.empty(n_src, dtype=float)
            for i in range(n_src):
                Csum = src_cov[i] + oth_cov[j]
                try:
                    invC = la.inv(Csum)
                except la.LinAlgError:
                    invC = la.inv(Csum + np.eye(2) * 1e-12)
                v = np.array([dra_vec[i], ddec_vec[i]])
                D2_vals[i] = float(v.dot(invC).dot(v))
            i_best = int(np.argmin(D2_vals)) if len(D2_vals) else None
            if len(D2_vals):
                D2_best_arr[j] = float(D2_vals[i_best])
                sep_best_arr[j] = float(np.hypot(dra_vec[i_best], ddec_vec[i_best]) * 3600.0)
                if D2_best_arr[j] <= 1.0:
                    assign_of_other[j] = i_best

    if pdf:
        # Helper: robust membership check for IDs (handles int/str mismatches)
        def _id_in_list(_ids, _val):
            try:
                if _val in _ids:
                    return True
            except Exception:
                pass
            try:
                sset = {str(x) for x in _ids}
                return str(_val) in sset
            except Exception:
                return False
        # Grid settings: 2 columns × 2 rows per PDF page
        ncols, nrows = 2, 2
        per_page = ncols * nrows
        total = len(other_df)
        # Precompute semi-axes for plotting according to plot_mode
        if plot_mode == 'match':
            _min_arcsec = float(min_radius.to(u.arcsec).value)
            oth_maj_plot = np.maximum(other_df['errMaj'].values * scale_factor, _min_arcsec)
            oth_min_plot = np.maximum(other_df['errMin'].values * scale_factor, _min_arcsec)
            src_maj_plot = np.maximum(source_df['errMaj'].values * scale_factor, _min_arcsec)
            src_min_plot = np.maximum(source_df['errMin'].values * scale_factor, _min_arcsec)
        else:
            oth_maj_plot = other_df['errMaj'].values
            oth_min_plot = other_df['errMin'].values
            src_maj_plot = source_df['errMaj'].values
            src_min_plot = source_df['errMin'].values
        # Build short-id maps per catalog for plotting (running numbers starting at 1)
        index_maps: dict[str, dict] = {}
        if all_catalogs is not None:
            for _k, _meta in all_catalogs.items():
                _df = _meta.get('data') if isinstance(_meta, dict) else None
                if _df is None:
                    continue
                _idx = list(_df.index)
                _map = {}
                for _i, _idval in enumerate(_idx, start=1):
                    _map[_idval] = _i
                    _map[str(_idval)] = _i
                index_maps[_k] = _map

        def _short_of(_key: str, _val):
            """Return running number for id `_val` within catalog `_key`, or None."""
            _m = index_maps.get(_key)
            if not _m:
                return None
            if _val in _m:
                return _m[_val]
            try:
                _v2 = int(_val)
                if _v2 in _m:
                    return _m[_v2]
            except Exception:
                pass
            return _m.get(str(_val))
        for start in range(0, total, per_page):
            end = min(start + per_page, total)
            indices = list(range(start, end))
            fig, axes = plt.subplots(nrows, ncols, figsize=(11, 8.5), constrained_layout=True)
            axes_flat = np.atleast_1d(axes).ravel()
            for ax in axes_flat[len(indices):]:
                fig.delaxes(ax)
            # Plot each source in this page
            hidden_texts = []
            for ax_idx, j in enumerate(indices):
                ax = axes_flat[ax_idx]
                # Reference position
                ra0 = other_df.iloc[j]['ra_deg']
                dec0 = other_df.iloc[j]['dec_deg']
                # Compute arcsec offsets for all new-catalog objects
                dx_oth = (other_df['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                dy_oth = (other_df['dec_deg'] - dec0) * 3600.0
                # Compute arcsec offsets for all master-catalog objects
                dx_src = (source_df['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                dy_src = (source_df['dec_deg'] - dec0) * 3600.0
                # Use precomputed nearest-master diagnostics
                D2_best = D2_best_arr[j]
                sep_best_arcsec = sep_best_arr[j]
                matched_flag = bool(D2_best <= 1.0)
                # Determine plotting margin in arcsec
                margin = max(other_df.iloc[j]['errMaj'], other_df.iloc[j]['errMin']) * 2
                # Enforce at least ±1 arcsec half‐width (5x5" area)
                margin = max(margin, 2.5)
                # Plot all new-catalog ellipses in light blue
                for k in range(len(other_df)):
                    if abs(dx_oth.iloc[k]) <= margin and abs(dy_oth.iloc[k]) <= margin:
                        e = Ellipse((dx_oth.iloc[k], dy_oth.iloc[k]),
                                    width=2*oth_min_plot[k],
                                    height=2*oth_maj_plot[k],
                                    angle=(90.0 - float(other_df.iloc[k]['errPA'])),
                                    edgecolor='lightblue', linestyle='-', fill=False)
                        ax.add_patch(e)
                # Plot all master-catalog ellipses in light red
                for m in range(len(source_df)):
                    if abs(dx_src.iloc[m]) <= margin and abs(dy_src.iloc[m]) <= margin:
                        e = Ellipse((dx_src.iloc[m], dy_src.iloc[m]),
                                    width=2*src_min_plot[m],
                                    height=2*src_maj_plot[m],
                                    angle=(90.0 - float(source_df.iloc[m]['errPA'])),
                                    edgecolor='lightcoral', linestyle='-', fill=False)
                        ax.add_patch(e)
                # Overlay and highlight the specific new source (blue solid)
                ell_new = Ellipse((0, 0),
                                  width=2*oth_min_plot[j],
                                  height=2*oth_maj_plot[j],
                                  angle=(90.0 - float(other_df.iloc[j]['errPA'])),
                                  edgecolor='blue', linestyle='-',
                                  fill=False, label=f"new: {other_id_col.replace('_id','')}")
                ax.add_patch(ell_new)
                # Overlay ALL cross‑identified master components (Gaia, 2MASS, WISE, Chandra, XMM)
                # Prefer the explicit reverse assignment based on Mahalanobis; fall back to matches scan
                master_idx = assign_of_other[j]
                if master_idx is None:
                    master_idx = next((i_idx for i_idx, ids in enumerate(matches)
                                       if _id_in_list(ids, other_df.iloc[j][other_id_col])), None)
                if master_idx is not None and all_catalogs is not None:
                    master_row = source_df.iloc[master_idx]
                    id_map = [
                        ('gaia_id', 'gaia', '-'),
                        ('2MASS', '2MASS', '--'),
                        ('wise_id', 'wise', ':'),
                        ('chandra_id', 'chandra', '-.'),
                        ('xmm_id', 'xmm', (0, (1, 1)))
                    ]
                    shown_labels = set()
                    for col, key, ls in id_map:
                        if col in source_df.columns and col in master_row:
                            val = master_row[col]
                            # skip missing/sentinel values
                            if pd.isna(val) or (isinstance(val, (int, np.integer, float)) and (val == -1 or (isinstance(val, float) and np.isnan(val)))):
                                continue
                            cat = all_catalogs.get(key)
                            if not cat:
                                continue
                            data = cat.get('data')
                            if data is None:
                                continue
                            # ensure index contains the id in its native type or as string
                            lookup_key = val
                            if lookup_key not in data.index:
                                try:
                                    lookup_key = int(val)
                                except Exception:
                                    lookup_key = str(val)
                            if lookup_key not in data.index:
                                continue
                            crow = data.loc[lookup_key]
                            dxm = (crow['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                            dym = (crow['dec_deg'] - dec0) * 3600.0
                            if plot_mode == 'match' and all_catalogs is not None and key in all_catalogs:
                                _f = float(all_catalogs[key].get('factor', 1.0))
                                _mr = float(all_catalogs[key].get('min_radius', 0.0))
                                _maj_m = max(float(crow['errMaj']) * _f, _mr)
                                _min_m = max(float(crow['errMin']) * _f, _mr)
                            else:
                                _maj_m = float(crow['errMaj'])
                                _min_m = float(crow['errMin'])
                            ell_m = Ellipse((dxm, dym),
                                            width=2*_min_m,
                                            height=2*_maj_m,
                                            angle=(90.0 - float(crow['errPA'])),
                                            edgecolor='red', linestyle='-',
                                            fill=False,
                                            label=(f'master:{key}' if f'master:{key}' not in shown_labels else None))
                            ax.add_patch(ell_m)
                            shown_labels.add(f'master:{key}')
                # Configure axes
                ax.set_xlim(-margin, margin)
                ax.set_ylim(-margin, margin)
                ax.set_aspect('equal', 'box')
                ax.invert_xaxis()
                ax.set_xlabel('ΔRA [arcsec]', labelpad=6)
                ax.set_ylabel('ΔDec [arcsec]', labelpad=6)
                ax.legend(loc='lower right', fontsize='small')
                # Build master_ids text using short running numbers; do not print the long catalog id inside the panel
                new_id = other_df.iloc[j][other_id_col]
                master_idx = assign_of_other[j]
                if master_idx is None:
                    master_idx = next((i_idx for i_idx, ids in enumerate(matches)
                                       if _id_in_list(ids, new_id)), None)
                ids_parts = []
                if master_idx is not None:
                    mrow = source_df.iloc[master_idx]
                    # use catalog keys matching those in `all_catalogs`
                    for col, key, label in [
                        ('gaia_id', 'gaia', 'gaia'),
                        ('2MASS',   '2MASS','2MASS'),
                        ('wise_id', 'wise', 'wise'),
                        ('chandra_id','chandra','chandra'),
                        ('xmm_id',  'xmm',  'xmm')
                    ]:
                        if col in source_df.columns and col in mrow:
                            val = mrow[col]
                            if pd.isna(val) or (isinstance(val, (int, np.integer, float)) and (val == -1 or (isinstance(val, float) and np.isnan(val)))):
                                continue
                            short = _short_of(key, val)
                            if short is not None:
                                ids_parts.append(f"{label}:{short}")
                ids_text = ' '.join(ids_parts) if ids_parts else 'none'
                # Only print master ids; omit the catalog id (now in the title)
                ax.text(0.05, 0.95,
                        f"master_ids: {ids_text}\nmin D^2: {D2_best:.2f}, sep: {sep_best_arcsec:.2f}\"",
                        transform=ax.transAxes, va='top', ha='left')

                # Title: show catalog name and running number within that catalog
                _cat_key = other_id_col.replace('_id', '')
                _short_new = _short_of(_cat_key, new_id)
                if _short_new is not None:
                    ax.set_title(f"{_cat_key} #{_short_new}", pad=10)
                else:
                    ax.set_title(f"{_cat_key}", pad=10)
            pdf.savefig(fig)
            plt.close(fig)
        pdf.close()

    return matches


def fetch_2mass_j_image(center: 'SkyCoord', width_deg: float, height_deg: float):
    """Fetch a 2MASS J-band cutout via SkyView (with retries),
    falling back to CDS hips2fits if needed. Images are cached on disk.

    Set environment variable TMASS_FETCH to control behaviour:
      - 'hips'     : use HiPS directly (default)
      - 'auto'     : try SkyView, then fall back to HiPS on failure
      - 'skyview'  : use SkyView only (no fallback)
      - 'off'      : disable image download entirely (no background)

    Returns (data, WCS) or (None, None) on failure.
    """
    # Cache path under images directory
    cache_name = os.path.join(
        IMAGES_DIR,
        f"2MASSJ_RA{center.ra.deg:.5f}_DEC{center.dec.deg:.5f}_W{width_deg:.4f}_H{height_deg:.4f}.fits"
    )
    refresh = globals().get("REFRESH", False)
    mode = os.environ.get("TMASS_FETCH", "hips").lower().strip()

    if (not refresh) and os.path.exists(cache_name):
        print(f"[2MASS J] cache hit → {cache_name}")
        try:
            with fits.open(cache_name, memmap=False) as hdul:
                data = hdul[0].data
                wcs = WCS(hdul[0].header)
            return data, wcs
        except Exception as e:
            print(f"[2MASS J] cache read failed ({e}), will re-download …")

    def _fetch_via_hips() -> tuple:
        """Final fallback: use CDS HiPS → FITS service."""
        try:
            print("[2MASS J] attempting hips2fits …")
            fov = max(width_deg, height_deg)
            # Aim for ~1"/pix, clamp size to [128, 2048]
            size = int(np.clip(round((fov * 3600.0) / 1.0), 128, 2048))
            params = {
                'hips': 'CDS/P/2MASS/J',
                'ra': center.ra.deg,
                'dec': center.dec.deg,
                'fov': fov,
                'width': size,
                'height': size,
                'projection': 'TAN',
                'format': 'fits',
            }
            endpoints = [
                'https://alasky.u-strasbg.fr/hips-image-services/hips2fits',
                'https://skyview.gsfc.nasa.gov/current/cgi/hips2fits.pl',
            ]
            sess = requests.Session()
            for base in endpoints:
                try:
                    print(f"[2MASS J] hips2fits GET {base} …")
                    resp = sess.get(base, params=params, timeout=30)
                    if resp.ok and resp.content:
                        with fits.open(io.BytesIO(resp.content), memmap=False) as hdul:
                            hdul.writeto(cache_name, overwrite=True)
                        print(f"[2MASS J] saved (hips2fits) → {cache_name}")
                        with fits.open(cache_name, memmap=False) as hdul:
                            data = hdul[0].data
                            wcs = WCS(hdul[0].header)
                        return data, wcs
                    else:
                        print(f"[2MASS J] hips2fits status {resp.status_code}")
                except Exception as e2:
                    print(f"[2MASS J] hips2fits error at {base}: {e2}")
        except Exception as e3:
            print(f"[2MASS J] hips2fits fallback error: {e3}")
        return None, None

    # If user forces HiPS, go straight there
    if mode == 'off':
        return None, None
    if mode == 'hips':
        return _fetch_via_hips()

    # Try multiple attempts with exponential backoff via SkyView
    surveys = ['2MASS J', '2MASS-J']
    max_tries = 4
    delays = [1, 3, 7, 15]
    last_error = None
    early_fallback = False

    # If user forces SkyView, we will not call HiPS later
    force_skyview = (mode == 'skyview')

    for attempt in range(max_tries):
        for survey in surveys:
            try:
                print(f"[2MASS J] attempt {attempt+1}/{max_tries} via SkyView.get_images(survey='{survey}') …")
                imgs = SkyView.get_images(position=center, survey=[survey],
                                          width=width_deg * u.deg, height=height_deg * u.deg)
                if imgs:
                    hdu = imgs[0][0]
                    hdu.writeto(cache_name, overwrite=True)
                    print(f"[2MASS J] saved → {cache_name}")
                    data = hdu.data
                    wcs = WCS(hdu.header)
                    return data, wcs
                else:
                    print(f"[2MASS J] SkyView.get_images returned no image for survey '{survey}'")
            except Exception as e:
                last_error = e
                msg = str(e)
                print(f"[2MASS J] get_images error (survey '{survey}'): {e}")
                # Early fallback if we hit a hard socket error
                if ('Connection aborted' in msg) or isinstance(e, (OSError, ConnectionError)):
                    early_fallback = True
                    break
        if early_fallback:
            break
        # Fallback within this attempt: try URL list then direct download via astropy
        try:
            print(f"[2MASS J] attempt {attempt+1}/{max_tries} via get_image_list + direct download …")
            urls = SkyView.get_image_list(position=center, survey=['2MASS J'],
                                          width=width_deg * u.deg, height=height_deg * u.deg)
            if not urls:
                urls = SkyView.get_image_list(position=center, survey=['2MASS-J'],
                                              width=width_deg * u.deg, height=height_deg * u.deg)
            if urls:
                url = urls[0]
                print(f"[2MASS J] downloading URL → {url}")
                tmp_path = download_file(url, cache=False, timeout=120)
                with fits.open(tmp_path, memmap=False) as hdul:
                    hdul.writeto(cache_name, overwrite=True)
                print(f"[2MASS J] saved → {cache_name}")
                with fits.open(cache_name, memmap=False) as hdul:
                    data = hdul[0].data
                    wcs = WCS(hdul[0].header)
                return data, wcs
            else:
                print("[2MASS J] get_image_list returned no URLs")
        except Exception as e:
            last_error = e
            print(f"[2MASS J] get_image_list/download error: {e}")
            if 'Connection aborted' in str(e):
                early_fallback = True
        # Backoff before next attempt unless we will fall back now
        if early_fallback:
            break
        if attempt < max_tries - 1:
            delay = delays[min(attempt, len(delays)-1)]
            print(f"[2MASS J] retrying in {delay}s …")
            try:
                time.sleep(delay)
            except KeyboardInterrupt:
                print("[2MASS J] interrupted during backoff")
                break

    if force_skyview:
        print("[2MASS J] ERROR: SkyView mode failed and HiPS fallback disabled (TMASS_FETCH=skyview)")
        return None, None

    # Fall back to HiPS if allowed
    data, wcs = _fetch_via_hips()
    if data is not None:
        return data, wcs

    # All attempts failed
    if last_error is not None:
        print(f"[2MASS J] ERROR: exhausted retries: {last_error}")
    else:
        print("[2MASS J] ERROR: exhausted retries with no specific exception")
    return None, None

# === Post-merge plotting so panels reflect final combined_df ===
def plot_after_merge(combined_df: pd.DataFrame,
                     other_df: pd.DataFrame,
                     id_col: str,
                     all_catalogs: dict,
                     pdf_path: str,
                     plot_mode: str = 'match',
                     ncols: int = 2,
                     nrows: int = 2,
                     invert_cmap: bool = False,
                     draw_images: bool = True,
                     aladin_dir: str | None = None,
                     samp_enabled: bool = False,
                     samp_addr: str | None = None) -> None:
    key = id_col.replace('_id', '')
    params = all_catalogs.get(key, {}) if isinstance(all_catalogs, dict) else {}
    factor = float(params.get('factor', 1.0))
    min_r = float(params.get('min_radius', 0.0))
    # Color palette per catalog (consistent across plot and table)
    CAT_COLORS = {
        'gaia':    '#1f77b4',  # blue
        '2MASS':   '#2ca02c',  # green
        'wise':    '#ff7f0e',  # orange
        'chandra': '#d62728',  # red
        'xmm':     '#9467bd',  # purple
    }
    cur_color = CAT_COLORS.get(key, '#17becf')  # color for the current (new) catalog

    # Avoid pyplot figure warning for many pages
    try:
        import matplotlib as _mpl
        _mpl.rcParams['figure.max_open_warning'] = 0
    except Exception:
        pass

    # Optional: establish SAMP connection once per figure to talk to Aladin Desktop
    samp_client = None
    samp_aladin_cids: list[str] = []
    if samp_enabled and ('SAMPIntegratedClient' in globals()) and (SAMPIntegratedClient is not None):
        try:
            if samp_addr:
                samp_client = SAMPIntegratedClient(name="NGC2264-Matcher", addr=samp_addr)
            else:
                samp_client = SAMPIntegratedClient(name="NGC2264-Matcher")
            # Bind to loopback to avoid OS routing issues
            try:
                samp_client.connect()
            except Exception:
                # Retry on IPv6 loopback
                samp_client = SAMPIntegratedClient(name="NGC2264-Matcher", addr="::1")
                samp_client.connect()
            # Discover Aladin-like clients
            for cid in samp_client.get_registered_clients():
                try:
                    meta = samp_client.get_metadata(cid)
                    name = ' '.join([str(meta.get(k, '')) for k in ('samp.name','application.name','author')]).lower()
                    if ('aladin' in name) or ('cds' in name):
                        samp_aladin_cids.append(cid)
                except Exception:
                    continue
            # Keep logs terse: only report once per call
            if not samp_aladin_cids:
                print('[SAMP] Warning: hub connected but no Aladin client found')
        except Exception as e:
            # Degrade quietly when SAMP is unavailable
            print(f"[SAMP] Connection failed: {e}")
            samp_client = None

    # Collect entries for an HTML index with per-panel links
    index_rows: list[dict] = []
    # Collect per-page HTML outputs
    page_files: list[str] = []

    # Helpers for epochs and PM-inflated ellipse drawing
    def _epoch_from_row(row) -> Optional[float]:
        for c in ['ref_epoch','epoch','EPOCH','Epoch']:
            if c in other_df.columns:
                try:
                    v = float(row[c])
                    if not np.isnan(v):
                        return v
                except Exception:
                    pass
        mjds = []
        for c in ['MJD','mjd','w1mjdmean','w2mjdmean','w3mjdmean','w4mjdmean']:
            if c in other_df.columns:
                try:
                    v = float(row[c])
                    if not np.isnan(v) and v > 0:
                        mjds.append(v)
                except Exception:
                    pass
        if mjds:
            mjd = float(np.nanmean(mjds))
            return 2000.0 + (mjd - 51544.5) / 365.25
        for c in ['JD','jd']:
            if c in other_df.columns:
                try:
                    jd = float(row[c])
                    if not np.isnan(jd) and jd > 0:
                        return 2000.0 + (jd - 2451545.0) / 365.25
                except Exception:
                    pass
        for c in ['date','dateobs','obs_date']:
            if c in other_df.columns:
                try:
                    from datetime import datetime
                    s = str(row[c])
                    dt = datetime.fromisoformat(s)
                    year_start = datetime(dt.year, 1, 1)
                    doy = (dt - year_start).days + 0.5
                    return dt.year + doy/365.25
                except Exception:
                    pass
        meta = all_catalogs.get(key, {}) if isinstance(all_catalogs, dict) else {}
        try:
            return float(meta.get('epoch')) if meta and (meta.get('epoch') is not None) else None
        except Exception:
            return None

    def _cov_from_axes_arcsec(maj_as, min_as, pa_deg):
        phi = np.deg2rad(pa_deg)
        sphi, cphi = np.sin(phi), np.cos(phi)
        vmaj = np.array([sphi, cphi])
        vmin = np.array([-cphi, sphi])
        return (maj_as*maj_as) * np.outer(vmaj, vmaj) + (min_as*min_as) * np.outer(vmin, vmin)

    def _ellipse_from_cov_arcsec(C):
        vals, vecs = np.linalg.eigh(C)
        order = np.argsort(vals)[::-1]
        vals = vals[order]
        vecs = vecs[:, order]
        maj = float(np.sqrt(max(vals[0], 0.0)))
        min_ = float(np.sqrt(max(vals[1], 0.0)))
        v = vecs[:, 0]
        # Position angle from North through East (PA_NE)
        pa_ne = float((np.degrees(np.arctan2(v[0], v[1])) % 360.0))
        # For plotting with ΔRA on x and inverted x-axis, use angle = PA - 90
        return maj, min_, (pa_ne - 90.0)

    # Build short-id maps from all catalogs
    index_maps: dict[str, dict] = {}
    for _k, _meta in all_catalogs.items():
        _df = _meta.get('data') if isinstance(_meta, dict) else None
        if _df is None:
            continue
        _idx = list(_df.index)
        _map = {}
        for _i, _idval in enumerate(_idx, start=1):
            _map[_idval] = _i
            _map[str(_idval)] = _i
        index_maps[_k] = _map

    def _short_of(_key: str, _val):
        _m = index_maps.get(_key)
        if not _m:
            return None
        if _val in _m:
            return _m[_val]
        try:
            _v2 = int(_val)
            if _v2 in _m:
                return _m[_v2]
        except Exception:
            pass
        return _m.get(str(_val))

    def _is_missing(val) -> bool:
        return (pd.isna(val)
                or (isinstance(val, (int, np.integer, float))
                    and ((isinstance(val, float) and np.isnan(val)) or val == -1)))

    pdf = PdfPages(pdf_path) if pdf_path else None
    per_page = ncols * nrows
    total = len(other_df)
    import math as _math
    total_pages = max(1, _math.ceil(total / max(1, per_page)))

    # Precompute new-catalog ellipse axes for plotting
    if plot_mode == 'match':
        oth_maj_plot = np.maximum(other_df['errMaj'].values * factor, min_r)
        oth_min_plot = np.maximum(other_df['errMin'].values * factor, min_r)
    else:
        oth_maj_plot = other_df['errMaj'].values
        oth_min_plot = other_df['errMin'].values

    # Allow restricting generation to a single page via env: ALADIN_PAGE_ONLY
    try:
        _page_only = int(str(os.environ.get('ALADIN_PAGE_ONLY', '0')))
    except Exception:
        _page_only = 0
    for start in range(0, total, per_page):
        page_num = (start // per_page) + 1
        if _page_only > 0 and page_num != _page_only:
            continue
        end = min(start + per_page, total)
        indices = list(range(start, end))
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 8.5), constrained_layout=True)
        axes_flat = np.atleast_1d(axes).ravel()
        for ax in axes_flat[len(indices):]:
            fig.delaxes(ax)
        # Collect clickable areas for this page (axes positions and links)
        page_areas: list[dict] = []
        for ax_idx, j in enumerate(indices):
            ax = axes_flat[ax_idx]
            new_id = other_df.iloc[j][id_col]
            ra0 = other_df.iloc[j]['ra_deg']
            dec0 = other_df.iloc[j]['dec_deg']
            epoch_row = _epoch_from_row(other_df.iloc[j])

            # Offsets for background
            dx_oth = (other_df['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
            dy_oth = (other_df['dec_deg'] - dec0) * 3600.0
            dx_src = (combined_df['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
            dy_src = (combined_df['dec_deg'] - dec0) * 3600.0
            # Global nearest Gaia master by angular separation (arcsec), excluding the synthetic self row
            if 'gaia_id' in combined_df.columns:
                valid_master_mask = (~combined_df['gaia_id'].isna()) & (combined_df['gaia_id'] != -1)
            else:
                valid_master_mask = np.ones(len(combined_df), dtype=bool)
            try:
                not_self_mask = combined_df[id_col].astype(str) != str(new_id)
            except Exception:
                not_self_mask = np.ones(len(combined_df), dtype=bool)
            nearest_mask = valid_master_mask & not_self_mask
            if np.any(nearest_mask):
                dist_all = np.hypot(dx_src[nearest_mask].values, dy_src[nearest_mask].values)
                nearest_sep_arcsec = float(np.min(dist_all))
            else:
                nearest_sep_arcsec = None
            # PM-corrected closest separation (computed later when we know Δt); used for title display
            closest_sep_pm_arcsec = None

            # Determine plotting margin in arcsec
            margin = max(oth_maj_plot[j], oth_min_plot[j]) * 2.0
            margin = max(margin, 5.0)  # at least ±2.5" half-width

            # If this is a 2MASS panel and allowed, draw a 2MASS J-band background aligned via WCS
            if draw_images and key == '2MASS':
                center_icrs = SkyCoord(ra=ra0 * u.deg, dec=dec0 * u.deg, frame='icrs')
                fov_deg = (2.0 * margin) / 3600.0  # full width/height in degrees
                data_img, wcs_img = fetch_2mass_j_image(center_icrs, fov_deg, fov_deg)
                if data_img is not None and wcs_img is not None:
                    ny, nx = data_img.shape[0], data_img.shape[1]
                    # Use mid-lines to anchor left/right and bottom/top explicitly
                    x_left, x_right = 0.0, float(nx - 1)
                    y_bot, y_top = 0.0, float(ny - 1)
                    x_mid, y_mid = (nx - 1) / 2.0, (ny - 1) / 2.0

                    # World at left/right midpoints
                    lr_pix = np.array([[x_left, y_mid], [x_right, y_mid]], dtype=float)
                    lr_world = wcs_img.wcs_pix2world(lr_pix, 0)
                    ra_left,  ra_right  = lr_world[0, 0], lr_world[1, 0]
                    # World at bottom/top midpoints
                    bt_pix = np.array([[x_mid, y_bot], [x_mid, y_top]], dtype=float)
                    bt_world = wcs_img.wcs_pix2world(bt_pix, 0)
                    dec_bot, dec_top = bt_world[0, 1], bt_world[1, 1]

                    # Convert to arcsec offsets relative to panel center (east-positive, north-positive)
                    dx_left  = (ra_left  - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                    dx_right = (ra_right - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                    dy_bot_o = (dec_bot  - dec0) * 3600.0
                    dy_top_o = (dec_top  - dec0) * 3600.0

                    # Preserve left→right and bottom→top orientation; do not sort
                    extent = [float(dx_left), float(dx_right), float(dy_bot_o), float(dy_top_o)]

                    # Contrast-friendly normalization: compress dynamic range so faint sources remain visible
                    # Controls (optional): TMASS_PERCENT (default 99.8), TMASS_ASINH_A (default 0.2)
                    try:
                        percent = float(os.environ.get('TMASS_PERCENT', '99.8'))
                    except Exception:
                        percent = 99.8
                    try:
                        a_param = float(os.environ.get('TMASS_ASINH_A', '0.2'))
                    except Exception:
                        a_param = 0.2

                    # Compute robust min/max using upper percentile; lower bound implied by (100 - percent)/2
                    interval = PercentileInterval(percent)
                    vmin, vmax = interval.get_limits(data_img)
                    norm = ImageNormalize(vmin=vmin, vmax=vmax,
                                           stretch=AsinhStretch(a=a_param), clip=True)

                    # Optional inversion of grayscale via CLI flag or env var TMASS_INVERT
                    invert_env = str(os.environ.get('TMASS_INVERT', '0')).lower() in ('1', 'true', 'yes', 'on')
                    cmap_name = 'gray_r' if (invert_cmap or invert_env) else 'gray'

                    ax.imshow(
                        data_img,
                        cmap=cmap_name,
                        origin='lower',
                        extent=extent,
                        zorder=0,
                        interpolation='nearest',
                        norm=norm
                    )

            # Plot new-catalog background (light blue)
            for k in range(len(other_df)):
                if abs(dx_oth.iloc[k]) <= margin and abs(dy_oth.iloc[k]) <= margin:
                    e = Ellipse((dx_oth.iloc[k], dy_oth.iloc[k]),
                                width=2*oth_min_plot[k],
                                height=2*oth_maj_plot[k],
                                angle=(90.0 - float(other_df.iloc[k]['errPA'])),
                                edgecolor='lightblue', linestyle='-', fill=False)
                    ax.add_patch(e)

            # Collect candidate master rows within the plotted area; per-catalog plotting happens below
            cand_indices = [int(m) for m in range(len(combined_df))
                            if abs(dx_src.iloc[m]) <= margin and abs(dy_src.iloc[m]) <= margin]
            # Suppress synthetic self-row (current catalog id with gaia_id == -1) from the table
            if (id_col in combined_df.columns) and ('gaia_id' in combined_df.columns):
                _new_id_str = str(new_id)
                filtered = []
                for m in cand_indices:
                    try:
                        same_id = (str(combined_df.iloc[m][id_col]) == _new_id_str)
                    except Exception:
                        same_id = False
                    is_synth = (pd.isna(combined_df.iloc[m].get('gaia_id', np.nan))
                                or combined_df.iloc[m].get('gaia_id', -1) == -1)
                    if same_id and is_synth:
                        continue
                    filtered.append(m)
                cand_indices = filtered

            # Overlay the specific new source
            e0 = Ellipse((0, 0),
                         width=2*oth_min_plot[j], height=2*oth_maj_plot[j],
                         angle=(90.0 - float(other_df.iloc[j]['errPA'])),
                         edgecolor=cur_color, linestyle='-', fill=False,
                         label=f"new: {key}")
            ax.add_patch(e0)

            # Find ALL combined_df rows with this id (no inference) and overlay all their components
            master_indices: list[int] = []
            if id_col in combined_df.columns:
                try:
                    mask = combined_df[id_col] == new_id
                except Exception:
                    mask = combined_df[id_col].astype(str) == str(new_id)
                master_indices = [int(ix) for ix in list(combined_df.index[mask])]

            ids_parts = []
            seen_ids = set()
            drawn_components = set()  # avoid drawing the same component twice in a panel
            if master_indices:
                for master_idx in master_indices:
                    mrow = combined_df.iloc[master_idx]
                    # Overlay all components listed in the combined row
                    for col, kcat, ls in [
                        ('gaia_id', 'gaia', '-'),
                        ('2MASS',   '2MASS','--'),
                        ('wise_id', 'wise', ':'),
                        ('chandra_id','chandra','-.'),
                        ('xmm_id',  'xmm',  (0,(1,1)))
                    ]:
                        # Skip overlay of the same catalog as the new one to avoid duplication
                        if kcat == key:
                            continue
                        if col in combined_df.columns and not _is_missing(mrow.get(col, None)):
                            val = mrow[col]
                            # collect short id for annotation (deduplicated)
                            short_val = _short_of(kcat, val)
                            if short_val is not None and (kcat, short_val) not in seen_ids:
                                ids_parts.append(f"{kcat}:{short_val}")
                                seen_ids.add((kcat, short_val))
                            # draw component ellipse (deduplicated by catalog+id)
                            cat = all_catalogs.get(kcat)
                            if not cat:
                                continue
                            data = cat.get('data')
                            if data is None:
                                continue
                            lookup_key = val
                            if lookup_key not in data.index:
                                try:
                                    lookup_key = int(val)
                                except Exception:
                                    lookup_key = str(val)
                            if lookup_key not in data.index:
                                continue
                            if (kcat, lookup_key) in drawn_components:
                                continue
                            crow = data.loc[lookup_key]
                            dxm = (crow['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                            dym = (crow['dec_deg'] - dec0) * 3600.0
                            if plot_mode == 'match' and kcat in all_catalogs:
                                _f = float(all_catalogs[kcat].get('factor', 1.0))
                                _mr = float(all_catalogs[kcat].get('min_radius', 0.0))
                                _maj_m = max(float(crow['errMaj']) * _f, _mr)
                                _min_m = max(float(crow['errMin']) * _f, _mr)
                            else:
                                _maj_m = float(crow['errMaj'])
                                _min_m = float(crow['errMin'])
                            lw = 1.8 if kcat == 'gaia' else 1.2
                            ell_m = Ellipse((dxm, dym), width=2*_min_m, height=2*_maj_m,
                                            angle=(90.0 - float(crow['errPA'])), edgecolor=CAT_COLORS.get(kcat, 'red'), linestyle='-',
                                            linewidth=lw, fill=False, label=None)
                            ax.add_patch(ell_m)
                            # Emphasize matched GAIA position with a colored plus + PM arrow/inflated ellipse
                            if kcat == 'gaia':
                                ax.plot(dxm, dym, marker='+', color=CAT_COLORS['gaia'],
                                        markersize=7, markeredgewidth=1.2, linestyle='None')
                                # If epoch is known and PM available, draw arrow and PM-inflated ellipse
                                try:
                                    if (epoch_row is not None
                                        and ('pmra' in crow.index) and ('pmdec' in crow.index)
                                        and ('pmra_error' in crow.index) and ('pmdec_error' in crow.index)
                                        and ('ref_epoch' in crow.index)):
                                        dt = float(epoch_row - float(crow['ref_epoch']))
                                        pm_dx = float(crow['pmra']) * 0.001 * dt   # arcsec along RA*
                                        pm_dy = float(crow['pmdec']) * 0.001 * dt  # arcsec Dec
                                        if np.hypot(pm_dx, pm_dy) > 0.05:
                                            ax.arrow(dxm, dym, pm_dx, pm_dy, width=0.0, head_width=0.03,
                                                     head_length=0.05, length_includes_head=True,
                                                     color=CAT_COLORS['gaia'], alpha=0.9)
                                        # PM-uncertainty inflation
                                        sig_ra = abs(float(crow['pmra_error'])) * 0.001 * abs(dt)
                                        sig_de = abs(float(crow['pmdec_error'])) * 0.001 * abs(dt)
                                        Cpos = _cov_from_axes_arcsec(_maj_m, _min_m, float(crow['errPA']))
                                        Cinf = Cpos + np.diag([sig_ra**2, sig_de**2])
                                        maj_i, min_i, pa_i = _ellipse_from_cov_arcsec(Cinf)
                                        ell_pm = Ellipse((dxm, dym), width=2*maj_i, height=2*min_i,
                                                         angle=pa_i, edgecolor=CAT_COLORS['gaia'], linestyle='--',
                                                         linewidth=1.2, fill=False, alpha=0.9)
                                        ax.add_patch(ell_pm)
                                except Exception:
                                    pass
                            drawn_components.add((kcat, lookup_key))

            # Background: draw components for UNMATCHED master rows in lighter/transparent colors
            # and draw a thin gray cross at every master position
            if 'gaia_id' in combined_df.columns:
                matched_set = {m for m in master_indices
                            if (not pd.isna(combined_df.iloc[m].get('gaia_id', np.nan))
                                and combined_df.iloc[m].get('gaia_id', -1) != -1)}
            else:
                matched_set = set(master_indices)                        
            for m in cand_indices:
                if m not in matched_set:
                    # gray cross at master position (unmatched only)
                    ax.plot(dx_src.iloc[m], dy_src.iloc[m], marker='+', color='gray', markersize=5,
                            markeredgewidth=0.6, linestyle='None')
                else:
                    # matched rows are emphasized elsewhere (e.g., GAIA plus & thicker ellipse)
                    pass
                if m in matched_set:
                    continue  # matched components already drawn above with solid colors
                mrow = combined_df.iloc[m]
                for colname, kcat, ls in [
                    ('gaia_id','gaia','-'),
                    ('2MASS','2MASS','--'),
                    ('wise_id','wise',':'),
                    ('chandra_id','chandra','-.'),
                    ('xmm_id','xmm',(0,(1,1)))
                ]:
                    if kcat == key:
                        continue  # skip the current catalog
                    if colname not in combined_df.columns or _is_missing(mrow.get(colname, None)):
                        continue
                    cat = all_catalogs.get(kcat)
                    if not cat or cat.get('data') is None:
                        continue
                    val = mrow[colname]
                    lookup_key = val if val in cat['data'].index else (int(val) if str(val).isdigit() else str(val))
                    if lookup_key not in cat['data'].index:
                        continue
                    crow = cat['data'].loc[lookup_key]
                    dxm = (crow['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                    dym = (crow['dec_deg'] - dec0) * 3600.0
                    if plot_mode == 'match' and kcat in all_catalogs:
                        _f = float(all_catalogs[kcat].get('factor', 1.0))
                        _mr = float(all_catalogs[kcat].get('min_radius', 0.0))
                        _maj_m = max(float(crow['errMaj']) * _f, _mr)
                        _min_m = max(float(crow['errMin']) * _f, _mr)
                    else:
                        _maj_m = float(crow['errMaj'])
                        _min_m = float(crow['errMin'])
                    ell_bg = Ellipse((dxm, dym), width=2*_maj_m, height=2*_min_m,
                                     angle=(float(crow['errPA']) - 90.0), edgecolor=CAT_COLORS.get(kcat, 'red'), linestyle='-',
                                     linewidth=0.8, fill=False, alpha=0.28)
                    ax.add_patch(ell_bg)

            # --- Compute compact diagnostics (best separation and Mahalanobis D²) for title ---
            def _cov_from_axes(maj_as, min_as, pa_deg):
                phi = np.deg2rad(pa_deg)
                sphi, cphi = np.sin(phi), np.cos(phi)
                vmaj = np.array([sphi, cphi])
                vmin = np.array([-cphi, sphi])
                return (maj_as*maj_as) * np.outer(vmaj, vmaj) + (min_as*min_as) * np.outer(vmin, vmin)

            best_sep = None
            best_d2 = None
            best_dt = None
            nearest_d2 = None
            nearest_dt = None
            if master_indices:
                # Keep only real master rows (with a Gaia id), excluding synthetic rows created from the new catalog
                if 'gaia_id' in combined_df.columns:
                    master_valid = [m for m in master_indices
                                    if (not pd.isna(combined_df.iloc[m].get('gaia_id', np.nan))
                                        and combined_df.iloc[m].get('gaia_id', -1) != -1)]
                else:
                    master_valid = master_indices
                if master_valid:
                    Cnew = _cov_from_axes(oth_maj_plot[j], oth_min_plot[j], float(other_df.iloc[j]['errPA']))
                    for m in master_valid:
                        # PM-corrected residual vector Gaia→catalog row (arcsec)
                        pm_dx = 0.0; pm_dy = 0.0
                        if (epoch_row is not None) and ('pmra' in combined_df.columns) and ('pmdec' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                            try:
                                ref_ep = float(combined_df.iloc[m]['ref_epoch'])
                                dt_tmp = float(epoch_row - ref_ep)
                                pm_dx = float(combined_df.iloc[m]['pmra']) * 0.001 * dt_tmp
                                pm_dy = float(combined_df.iloc[m]['pmdec']) * 0.001 * dt_tmp
                            except Exception:
                                pm_dx = 0.0; pm_dy = 0.0
                        vec = np.array([float(dx_src.iloc[m]) + pm_dx, float(dy_src.iloc[m]) + pm_dy])
                        mrow = combined_df.iloc[m]
                        Cmas = _cov_from_axes(float(mrow['errMaj']), float(mrow['errMin']), float(mrow['errPA']))
                        # PM-uncertainty inflation for this candidate baseline
                        if (epoch_row is not None) and ('pmra_error' in combined_df.columns) and ('pmdec_error' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                            try:
                                sig_ra = abs(float(mrow['pmra_error'])) * 0.001 * abs(dt_tmp)
                                sig_de = abs(float(mrow['pmdec_error'])) * 0.001 * abs(dt_tmp)
                                Cmas = Cmas + np.diag([sig_ra**2, sig_de**2])
                            except Exception:
                                pass
                        Csum = Cnew + Cmas
                        try:
                            invC = la.inv(Csum)
                        except la.LinAlgError:
                            invC = la.inv(Csum + np.eye(2) * 1e-6)
                        d2 = float(vec.dot(invC).dot(vec))
                        sep = float(np.hypot(*vec))
                        if best_d2 is None or d2 < best_d2:
                            best_d2, best_sep = d2, sep
                            if epoch_row is not None and ('ref_epoch' in combined_df.columns):
                                try:
                                    best_dt = float(epoch_row - float(mrow['ref_epoch']))
                                except Exception:
                                    best_dt = None
                else:
                    # There is a synthetic self-row but no valid Gaia master: compute D² to the nearest valid master
                    if 'gaia_id' in combined_df.columns:
                        valid_master_mask = (~combined_df['gaia_id'].isna()) & (combined_df['gaia_id'] != -1)
                    else:
                        valid_master_mask = np.ones(len(combined_df), dtype=bool)
                    try:
                        not_self_mask = combined_df[id_col].astype(str) != str(new_id)
                    except Exception:
                        not_self_mask = np.ones(len(combined_df), dtype=bool)
                    nearest_mask = valid_master_mask & not_self_mask
                    if np.any(nearest_mask):
                        idxs = np.where(nearest_mask)[0]
                        dist_all = np.hypot(dx_src.iloc[idxs].values, dy_src.iloc[idxs].values)
                        k = int(np.argmin(dist_all))
                        mnear = int(idxs[k])
                        Cnew = _cov_from_axes(oth_maj_plot[j], oth_min_plot[j], float(other_df.iloc[j]['errPA']))
                        mrow = combined_df.iloc[mnear]
                        Cmas = _cov_from_axes(float(mrow['errMaj']), float(mrow['errMin']), float(mrow['errPA']))
                        if (epoch_row is not None) and ('pmra_error' in combined_df.columns) and ('pmdec_error' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                            try:
                                ref_ep = float(mrow['ref_epoch']); dt_tmp = float(epoch_row - ref_ep)
                                sig_ra = abs(float(mrow['pmra_error'])) * 0.001 * abs(dt_tmp)
                                sig_de = abs(float(mrow['pmdec_error'])) * 0.001 * abs(dt_tmp)
                                Cmas = Cmas + np.diag([sig_ra**2, sig_de**2])
                            except Exception:
                                pass
                        Csum = Cnew + Cmas
                        try:
                            invC = la.inv(Csum)
                        except la.LinAlgError:
                            invC = la.inv(Csum + np.eye(2) * 1e-6)
                        # PM-corrected residual vector to this nearest master
                        pm_dx = 0.0; pm_dy = 0.0
                        if (epoch_row is not None) and ('pmra' in combined_df.columns) and ('pmdec' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                            try:
                                ref_ep = float(mrow['ref_epoch']); dt_tmp = float(epoch_row - ref_ep)
                                pm_dx = float(mrow['pmra']) * 0.001 * dt_tmp
                                pm_dy = float(mrow['pmdec']) * 0.001 * dt_tmp
                            except Exception:
                                pm_dx = 0.0; pm_dy = 0.0
                        vec = np.array([float(dx_src.iloc[mnear]) + pm_dx, float(dy_src.iloc[mnear]) + pm_dy])
                        nearest_d2 = float(vec.dot(invC).dot(vec))
                        closest_sep_pm_arcsec = float(np.hypot(*vec))
                        try:
                            nearest_dt = float(epoch_row - float(mrow['ref_epoch'])) if (epoch_row is not None) else None
                        except Exception:
                            nearest_dt = None
            else:
                # No master rows with the same id at all: compute D² to the nearest valid master by angular separation
                if 'gaia_id' in combined_df.columns:
                    valid_master_mask = (~combined_df['gaia_id'].isna()) & (combined_df['gaia_id'] != -1)
                else:
                    valid_master_mask = np.ones(len(combined_df), dtype=bool)
                try:
                    not_self_mask = combined_df[id_col].astype(str) != str(new_id)
                except Exception:
                    not_self_mask = np.ones(len(combined_df), dtype=bool)
                nearest_mask = valid_master_mask & not_self_mask
                if np.any(nearest_mask):
                    idxs = np.where(nearest_mask)[0]
                    dist_all = np.hypot(dx_src.iloc[idxs].values, dy_src.iloc[idxs].values)
                    k = int(np.argmin(dist_all))
                    mnear = int(idxs[k])
                    Cnew = _cov_from_axes(oth_maj_plot[j], oth_min_plot[j], float(other_df.iloc[j]['errPA']))
                    mrow = combined_df.iloc[mnear]
                    Cmas = _cov_from_axes(float(mrow['errMaj']), float(mrow['errMin']), float(mrow['errPA']))
                    if (epoch_row is not None) and ('pmra_error' in combined_df.columns) and ('pmdec_error' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                        try:
                            ref_ep = float(mrow['ref_epoch']); dt_tmp = float(epoch_row - ref_ep)
                            sig_ra = abs(float(mrow['pmra_error'])) * 0.001 * abs(dt_tmp)
                            sig_de = abs(float(mrow['pmdec_error'])) * 0.001 * abs(dt_tmp)
                            Cmas = Cmas + np.diag([sig_ra**2, sig_de**2])
                        except Exception:
                            pass
                    Csum = Cnew + Cmas
                    try:
                        invC = la.inv(Csum)
                    except la.LinAlgError:
                        invC = la.inv(Csum + np.eye(2) * 1e-6)
                    # PM-corrected residual vector to this nearest master
                    pm_dx = 0.0; pm_dy = 0.0
                    if (epoch_row is not None) and ('pmra' in combined_df.columns) and ('pmdec' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                        try:
                            ref_ep = float(mrow['ref_epoch']); dt_tmp = float(epoch_row - ref_ep)
                            pm_dx = float(mrow['pmra']) * 0.001 * dt_tmp
                            pm_dy = float(mrow['pmdec']) * 0.001 * dt_tmp
                        except Exception:
                            pm_dx = 0.0; pm_dy = 0.0
                    vec = np.array([float(dx_src.iloc[mnear]) + pm_dx, float(dy_src.iloc[mnear]) + pm_dy])
                    nearest_d2 = float(vec.dot(invC).dot(vec))
                    closest_sep_pm_arcsec = float(np.hypot(*vec))
                    try:
                        nearest_dt = float(epoch_row - float(mrow['ref_epoch'])) if (epoch_row is not None) else None
                    except Exception:
                        nearest_dt = None

            # === Dynamic association table across available catalogs (exclude current catalog) ===
            # Determine which catalogs to include based on the candidate rows; exclude current catalog `key`
            cat_col_map = [
                ('gaia',    'gaia_id'),
                ('2MASS',   '2MASS'),
                ('wise',    'wise_id'),
                ('chandra', 'chandra_id'),
                ('xmm',     'xmm_id'),
            ]
            col_for = dict(cat_col_map)
            include = set()
            for m in cand_indices:
                mrow = combined_df.iloc[m]
                for kcat, colname in cat_col_map:
                    if kcat == key:
                        continue
                    if colname in combined_df.columns and not _is_missing(mrow.get(colname, None)):
                        include.add(kcat)
            display_cats = [k for k, _ in cat_col_map if (k in include and k != key)]

            # Build rows sorted by distance from the current catalog source.
            # Prefer Gaia (master) positions when available to avoid shifts from
            # earlier updates of combined_df coords.
            if 'gaia_id' in combined_df.columns:
                matched_set = {m for m in master_indices
                            if (not pd.isna(combined_df.iloc[m].get('gaia_id', np.nan))
                                and combined_df.iloc[m].get('gaia_id', -1) != -1)}
            else:
                matched_set = set(master_indices)
            try:
                _gaia_df = all_catalogs.get('gaia', {}).get('data') if isinstance(all_catalogs, dict) else None
            except Exception:
                _gaia_df = None
            cand_with_dist = []
            for m in cand_indices:
                try:
                    if (_gaia_df is not None) and ('gaia_id' in combined_df.columns):
                        gid = combined_df.iloc[m].get('gaia_id', None)
                        if (gid is not None) and (gid in _gaia_df.index):
                            grow = _gaia_df.loc[gid]
                            dxm = (float(grow['ra_deg']) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                            dym = (float(grow['dec_deg']) - dec0) * 3600.0
                            dist_as = float(np.hypot(dxm, dym))
                        else:
                            dist_as = float(np.hypot(dx_src.iloc[m], dy_src.iloc[m]))
                    else:
                        dist_as = float(np.hypot(dx_src.iloc[m], dy_src.iloc[m]))
                except Exception:
                    dist_as = float(np.hypot(dx_src.iloc[m], dy_src.iloc[m]))
                cand_with_dist.append((int(m), dist_as))
            cand_with_dist.sort(key=lambda t: t[1])

            def _fmt(val):
                return str(val) if val is not None else '—'

            # Build table columns and headers; draw G magnitude near the right panel edge (not in table)
            col_keys = list(display_cats)
            add_gmag = ('gaia' in display_cats) and ('Gmag' in combined_df.columns)
            headers = [k.upper() for k in col_keys]
            raw_rows = []  # list of (cells_list, is_matched)
            g_values = []  # per-row G values (strings) to render at right edge
            for m, _dist in cand_with_dist:
                mrow = combined_df.iloc[m]
                cells = []
                for kcat in col_keys:
                    colname = col_for[kcat]
                    if colname in combined_df.columns and not _is_missing(mrow.get(colname, None)):
                        cells.append(_fmt(_short_of(kcat, mrow[colname])))
                    else:
                        cells.append(_fmt(None))
                raw_rows.append((cells, (m in matched_set)))
                # Right-edge G value
                if add_gmag and ('gaia_id' in combined_df.columns) and not _is_missing(mrow.get('gaia_id', None)) \
                        and pd.notna(mrow.get('Gmag', np.nan)):
                    try:
                        gval = float(mrow['Gmag'])
                        g_values.append(f"{gval:.2f}")
                    except Exception:
                        g_values.append(_fmt(None))
                else:
                    g_values.append(_fmt(None))

            if headers:
                # Compute column widths from headers and rows (in characters)
                table_data = [headers] + [cells for (cells, _) in raw_rows]
                col_w = [max(len(str(r[i])) for r in table_data) for i in range(len(headers))]

                # --- Scale table size with panel grid ---
                base_font = 9
                base_line_h = 0.040  # slightly tighter than before
                # Scale roughly with panel area (relative to 2x2)
                area_scale = (4.0 / float(max(1, ncols * nrows))) ** 0.5
                font_size = int(round(np.clip(base_font * area_scale, 6, 12)))
                line_h = base_line_h * (font_size / base_font)

                # Adjust top anchor a bit when many panels are stacked
                y0 = 0.95 - 0.02 * max(0, max(nrows, ncols) - 2)
                y0 = max(0.85, min(0.95, y0))

                # Checkmark space in first column if any matched rows
                has_match_row = any(is_m for (_, is_m) in raw_rows)
                if has_match_row and len(col_w) > 0:
                    col_w[0] += 2  # reserve space for '✓ '

                # Compute how many rows fit (leave a small bottom margin)
                bottom_margin = 0.10
                # Minus only the header line now (no separator)
                max_rows = max(1, int((y0 - bottom_margin) / line_h) - 1)
                rows_limited = raw_rows[:max_rows]

                # Character width in points for monospace font (approx)
                char_w_pts = font_size * 0.6

                # Compute x-offsets in points for each column (include two spaces between columns)
                x_offsets = []
                acc = 0.0
                for i in range(len(headers)):
                    x_offsets.append(acc)
                    acc += (col_w[i] + 2) * char_w_pts

                # Draw header per column with catalog colors (no separator line)
                for i, head in enumerate(headers):
                    trans = offset_copy(ax.transAxes, fig=fig, x=x_offsets[i], y=0, units='points')
                    key_i = col_keys[i]
                    key_color = CAT_COLORS.get(key_i, 'black')
                    ax.text(0.05, y0, head, transform=trans,
                            va='top', ha='left', family='monospace', fontsize=font_size,
                            fontweight='bold', color=key_color)
                # Right-edge G header
                if add_gmag:
                    ax.text(0.95, y0, 'G', transform=ax.transAxes,
                            va='top', ha='right', family='monospace', fontsize=font_size,
                            fontweight='bold', color=CAT_COLORS.get('gaia','black'))

                # Draw rows per column, color-coded and bold if matched; add '✓ ' on matched rows
                for r_idx, (cells, is_matched) in enumerate(rows_limited):
                    y = y0 - (r_idx + 1) * line_h
                    for i, cell in enumerate(cells):
                        text_cell = str(cell)
                        if i == 0:
                            prefix = '✓ ' if is_matched else '  '
                            text_cell = prefix + text_cell
                        trans = offset_copy(ax.transAxes, fig=fig, x=x_offsets[i], y=0, units='points')
                        key_i = col_keys[i]
                        key_color = CAT_COLORS.get(key_i, 'black')
                        ax.text(0.05, y, text_cell.ljust(col_w[i]), transform=trans,
                                va='top', ha='left', family='monospace', fontsize=font_size,
                                fontweight=('bold' if is_matched else 'normal'),
                                color=key_color)
                    # Right-edge G value
                    if add_gmag and r_idx < len(g_values):
                        ax.text(0.95, y, str(g_values[r_idx]), transform=ax.transAxes,
                                va='top', ha='right', family='monospace', fontsize=font_size,
                                fontweight=('bold' if is_matched else 'normal'),
                                color=CAT_COLORS.get('gaia','black'))

                # If there are more rows than fit, add a small ellipsis line
                if len(rows_limited) < len(raw_rows):
                    y = y0 - (len(rows_limited) + 1) * line_h
                    trans = offset_copy(ax.transAxes, fig=fig, x=x_offsets[0], y=0, units='points')
                    ax.text(0.05, y, f"… (+{len(raw_rows) - len(rows_limited)} more)", transform=trans,
                            va='top', ha='left', family='monospace', fontsize=font_size, color='gray')


            # Title with short id for the current catalog + compact diagnostics
            short_new = _short_of(key, new_id)
            title = f"{key} #{short_new}" if short_new is not None else f"{key}"
            # Mark potential blends when multiple valid Gaia masters share this 2MASS id
            if key == '2MASS' and 'gaia_id' in combined_df.columns:
                try:
                    same_id_mask = (combined_df[id_col].astype(str) == str(new_id))
                except Exception:
                    same_id_mask = np.zeros(len(combined_df), dtype=bool)
                if np.any(same_id_mask):
                    valid_gaia = same_id_mask & (~combined_df['gaia_id'].isna()) & (combined_df['gaia_id'] != -1)
                    if int(np.sum(valid_gaia)) >= 2:
                        title += " — blend"
            if (best_sep is not None) and (best_d2 is not None):
                # Use Unicode double-prime for arcsec and superscript 2 for D²
                title += f" — sep={best_sep:.2f}″, D²={best_d2:.2f}"
                if best_dt is not None:
                    title += f", Δt={best_dt:.1f}y"
            elif (nearest_sep_arcsec is not None) or (closest_sep_pm_arcsec is not None):
                # No matches: report separation to the closest master source
                if nearest_d2 is not None:
                    sep_to_show = (closest_sep_pm_arcsec
                                   if closest_sep_pm_arcsec is not None else nearest_sep_arcsec)
                    title += f" — closest={sep_to_show:.2f}″, D²={nearest_d2:.2f}"
                    if nearest_dt is not None:
                        title += f", Δt={nearest_dt:.1f}y"
                else:
                    sep_to_show = (closest_sep_pm_arcsec
                                   if closest_sep_pm_arcsec is not None else nearest_sep_arcsec)
                    title += f" — closest={sep_to_show:.2f}″"
            # Decide if this panel is problematic and mark it visibly
            try:
                if 'gaia_id' in combined_df.columns:
                    num_valid_master = int(np.sum((combined_df.loc[master_indices, 'gaia_id'].astype(float) != -1)
                                                  & (~combined_df.loc[master_indices, 'gaia_id'].isna()))) if master_indices else 0
                else:
                    num_valid_master = len(master_indices)
            except Exception:
                num_valid_master = len(master_indices) if master_indices else 0

            try:
                near_thr = float(max(3.0, float(oth_min_plot[j])))
            except Exception:
                near_thr = 3.0
            try:
                near_sep = float(closest_sep_pm_arcsec) if (closest_sep_pm_arcsec is not None) else (float(nearest_sep_arcsec) if (nearest_sep_arcsec is not None) else None)
            except Exception:
                near_sep = None
            close_count = 0
            try:
                if master_indices:
                    for _m in master_indices:
                        if 'gaia_id' in combined_df.columns and (pd.isna(combined_df.iloc[_m].get('gaia_id', np.nan))
                                                                or combined_df.iloc[_m].get('gaia_id', -1) == -1):
                            continue
                        d_arcsec = float(np.hypot(dx_src.iloc[_m], dy_src.iloc[_m]))
                        if d_arcsec <= near_thr:
                            close_count += 1
            except Exception:
                pass
            is_problem = False
            if (num_valid_master == 0 and (near_sep is not None) and (near_sep <= near_thr)):
                is_problem = True
            elif (num_valid_master >= 1 and close_count >= 2):
                is_problem = True
            if is_problem:
                try:
                    for sp in ax.spines.values():
                        sp.set_linewidth(2.0)
                        sp.set_edgecolor('#d62728')
                except Exception:
                    pass
                try:
                    ax.text(0.02, 0.02, '⚠', transform=ax.transAxes, ha='left', va='bottom',
                            fontsize=18, color='#d62728', fontweight='bold', zorder=20,
                            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='#d62728', lw=1.2, alpha=0.9))
                except Exception:
                    pass
            ax.set_title(title, pad=10)

            ax.set_xlim(-margin, margin)
            ax.set_ylim(-margin, margin)
            ax.set_aspect('equal', 'box')
            ax.invert_xaxis()
            ax.set_xlabel('ΔRA [arcsec]', labelpad=6)
            ax.set_ylabel('ΔDec [arcsec]', labelpad=6)
            legend_handles = []
            # Current catalog (new source)
            legend_handles.append(Line2D([], [], color=cur_color, label=f'{key} (new)'))
            # Faint style indicator for unmatched
            legend_handles.append(Line2D([], [], color='black', alpha=0.28, label='unmatched (faint)'))
            # Master positions
            legend_handles.append(Line2D([], [], color='gray', marker='+', linestyle='None', label='master pos'))
            ax.legend(handles=legend_handles, loc='lower right', fontsize='small', framealpha=0.9)

            # Add a clickable link to open the same view in Aladin (Lite link),
            # and optionally write a Desktop Aladin script + region file with a file:// link
            try:
                from urllib.parse import quote_plus
                # Center and FOV in degrees consistent with this panel
                fov_deg = (2.0 * margin) / 3600.0
                target = f"{ra0:.7f} {dec0:.7f}"
                # Choose a survey close to the current catalog for visual consistency
                survey_map = {
                    '2MASS': 'P/2MASS/J',
                    'wise': 'P/AllWISE/color',
                    'gaia': 'P/DSS2/color',
                    'chandra': 'P/DSS2/color',
                    'xmm': 'P/DSS2/color',
                }
                survey = survey_map.get(key, 'P/DSS2/color')
                url = (
                    'https://aladin.u-strasbg.fr/AladinLite/?'
                    f'target={quote_plus(target)}&fov={fov_deg:.6f}&survey={quote_plus(survey)}'
                )
                text_lite = ax.text(
                    0.98, 0.98, 'Aladin Lite', transform=ax.transAxes,
                    ha='right', va='top', fontsize=8, color='#1f77b4',
                    url=url, zorder=10,
                    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#1f77b4', lw=0.8, alpha=0.95)
                )
                try:
                    hidden_texts.append(text_lite)
                except Exception:
                    pass

                # If a directory is provided, create an Aladin Desktop script (.ajs) and a DS9 region file
                if aladin_dir:
                    os.makedirs(aladin_dir, exist_ok=True)
                    # Sanitize an id for filenames
                    short_new = _short_of(key, new_id)
                    tag = f"{key}_{short_new}" if short_new is not None else f"{key}"
                    base = f"aladin_{tag}_c{start+ax_idx+1}"
                    reg_path = os.path.abspath(os.path.join(aladin_dir, base + '.reg'))
                    ajs_path = os.path.abspath(os.path.join(aladin_dir, base + '.ajs'))

                    # Build a DS9 region content for the panel overlays (fk5)
                    reg = []
                    reg.append('# Region file format: DS9')
                    reg.append('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')
                    reg.append('fk5')
                    # New catalog ellipse (match-mode scaled)
                    reg.append(
                        'ellipse(%.7f,%.7f,%.3f\",%.3f\",%.1f) # color=%s width=2' % (
                            ra0, dec0, oth_maj_plot[j], oth_min_plot[j], float(other_df.iloc[j]['errPA']),
                            'green'
                        )
                    )
                    # Overlay Gaia and others in the panel area
                    def add_point(r, d, color='blue', shape='cross', size=8):
                        reg.append('point(%.7f,%.7f) # point=%s color=%s size=%d' % (r, d, shape, color, size))
                    def add_line(r1, d1, r2, d2, color='blue'):
                        reg.append('line(%.7f,%.7f,%.7f,%.7f) # color=%s' % (r1, d1, r2, d2, color))

                    # Add Gaia points and PM arrows where available
                    if 'gaia' in all_catalogs:
                        gaia_df = all_catalogs['gaia']['data']
                        # reuse cand_indices to restrict nearby content
                        for m in cand_indices:
                            if 'gaia_id' not in combined_df.columns:
                                break
                            if pd.isna(combined_df.iloc[m].get('gaia_id', np.nan)) or combined_df.iloc[m].get('gaia_id', -1) == -1:
                                continue
                            gid = combined_df.iloc[m]['gaia_id']
                            if gid in gaia_df.index:
                                grow = gaia_df.loc[gid]
                                add_point(float(grow['ra_deg']), float(grow['dec_deg']), color='blue', shape='cross', size=8)
                                # PM arrow if epoch known
                                if epoch_row is not None and all(k in grow for k in ['pmra','pmdec','ref_epoch']):
                                    try:
                                        dt = float(epoch_row - float(grow['ref_epoch']))
                                        dra_deg = (float(grow['pmra']) / 3600e3) * (1.0 / np.cos(np.deg2rad(float(grow['dec_deg'])))) * dt
                                        ddec_deg = (float(grow['pmdec']) / 3600e3) * dt
                                        add_line(float(grow['ra_deg']), float(grow['dec_deg']),
                                                 float(grow['ra_deg']) + dra_deg, float(grow['dec_deg']) + ddec_deg,
                                                 color='blue')
                                    except Exception:
                                        pass

                    # Write region file
                    with open(reg_path, 'w') as f:
                        f.write('\n'.join(reg) + '\n')

                    # Build a minimal Aladin script: center, FoV, and overlay planes
                    survey_aladin = {
                        '2MASS': 'CDS/P/2MASS/J',
                        'wise': 'CDS/P/AllWISE/Color',
                    }.get(key, 'CDS/P/DSS2/color')
                    fov_arcmin = fov_deg * 60.0
                    ajs = [
                        f'info "PDF script: {base}"',
                        # Keep current background; just move and adjust FoV
                        f'{ra0:.7f} {dec0:.7f}',
                        f'zoom {fov_arcmin:.3f} arcmin',
                    ]

                    # Build colorized overlays via Aladin draw commands in separate planes
                    # First, remove previous planes (if any) for idempotent refresh
                    # Remove any overlays created by previous clicks
                    ajs.append('rm "aladin_*_ellipse"')
                    ajs.append('rm "aladin_*_gaia"')
                    ajs.append('rm "aladin_*"')
                    ajs.append(f'rm "{base}_ellipse"')
                    ajs.append(f'rm "{base}_gaia"')
                    # Also try to remove a previously loaded region file plane (both full path and basename)
                    try:
                        _reg_base = os.path.basename(reg_path)
                    except Exception:
                        _reg_base = ''
                    if _reg_base:
                        ajs.append(f'rm "{_reg_base}"')
                    ajs.append(f'rm "{reg_path}"')
                    # Plane 1: ellipse for the new/other catalog source (green)
                    ajs.append(f'draw newtool("{base}_ellipse")')
                    ajs.append(
                        'draw green ellipse(%.7f %.7f %.3farcsec %.3farcsec %.1f)'
                        % (ra0, dec0, oth_maj_plot[j], oth_min_plot[j], float(other_df.iloc[j]['errPA']))
                    )

                    # Plane 2: Gaia candidates and proper motion arrows (blue)
                    if 'gaia' in all_catalogs:
                        ajs.append(f'draw newtool("{base}_gaia")')
                        gaia_df = all_catalogs['gaia']['data']
                        g_factor = float(all_catalogs['gaia'].get('factor', 1.0))
                        g_minr   = float(all_catalogs['gaia'].get('min_radius', 0.0))
                        for m in cand_indices:
                            if 'gaia_id' not in combined_df.columns:
                                break
                            gid = combined_df.iloc[m].get('gaia_id', -1)
                            if pd.isna(gid) or gid == -1:
                                continue
                            if gid in gaia_df.index:
                                grow = gaia_df.loc[gid]
                                try:
                                    r_g = float(grow['ra_deg']); d_g = float(grow['dec_deg'])
                                except Exception:
                                    continue
                                # Identification ellipse for Gaia, scaled as in the PDF panels
                                try:
                                    maj_g = max(float(grow['errMaj']) * g_factor, g_minr)
                                    min_g = max(float(grow['errMin']) * g_factor, g_minr)
                                    pa_g  = float(grow['errPA'])
                                except Exception:
                                    maj_g = float(grow.get('errMaj', 0.5))
                                    min_g = float(grow.get('errMin', 0.5))
                                    pa_g  = float(grow.get('errPA', 0.0))
                                ajs.append('draw blue ellipse(%.7f %.7f %.3farcsec %.3farcsec %.1f)'
                                           % (r_g, d_g, maj_g, min_g, pa_g))
                                # PM arrow if epoch known
                                if epoch_row is not None and all(k in grow for k in ['pmra','pmdec','ref_epoch']):
                                    try:
                                        dt = float(epoch_row - float(grow['ref_epoch']))
                                        dra_deg = (float(grow['pmra']) / 3600e3) * (1.0 / max(1e-8, np.cos(np.deg2rad(d_g)))) * dt
                                        ddec_deg = (float(grow['pmdec']) / 3600e3) * dt
                                        ajs.append('draw blue line(%.7f %.7f %.7f %.7f)'
                                                   % (r_g, d_g, r_g + dra_deg, d_g + ddec_deg))
                                    except Exception:
                                        pass
                    # Optional: load the DS9 region file (can re-introduce accumulation in some Aladin versions)
                    if str(os.environ.get('ALADIN_LOAD_REG', '0')).lower() in ('1','true','yes','on'):
                        ajs.append(f'load "{reg_path}"')
                        ajs.append('sync')
                    with open(ajs_path, 'w') as f:
                        f.write('\n'.join(ajs) + '\n')

                    # Try to create a small launcher for Desktop Aladin, to avoid PDF/file association issues
                    launcher_path = None
                    launcher_open_path = None
                    launcher_openfile_path = None
                    launcher_samp_path = None
                    launcher_samp_app_path = None
                    # 1) Try JAR-based launcher
                    try:
                        jar_env = os.environ.get('ALADIN_JAR', '').strip()
                        jar_candidates = [p for p in [jar_env] if p]
                        # Common macOS bundle locations (newer and older)
                        jar_candidates += [
                            '/Applications/Aladin.app/Contents/app/Aladin.jar',
                            '/Applications/Aladin.app/Contents/Java/Aladin.jar',
                        ]
                        jar_path = next((p for p in jar_candidates if os.path.exists(p)), '')
                        if jar_path:
                            if sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
                                ext = '.command' if sys.platform.startswith('darwin') else '.sh'
                                launcher_path = os.path.abspath(os.path.join(aladin_dir, base + ext))
                                script = [
                                    '#!/usr/bin/env bash',
                                    'set -euo pipefail',
                                    f'JAR="{jar_path}"',
                                    f'exec java ${ALADIN_JAVA_OPTS:-} -jar "$JAR" -script "{ajs_path}"'
                                ]
                                with open(launcher_path, 'w') as f:
                                    f.write('\n'.join(script) + '\n')
                                try:
                                    os.chmod(launcher_path, 0o755)
                                except Exception:
                                    pass
                                # Also create macOS launchers that target the running app via open -a
                                if sys.platform.startswith('darwin'):
                                    launcher_open_path = os.path.abspath(os.path.join(aladin_dir, base + '_open.command'))
                                    aladin_app = os.environ.get('ALADIN_APP', '/Applications/Aladin.app')
                                    script_open = [
                                        '#!/usr/bin/env bash',
                                        'set -euo pipefail',
                                        f'APP="{aladin_app}"',
                                        f'open -a "$APP" --args -script "{ajs_path}"',
                                        'sleep 0.2 || true',
                                        "osascript -e 'tell application \"Aladin\" to activate' || true",
                                    ]
                                    with open(launcher_open_path, 'w') as f:
                                        f.write('\n'.join(script_open) + '\n')
                                    try:
                                        os.chmod(launcher_open_path, 0o755)
                                    except Exception:
                                        pass
                                    # Alternative: pass the .ajs as a document (relies on .ajs association)
                                    launcher_openfile_path = os.path.abspath(os.path.join(aladin_dir, base + '_openfile.command'))
                                    script_openfile = [
                                        '#!/usr/bin/env bash',
                                        'set -euo pipefail',
                                        f'APP="{aladin_app}"',
                                        f'open -a "$APP" "{ajs_path}"',
                                        'sleep 0.2 || true',
                                        "osascript -e 'tell application \"Aladin\" to activate' || true",
                                    ]
                                    with open(launcher_openfile_path, 'w') as f:
                                        f.write('\n'.join(script_openfile) + '\n')
                                    try:
                                        os.chmod(launcher_openfile_path, 0o755)
                                    except Exception:
                                        pass
                            elif sys.platform.startswith('win'):
                                ext = '.bat'
                                launcher_path = os.path.abspath(os.path.join(aladin_dir, base + ext))
                                script = [
                                    '@echo off',
                                    'setlocal',
                                    f'set JAR={jar_path}',
                                    f'java %ALADIN_JAVA_OPTS% -jar "%JAR%" -script "{ajs_path}"'
                                ]
                                with open(launcher_path, 'w') as f:
                                    f.write('\r\n'.join(script) + '\r\n')
                    except Exception:
                        launcher_path = None

                    # 2) macOS fallback: open -a Aladin.app
                    if (launcher_path is None) and sys.platform.startswith('darwin'):
                        try:
                            ext = '.command'
                            launcher_path = os.path.abspath(os.path.join(aladin_dir, base + ext))
                            aladin_app = os.environ.get('ALADIN_APP', '/Applications/Aladin.app')
                            script = [
                                '#!/usr/bin/env bash',
                                'set -euo pipefail',
                                f'APP="{aladin_app}"',
                                f'open -a "$APP" --args -script "{ajs_path}"',
                                'sleep 0.2 || true',
                                "osascript -e 'tell application \"Aladin\" to activate' || true",
                            ]
                            with open(launcher_path, 'w') as f:
                                f.write('\n'.join(script) + '\n')
                            try:
                                os.chmod(launcher_path, 0o755)
                            except Exception:
                                pass
                        except Exception:
                            launcher_path = None

                    # 3) SAMP launcher to drive the running Aladin instance (preferred for in-place control)
                    try:
                        samp_py = os.path.abspath(os.path.join(os.getcwd(), 'aladin_samp_test.py'))
                        if os.path.exists(samp_py):
                            if sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
                                ext = '.command' if sys.platform.startswith('darwin') else '.sh'
                                launcher_samp_path = os.path.abspath(os.path.join(aladin_dir, base + '_samp' + ext))
                                script = [
                                    '#!/usr/bin/env bash',
                                    'set -euo pipefail',
                                    'PYTHON_BIN="${PYTHON:-python3}"',
                                    'ADDR="${SAMP_ADDR:-127.0.0.1}"',
                                    f'"$PYTHON_BIN" "{samp_py}" --mode script --send-ajs "{ajs_path}" --addr "$ADDR"',
                                ]
                                with open(launcher_samp_path, 'w') as f:
                                    f.write('\n'.join(script) + '\n')
                                try:
                                    os.chmod(launcher_samp_path, 0o755)
                                except Exception:
                                    pass
                                # On macOS, also try to generate a small AppleScript .app to avoid opening Terminal
                                if sys.platform.startswith('darwin') and os.path.exists('/usr/bin/osacompile'):
                                    try:
                                        launcher_samp_app_path = os.path.abspath(os.path.join(aladin_dir, base + '_samp.app'))
                                        applescript_src = os.path.abspath(os.path.join(aladin_dir, base + '_samp.applescript'))
                                        pybin = os.environ.get('PYTHON', sys.executable if sys.executable else 'python3')
                                        # Build a background shell command; escape quotes for AppleScript
                                        cmd = f'{pybin} "{samp_py}" --mode script --send-ajs "{ajs_path}" --addr 127.0.0.1 >/dev/null 2>&1 &'
                                        cmd_esc = cmd.replace('\\', '\\\\').replace('"', '\\"')
                                        content = 'on run\n    do shell script "' + cmd_esc + '"\nend run\n'
                                        with open(applescript_src, 'w') as f:
                                            f.write(content)
                                        # Compile AppleScript into an app bundle
                                        os.system(f'/usr/bin/osacompile -o "{launcher_samp_app_path}" "{applescript_src}" >/dev/null 2>&1')
                                        # Clean up source file if compile succeeded
                                        if os.path.exists(launcher_samp_app_path):
                                            try:
                                                os.remove(applescript_src)
                                            except Exception:
                                                pass
                                    except Exception:
                                        launcher_samp_app_path = None
                            elif sys.platform.startswith('win'):
                                ext = '.bat'
                                launcher_samp_path = os.path.abspath(os.path.join(aladin_dir, base + '_samp' + ext))
                                script = [
                                    '@echo off',
                                    'setlocal',
                                    'set PYTHON_BIN=%PYTHON%',
                                    'if "%PYTHON_BIN%"=="" set PYTHON_BIN=python',
                                    'set ADDR=%SAMP_ADDR%',
                                    'if "%ADDR%"=="" set ADDR=127.0.0.1',
                                    f'"%PYTHON_BIN%" "{samp_py}" --mode script --send-ajs "{ajs_path}" --addr "%ADDR%"'
                                ]
                                with open(launcher_samp_path, 'w') as f:
                                    f.write('\r\n'.join(script) + '\r\n')
                    except Exception:
                        launcher_samp_path = None

                    # Left-side button: keep only one SAMP link to control the running Aladin
                    if launcher_samp_app_path and os.path.exists(launcher_samp_app_path):
                        file_url = 'file://' + launcher_samp_app_path
                        primary_label = 'Aladin'
                    elif launcher_samp_path and os.path.exists(launcher_samp_path):
                        file_url = 'file://' + launcher_samp_path
                        primary_label = 'Aladin'
                    else:
                        file_url = 'file://' + os.path.abspath(aladin_dir)
                        primary_label = 'Open Scripts Folder'
                    text_aladin = ax.text(
                        0.02, 0.98, primary_label, transform=ax.transAxes,
                        ha='left', va='top', fontsize=8, color='#2ca02c',
                        url=file_url, zorder=10,
                        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#2ca02c', lw=0.8, alpha=0.95)
                    )
                    try:
                        hidden_texts.append(text_aladin)
                    except Exception:
                        pass

                    # Determine if this panel is potentially problematic
                    try:
                        # Count valid master components (with Gaia id)
                        if 'gaia_id' in combined_df.columns:
                            num_valid_master = int(np.sum((combined_df.loc[master_indices, 'gaia_id'].astype(float) != -1) & (~combined_df.loc[master_indices, 'gaia_id'].isna()))) if master_indices else 0
                        else:
                            num_valid_master = len(master_indices)
                    except Exception:
                        num_valid_master = len(master_indices) if master_indices else 0
                    # Use nearest separation (PM-corrected if available)
                    try:
                        near_sep = float(closest_sep_pm_arcsec) if (closest_sep_pm_arcsec is not None) else (float(nearest_sep_arcsec) if (nearest_sep_arcsec is not None) else None)
                    except Exception:
                        near_sep = None
                    # Threshold: larger of 3" and the new-source minor axis
                    try:
                        near_thr = float(max(3.0, float(oth_min_plot[j])))
                    except Exception:
                        near_thr = 3.0
                    # If there are multiple nearby masters within threshold, count them
                    close_count = 0
                    try:
                        if master_indices:
                            for _m in master_indices:
                                # Skip synthetic rows if present
                                if 'gaia_id' in combined_df.columns and (pd.isna(combined_df.iloc[_m].get('gaia_id', np.nan)) or combined_df.iloc[_m].get('gaia_id', -1) == -1):
                                    continue
                                d_arcsec = float(np.hypot(dx_src.iloc[_m], dy_src.iloc[_m]))
                                if d_arcsec <= near_thr:
                                    close_count += 1
                    except Exception:
                        pass
                    is_problem = False
                    # Problematic if: (no identifications and a nearby master) or (an identification exists and another close-by master is present)
                    if (num_valid_master == 0 and (near_sep is not None) and (near_sep <= near_thr)):
                        is_problem = True
                    elif (num_valid_master >= 1 and close_count >= 2):
                        is_problem = True

                    # Save an entry for the HTML index
                    try:
                        index_rows.append({
                            'base': base,
                            'aladin': file_url,
                            'lite': url,
                        })
                    except Exception:
                        pass

                    # Save area info for page-level HTML image map
                    try:
                        page_areas.append({
                            'base': base,
                            'ax_index': ax_idx,
                            'aladin': file_url,
                            'lite': url,
                            'problem': 1 if is_problem else 0,
                        })
                    except Exception:
                        pass

                    # Also prepare a simple VOTable overlay and send via SAMP (safe, no scripts)
                    if samp_client is not None and samp_aladin_cids:
                        try:
                            vot_path = os.path.abspath(os.path.join(aladin_dir, base + '.vot'))
                            rows = [("new", f"{ra0:.7f}", f"{dec0:.7f}")]
                            # Add visible Gaia candidates in the panel area
                            if 'gaia' in all_catalogs and 'gaia_id' in combined_df.columns:
                                for m in cand_indices:
                                    gid = combined_df.iloc[m].get('gaia_id', -1)
                                    if pd.isna(gid) or gid == -1:
                                        continue
                                    # skip duplicates of the same Gaia id
                                    rows.append((f"gaia_{int(gid)}",
                                                 f"{float(combined_df.iloc[m]['ra_deg']):.7f}",
                                                 f"{float(combined_df.iloc[m]['dec_deg']):.7f}"))
                            # Minimal VOTable
                            vot = [
                                "<?xml version=\"1.0\" encoding=\"utf-8\"?>",
                                "<VOTABLE version=\"1.3\" xmlns=\"http://www.ivoa.net/xml/VOTable/v1.3\">",
                                " <RESOURCE>",
                                "  <TABLE>",
                                "   <FIELD name=\"name\" datatype=\"char\" arraysize=\"*\"/>",
                                "   <FIELD name=\"RAJ2000\" ucd=\"pos.eq.ra;meta.main\" unit=\"deg\" datatype=\"double\"/>",
                                "   <FIELD name=\"DEJ2000\" ucd=\"pos.eq.dec;meta.main\" unit=\"deg\" datatype=\"double\"/>",
                                "   <DATA>",
                                "    <TABLEDATA>",
                            ]
                            # Limit duplicates by keeping order but unique names
                            seen = set()
                            for name, ra_s, de_s in rows:
                                if name in seen:
                                    continue
                                vot.append(f"     <TR><TD>{name}</TD><TD>{ra_s}</TD><TD>{de_s}</TD></TR>")
                                seen.add(name)
                            vot += [
                                "    </TABLEDATA>",
                                "   </DATA>",
                                "  </TABLE>",
                                " </RESOURCE>",
                                "</VOTABLE>",
                            ]
                            with open(vot_path, 'w') as f:
                                f.write('\n'.join(vot) + '\n')
                            # Center then load VOTable (no prompts in Aladin)
                            center_msg = {'samp.mtype': 'coord.pointAt.sky',
                                          'samp.params': {'ra': f'{ra0:.7f}', 'dec': f'{dec0:.7f}'}}
                            table_msg = {'samp.mtype': 'table.load.votable',
                                         'samp.params': {'url': 'file://' + vot_path}}
                            for cid in samp_aladin_cids:
                                try:
                                    samp_client.notify(cid, center_msg)
                                    samp_client.notify(cid, table_msg)
                                except Exception:
                                    continue
                        except Exception as e:
                            print(f"[SAMP] overlay failed: {e}")
            except Exception:
                pass
        # End panel loop for this page: save the figure
        # Save this page as PNG for HTML rendering
        try:
            scripts_dir = os.environ.get('ALADIN_SCRIPTS_DIR', 'aladin_scripts')
            html_dir = os.path.join(scripts_dir, 'html')
            os.makedirs(html_dir, exist_ok=True)
            try:
                nav_js_path = os.path.join(html_dir, 'nav_keys.js')
                nav_js = (
                    "// Keyboard navigation: ←/→ pages; 'a' to trigger Aladin.\n"
                    "(function(){\n"
                    "  var _mx, _my; document.addEventListener('mousemove', function(ev){ _mx=ev.clientX; _my=ev.clientY; }, {passive:true});\n"
                    "  function isEditable(el){ if(!el) return false; var tag=(el.tagName||'').toUpperCase(); return el.isContentEditable||tag==='INPUT'||tag==='TEXTAREA'||tag==='SELECT'; }\n"
                    "  function closest(el, sel){ while(el && el.nodeType===1){ if(el.matches && el.matches(sel)) return el; el = el.parentElement; } return null; }\n"
                    "  function findNav(which){\n"
                    "    var navs = document.querySelectorAll('.nav');\n"
                    "    for(var i=0;i<navs.length;i++){\n"
                    "      var as = navs[i].querySelectorAll('a[href]');\n"
                    "      for(var j=0;j<as.length;j++){\n"
                    "        var a = as[j]; var t = (a.textContent||'').trim();\n"
                    "        if(which==='prev' && t.indexOf('Prev Page') !== -1) return a;\n"
                    "        if(which==='next' && t.indexOf('Next Page') !== -1) return a;\n"
                    "      }\n"
                    "    }\n"
                    "    var m = (location.pathname||'').match(/([^\\/]+)_page(\\d+)\\.html$/);\n"
                    "    if(m){ var base = m[1]; var n = parseInt(m[2],10); var target = which==='prev' ? (n-1) : (n+1); if(target>=1){ var a = document.createElement('a'); a.setAttribute('href', base + '_page' + target + '.html'); return a; } }\n"
                    "    return null;\n"
                    "  }\n"
                    "  function triggerAladin(){\n"
                    "    var cand=null;\n"
                    "    if(typeof _mx==='number' && typeof _my==='number'){ var el=document.elementFromPoint(_mx,_my); cand = closest(el, 'a.btn.al[data-ajs]'); }\n"
                    "    if(!cand) cand = document.querySelector('a.btn.al[data-ajs]');\n"
                    "    if(cand){ try{ cand.click(); return true; }catch(e){} }\n"
                    "    try{ if(typeof ALADIN_ITEMS!=='undefined' && ALADIN_ITEMS.length){ var idx=(window._srcIdx===undefined?0:window._srcIdx); var ajs = ALADIN_ITEMS[idx] || ALADIN_ITEMS[0]; var u=(typeof SAMP_URL!=='undefined'?SAMP_URL:'http://127.0.0.1:8765'); var base=(u && u.endsWith && u.endsWith('/'))?u.slice(0,-1):u; var url = base + '/run_samp?file=' + encodeURIComponent(ajs) + '&_t=' + Date.now(); var sink=document.getElementsByName('aladin_sink')[0]; if(sink){ try{ sink.src = url; return true; }catch(e){} } try{ window.open(url, '_blank'); return true; }catch(e){} } }catch(e){}\n"
                    "    return false;\n"
                    "  }\n"
                    "  function onKey(e){\n"
                    "    if(isEditable(e.target)) return;\n"
                    "    var k = (e.key||'');\n"
                    "    if(k === 'ArrowLeft'){ var a = findNav('prev'); if(a && a.getAttribute('href')){ e.preventDefault(); location.href = a.getAttribute('href'); } }\n"
                    "    else if(k === 'ArrowRight'){ var b = findNav('next'); if(b && b.getAttribute('href')){ e.preventDefault(); location.href = b.getAttribute('href'); } }\n"
                    "    else if(k && k.toLowerCase() === 'a' && !e.ctrlKey && !e.metaKey && !e.altKey){ if(triggerAladin()){ e.preventDefault(); } }\n"
                    "  }\n"
                    "  document.addEventListener('keydown', onKey, true);\n"
                    "})();\n"
                )
                if not os.path.exists(nav_js_path):
                    with open(nav_js_path, 'w') as _f:
                        _f.write(nav_js)
            except Exception:
                pass
            page_num = (start // per_page) + 1
            # HTML emission mode: all | problems | none (env: ALADIN_HTML_MODE)
            try:
                _html_mode = str(os.environ.get('ALADIN_HTML_MODE', 'all')).strip().lower()
            except Exception:
                _html_mode = 'all'
            if _html_mode not in ('all', 'problems', 'problem', 'problem-only', 'none'):
                _html_mode = 'all'
            page_has_problem = any(bool(a.get('problem')) for a in page_areas) if page_areas else False
            skip_html = (_html_mode == 'none') or (_html_mode in ('problems', 'problem', 'problem-only') and not page_has_problem)

            # Choose raster format for HTML image (env: HTML_IMAGE_FORMAT = png|jpg|webp)
            try:
                _imgfmt = str(os.environ.get('HTML_IMAGE_FORMAT', 'jpg')).strip().lower()
            except Exception:
                _imgfmt = 'jpg'
            if _imgfmt in ('jpg', 'jpeg'):
                _ext = 'jpg'
                _imgfmt = 'jpg'
            elif _imgfmt in ('webp',):
                _ext = 'webp'
                _imgfmt = 'webp'
            else:
                _ext = 'png'
                _imgfmt = 'png'
            png_name = f'{key}_page{page_num}.{_ext}'
            png_path = os.path.join(html_dir, png_name)
            html_name = f'{key}_page{page_num}.html'
            html_path = os.path.join(html_dir, html_name)
            # Save figure to raster at a consistent DPI
            HTML_DPI = int(str(os.environ.get('HTML_DPI', '150')))
            # Display scale for HTML (e.g., 0.7 makes it 30% smaller)
            try:
                HTML_SCALE = float(str(os.environ.get('HTML_SCALE', '0.7')))
            except Exception:
                HTML_SCALE = 0.7
            HTML_SCALE = max(0.2, min(2.0, HTML_SCALE))
            # Additional raster downscale factor relative to on-page display size.
            # For example, 0.5 saves the image at half the displayed resolution
            # (browser will upscale to display size), greatly reducing file size.
            try:
                HTML_RASTER_SCALE = float(str(os.environ.get('HTML_RASTER_SCALE', '0.8')))
            except Exception:
                HTML_RASTER_SCALE = 0.8
            HTML_RASTER_SCALE = max(0.2, min(1.0, HTML_RASTER_SCALE))
            # Ensure layout is finalized first (so text bboxes match saved PNG), then save
            try:
                # Hide in-figure link labels for PNG (HTML uses overlay buttons)
                for _t in hidden_texts:
                    try:
                        _t.set_visible(False)
                    except Exception:
                        pass
                fig.canvas.draw()
            except Exception:
                pass
            # Match renderer DPI to output to avoid bbox scaling mismatch
            try:
                fig.set_dpi(HTML_DPI)
            except Exception:
                pass
            # Save raster (respect format/quality); skip if HTML is disabled for this page
            if not skip_html:
                try:
                    # Compute the actual save DPI to cap resolution:
                    # save_dpi = HTML_DPI * HTML_SCALE (display size) * HTML_RASTER_SCALE (extra downscale)
                    dpi_out = int(max(50, round(HTML_DPI * HTML_SCALE * HTML_RASTER_SCALE)))
                    if _imgfmt == 'jpg':
                        try:
                            _q = int(str(os.environ.get('HTML_JPEG_QUALITY', '80')))
                        except Exception:
                            _q = 80
                        _q = max(20, min(95, _q))
                        fig.savefig(png_path, dpi=dpi_out, format='jpg', pil_kwargs={'quality': _q, 'optimize': True, 'progressive': True})
                    elif _imgfmt == 'webp':
                        try:
                            _q = int(str(os.environ.get('HTML_WEBP_QUALITY', '80')))
                        except Exception:
                            _q = 80
                        _q = max(30, min(95, _q))
                        fig.savefig(png_path, dpi=dpi_out, format='webp', pil_kwargs={'quality': _q, 'method': 6})
                    else:
                        fig.savefig(png_path, dpi=dpi_out)
                except TypeError:
                    # Older Matplotlib without pil_kwargs support
                    if _imgfmt == 'jpg':
                        fig.savefig(png_path, dpi=dpi_out, format='jpg')
                    elif _imgfmt == 'webp':
                        fig.savefig(png_path, dpi=dpi_out, format='webp')
                    else:
                        fig.savefig(png_path, dpi=dpi_out)
            fig_w = int(fig.get_size_inches()[0] * HTML_DPI)
            fig_h = int(fig.get_size_inches()[1] * HTML_DPI)
            disp_w = int(fig_w * HTML_SCALE)
            disp_h = int(fig_h * HTML_SCALE)
            xscale = disp_w / max(1, fig_w)
            yscale = disp_h / max(1, fig_h)
            # Button size heuristics (pixels)
            BTN_W = int(110 * HTML_SCALE)
            BTN_H = int(24 * HTML_SCALE)
            areas = []
            renderer = getattr(fig.canvas, 'get_renderer', lambda: None)()
            for area in page_areas:
                # Compute overlays from axes window extents with fixed insets
                ax_index = area.get('ax_index')
                try:
                    ax = axes_flat[ax_index]
                except Exception:
                    continue
                try:
                    win = ax.get_window_extent(renderer=renderer)
                    x0p, y0p, x1p, y1p = float(win.x0), float(win.y0), float(win.x1), float(win.y1)
                except Exception:
                    pos = ax.get_position()
                    x0p, y0p, x1p, y1p = pos.x0 * fig_w, pos.y0 * fig_h, pos.x1 * fig_w, pos.y1 * fig_h
                ax_l = int(x0p * xscale)
                ax_r = int(x1p * xscale)
                ax_top = int((fig_h - y1p) * yscale)
                ax_bot = int((fig_h - y0p) * yscale)
                ax_w = max(0, ax_r - ax_l)
                ax_h = max(0, ax_bot - ax_top)
                # Use fixed-pixel inset (scaled) for consistent alignment; user-tunable via HTML_BTN_INSET
                try:
                    INSET_PX = int(str(os.environ.get('HTML_BTN_INSET', '16')))
                except Exception:
                    INSET_PX = 8
                inset = max(1, int(INSET_PX * HTML_SCALE))
                # Top-left button (Aladin)
                left_x1 = ax_l + inset
                left_y1 = ax_top + inset
                left_x2 = min(ax_r, left_x1 + BTN_W)
                left_y2 = min(ax_bot, left_y1 + BTN_H)
                # Top-right button (Aladin Lite)
                right_x2 = ax_r - inset
                right_x1 = max(ax_l, right_x2 - BTN_W)
                right_y1 = ax_top + inset
                right_y2 = min(ax_bot, right_y1 + BTN_H)
                al_rect = (left_x1, left_y1, left_x2, left_y2)
                lite_rect = (right_x1, right_y1, right_x2, right_y2)
                # Full panel rectangle (entire axes region)
                panel_rect = (ax_l, ax_top, ax_r, ax_bot)

                areas.append({
                    'base': area['base'],
                    'type': 'aladin',
                    'coords': al_rect,
                    'href': area['aladin'],
                })
                areas.append({
                    'base': area['base'],
                    'type': 'lite',
                    'coords': lite_rect,
                    'href': area['lite'],
                })
                # Add a panel overlay for skip/filter UI and problem flag
                try:
                    areas.append({
                        'base': area.get('base',''),
                        'type': 'panel',
                        'coords': panel_rect,
                        'href': '',
                        'problem': int(area.get('problem', 0)),
                    })
                except Exception:
                    pass
            # Write page HTML
            map_name = f"map_{key}_{page_num}"
            lines = [
                '<!doctype html>',
                '<meta charset="utf-8">',
                f'<title>{key} — Page {page_num}</title>',
                # Responsive layout: image scales to container width, overlays use % positions
                '<style>'
                'body{margin:10px;background:#fafafa}'
                f'.wrap{{position:relative;max-width:100%;width:{disp_w}px}}'
                'img{display:block;width:100%;height:auto;border:0;box-shadow:0 1px 6px rgba(0,0,0,.1)}'
                '.btn{position:absolute;display:flex;align-items:center;justify-content:center;padding:2px 6px;'
                'border-radius:4px;font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Arial;font-size:calc(10px + 0.2vw);'
                'color:#fff;text-decoration:none;cursor:pointer;z-index:10;opacity:.92}'
                '.btn:hover{opacity:1}'
                '.btn.al{background:#2ca02c}'
                '.btn.lite{background:#1f77b4}'
                '.nav{margin-bottom:8px;font:14px -apple-system,BlinkMacSystemFont,Segoe UI,Arial}'
                '</style>',
                # Navigation links
                '<div class="nav">'
            ]
            prev_name = f'{key}_page{page_num-1}.html' if page_num>1 else ''
            next_name = f'{key}_page{page_num+1}.html' if page_num<total_pages else ''
            if prev_name:
                lines.append(f'<a href="{prev_name}">« Prev Page</a>')
            if next_name:
                if prev_name:
                    lines.append(' | ')
                lines.append(f'<a href="{next_name}">Next Page »</a>')
            # Keyboard usage hint in the top nav + filter toggle
            lines.append('<span style="opacity:.7;margin-left:8px">Keyboard: ←/→ pages, "a" Aladin</span>')
            lines.append(' | <label><input type="checkbox" id="only_prob"> Show only to-check</label>')
            # Top-right link to go back one level (master index)
            lines.append('<span style="float:right"><a href="../../aladin_index.html" title="Back to index">Back</a></span>')
            # Jump-to-page form
            lines.append(
                f'<form class="jump" onsubmit="return gotoPageFromInput(this);" '
                f'style="display:inline;margin-left:10px">'
                f'<label>Page <input type="number" min="1" max="{total_pages}" value="{page_num}" '
                f'style="width:4em" /> of {total_pages}</label> '
                f'<button type="submit">Go</button>'
                f'</form>'
            )
            lines += [
                '</div>',
                # hidden sink frame so link navigations don't replace the page
                '<iframe name="aladin_sink" width="0" height="0" style="display:none;border:0;position:absolute;left:-9999px;"></iframe>',
                f'<div class="wrap"><img src="{png_name}" alt="">'
            ]
            # Overlay anchors
            aladin_items_js = []
            for a in areas:
                x1, y1, x2, y2 = a['coords']
                w = max(1, x2 - x1); h = max(1, y2 - y1)
                # Convert to percentages for responsive scaling
                lp = (x1 / disp_w) * 100.0
                tp = (y1 / disp_h) * 100.0
                wp = (w  / disp_w) * 100.0
                hp = (h  / disp_h) * 100.0
                href = a['href'] if a['href'] else ''
                if a['type'] == 'aladin':
                    base_ajs = a['base'] + '.ajs'
                    default_http = f'http://127.0.0.1:8765/run_samp?file={base_ajs}'
                    lines.append(
                        f'<a class="btn al" style="left:{lp:.3f}%;top:{tp:.3f}%;width:{wp:.3f}%;height:{hp:.3f}%;" '
                        f'href="{default_http}" target="aladin_sink" data-ajs="{base_ajs}" onclick="return runSamp(event,this);" title="Aladin">Aladin</a>'
                    )
                    aladin_items_js.append(base_ajs)
                elif a['type'] == 'lite':
                    lines.append(
                        f'<a class="btn lite" style="left:{lp:.3f}%;top:{tp:.3f}%;width:{wp:.3f}%;height:{hp:.3f}%;" '
                        f'href="{href}" title="Aladin Lite" target="_blank" rel="noopener">Aladin&nbsp;Lite</a>'
                    )
                elif a['type'] == 'panel':
                    # Add panel overlay container with skip checkbox and mask for filtering
                    base_id = a.get('base', '')
                    prob = int(a.get('problem', 0))
                    lines.append(
                        f'<div class="panelbox" data-base="{base_id}" data-problem="{prob}" '
                        f'style="position:absolute;left:{lp:.3f}%;top:{tp:.3f}%;width:{wp:.3f}%;height:{hp:.3f}%;z-index:8;">'
                        f'<input type="checkbox" class="skipbox" data-base="{base_id}" title="Skip" '
                        f'style="position:absolute;left:50%;top:4px;transform:translateX(-50%);z-index:12;" onclick="return toggleSkip(this);" />'
                        f'<div class="mask" style="display:none;position:absolute;left:0;top:0;width:100%;height:100%;background:rgba(255,255,255,0.85);"></div>'
                        f'</div>'
                    )
            lines.append('</div>')
            # JS helper to call local SAMP link server when available
            lines += [
                '<script>\n'
                f'const ALADIN_ITEMS = {aladin_items_js!s};\n'
                'const SAMP_URL = (window.localStorage && localStorage.getItem("ALADIN_LINK_URL")) || "http://127.0.0.1:8765";\n'
                f'const PAGE_KEY = {key!r}; const PAGE_NUM = {page_num}; const TOTAL_PAGES = {total_pages};\n'
                'try{ window.PAGE_KEY = PAGE_KEY; window.PAGE_NUM = PAGE_NUM; window.TOTAL_PAGES = TOTAL_PAGES; }catch(e){}\n'
                'function _baseUrl(u){ try{ return (u && u.endsWith && u.endsWith("/")) ? u.slice(0,-1) : u; }catch(e){ return u; } }\n'
                'function runSamp(ev, el){\n'
                '  try{\n'
                '    const ajs = el.getAttribute("data-ajs");\n'
                '    const base = _baseUrl(SAMP_URL);\n'
                '    const url = base + "/run_samp?file=" + encodeURIComponent(ajs) + "&_t=" + Date.now();\n'
                '    // Update link target/href to use the configured helper, let browser navigate in hidden frame\n'
                '    el.href = url;\n'
                '    el.target = "aladin_sink";\n'
                '  }catch(e){}\n'
                '  return true;\n'
                '}\n'
                'function gotoPage(n){ try{ n = parseInt(n,10)||1; }catch(e){ n=1;} if(n<1) n=1; if(n>TOTAL_PAGES) n=TOTAL_PAGES; window.location.href = PAGE_KEY + "_page" + n + ".html"; return false; }\n'
                'function gotoPageFromInput(form){ var inp = form.querySelector("input[type=number]"); if(!inp) return false; return gotoPage(inp.value); }\n'
                'async function nextSrc(dir){\n'
                '  if(!ALADIN_ITEMS.length) return false;\n'
                '  window._srcIdx = (window._srcIdx===undefined? -1 : window._srcIdx) + (dir||1);\n'
                '  if(window._srcIdx < 0) window._srcIdx = ALADIN_ITEMS.length-1;\n'
                '  if(window._srcIdx >= ALADIN_ITEMS.length) window._srcIdx = 0;\n'
                '  const ajs = ALADIN_ITEMS[window._srcIdx];\n'
                '  try{ const base = _baseUrl(SAMP_URL); const url = base + "/run_samp?file=" + encodeURIComponent(ajs); await fetch(url); }catch(e){}\n'
                '  return false;\n'
                '}\n'
                '// Panel filtering and skip controls\n'
                'function _readSkip(){ try{ const k = "SKIP_"+PAGE_KEY; const s = localStorage.getItem(k); return s? JSON.parse(s) : {}; }catch(e){ return {}; } }\n'
                'function _writeSkip(obj){ try{ const k = "SKIP_"+PAGE_KEY; localStorage.setItem(k, JSON.stringify(obj)); }catch(e){} }\n'
                'function _readOnly(){ try{ const k = "ONLY_PROB_"+PAGE_KEY; return localStorage.getItem(k)==='+'"1"'+'; }catch(e){ return false; } }\n'
                'function _writeOnly(v){ try{ const k = "ONLY_PROB_"+PAGE_KEY; localStorage.setItem(k, v?"1":"0"); }catch(e){} }\n'
                'function toggleSkip(cb){ try{ const base = cb.getAttribute("data-base"); if(!base) return false; const sk=_readSkip(); if(cb.checked){ sk[base]=1; } else { delete sk[base]; } _writeSkip(sk); const pb = cb.closest(".panelbox"); if(pb){ pb.setAttribute("data-skip", cb.checked?"1":"0"); } updatePanels(); }catch(e){} return true; }\n'
                'function updatePanels(){ try{ const only = document.getElementById("only_prob"); const onlyOn = !!(only && only.checked); const pb = document.querySelectorAll(".panelbox"); pb.forEach(function(p){ const prob = p.getAttribute("data-problem")==="1"; const skip = p.getAttribute("data-skip")==="1"; const mask = p.querySelector(".mask"); let show = true; if(onlyOn){ show = prob && !skip; } else { show = !skip; } if(mask){ mask.style.display = show? "none" : "block"; } }); }catch(e){} return false; }\n'
                'function _hasUnskippedProblems(){ try{ const pb = document.querySelectorAll(".panelbox"); for(let i=0;i<pb.length;i++){ const p=pb[i]; if(p.getAttribute("data-problem")==="1" && p.getAttribute("data-skip")!=="1"){ return true; } } }catch(e){} return false; }\n'
                'function _getSeek(){ try{ const h=String(location.hash||""); const m=h.match(/seek=(next|prev)/); return m? m[1] : "next"; }catch(e){ return "next"; } }\n'
                'function _navToPage(n, dir){ if(typeof TOTAL_PAGES!=="number") return; if(n<1||n>TOTAL_PAGES) return; const href = PAGE_KEY + "_page" + n + ".html#seek=" + (dir||"next"); location.href = href; }\n'
                'function _maybeAutoAdvance(){ try{ const onlyOn = _readOnly(); if(!onlyOn) return; if(_hasUnskippedProblems()) return; const dir = _getSeek(); if(dir==="prev"){ if(PAGE_NUM>1){ _navToPage(PAGE_NUM-1, "prev"); } } else { if(PAGE_NUM<TOTAL_PAGES){ _navToPage(PAGE_NUM+1, "next"); } } }catch(e){} }\n'
                'document.addEventListener("DOMContentLoaded", function(){ try{ const sk=_readSkip(); document.querySelectorAll(".panelbox").forEach(function(p){ const b=p.getAttribute("data-base"); if(sk[b]){ p.setAttribute("data-skip","1"); const cb=p.querySelector(".skipbox"); if(cb) cb.checked=true; } }); const only=document.getElementById("only_prob"); if(only){ only.checked = _readOnly(); only.addEventListener("change", function(){ _writeOnly(only.checked); updatePanels(); if(only.checked){ _navToPage(PAGE_NUM, _getSeek()); } }); } updatePanels(); // tag Prev/Next anchors with seek\n'
                '  document.querySelectorAll(".nav a").forEach(function(a){ const t=(a.textContent||"").toLowerCase(); if(t.indexOf("next page")!==-1){ if(a.href.indexOf("#seek=")===-1) a.href += "#seek=next"; } if(t.indexOf("prev page")!==-1 || t.indexOf("« prev page")!==-1){ if(a.href.indexOf("#seek=")===-1) a.href += "#seek=prev"; } }); _maybeAutoAdvance(); }catch(e){} });\n'
                '</script>'
            ]
            # Remove bottom per-source nav (was unreliable without local server)
            lines.append('<script src="nav_keys.js"></script>')
            # Only write HTML if not skipped
            if not skip_html:
                with open(html_path, 'w') as f:
                    f.write('\n'.join(lines) + '\n')
                page_files.append(os.path.join('html', html_name))
        except Exception as e:
            print(f"[html] Failed to write page HTML: {e}")

        if pdf:
            # Restore labels for PDF version
            try:
                for _t in hidden_texts:
                    try:
                        _t.set_visible(True)
                    except Exception:
                        pass
                fig.canvas.draw()
            except Exception:
                pass
            pdf.savefig(fig)
        plt.close(fig)
    # Close SAMP connection
    if samp_client is not None:
        try:
            samp_client.disconnect()
        except Exception:
            pass
    if pdf:
        pdf.close()

    # Skip writing per-catalog index pages; we'll link directly to page 1 from the master index


# Modify main() to return combined_df and catalogs for interactive use
def main() -> None:
    parser = argparse.ArgumentParser(description="Build NGC2264 multi-wavelength catalog")
    parser.add_argument('-r', '--refresh', action='store_true', default=False,
                        help='Re-download all catalogs, ignoring cache')
    parser.add_argument('--pdf', action='store_true', default=True,
                        help='Generate PDF of match ellipses per catalog')
    parser.add_argument('--grid', default='1x1',
                        help="Panels per PDF page as CxR, e.g. '3x2' for 3 columns × 2 rows")
    parser.add_argument('--invert-cmap', action='store_true', default=True,
                        help='Invert grayscale colormap for background images (2MASS J)')
    parser.add_argument('--no-images', action='store_true', default=False,
                        help='Skip downloading 2MASS J images in PDFs (faster, avoids network stalls)')
    parser.add_argument('--samp', action='store_true', default=False,
                        help='If an Aladin SAMP hub is running, send panel views to Aladin Desktop')
    parser.add_argument('--samp-addr', default='127.0.0.1',
                        help='Local address to bind the SAMP client to (127.0.0.1 or ::1)')
    parser.add_argument('--chandra-csv-path', default='/Users/ettoref/ASTRONOMY/DATA/N2264_XMM_alt/N2264_acis12.csv',
                        help='Path to pre-exported Chandra CSV with columns ra_deg,dec_deg,epos/pub_n')
    args, _ = parser.parse_known_args()
    global REFRESH
    REFRESH = args.refresh
    # Parse grid argument (columns x rows)
    try:
        _gc, _gr = args.grid.lower().split('x')
        GRID_COLS, GRID_ROWS = max(1, int(_gc)), max(1, int(_gr))
    except Exception:
        print("[warn] --grid format invalid; using default 2x2")
        GRID_COLS, GRID_ROWS = 2, 2
    global INVERT_CMAP
    INVERT_CMAP = args.invert_cmap
    # Per-catalog parameters: factor and min_radius (arcsec)
    catalog_params = {
        'gaia':    {'factor': 2.0, 'min_radius': 0.05, 'epoch': 2016.0},
        '2MASS':   {'factor': 2.0, 'min_radius': 0.10, 'epoch': 2000.0},
        'wise':    {'factor': 2.0, 'min_radius': 1.00, 'epoch': 2010.5},
        'chandra': {'factor': 2.0, 'min_radius': 5.00, 'epoch': 2002.0},
        'xmm':     {'factor': 2.0, 'min_radius': 5.00, 'epoch': 2003.0}
    }
    # Define the central coordinate of NGC 2264 and search radius
    center = SkyCoord(ra=100.25 * u.deg, dec=9.883333 * u.deg, frame='icrs')
    radius = 1.0 * u.deg

    print("Querying Gaia DR3 ...")
    gaia = query_gaia(center, radius)
    print(f"Retrieved {len(gaia)} Gaia sources (raw)")
    gaia = filter_gaia_quality(gaia)

    print("Querying 2MASS ...")
    tmass = query_2mass(center, radius)
    print(f"Retrieved {len(tmass)} 2MASS sources (raw)")
    tmass = filter_2mass_quality(tmass)

    print("Querying WISE ...")
    wise = query_wise(center, radius)
    print(f"Retrieved {len(wise)} WISE sources (raw)")
    wise = filter_wise_quality(wise)

    print("Querying Chandra from CSV ...")
    try:
        chan = query_chandra_csv(center, radius, csv_path=args.chandra_csv_path)
    except Exception as e:
        raise RuntimeError(f"Failed to load Chandra CSV ('{args.chandra_csv_path}'): {e}")
    print(f"Retrieved {len(chan)} Chandra sources (raw)")
    chan = filter_chandra_quality(chan)

    print("Querying XMM‑Newton ...")
    xmm = query_xmm(center, radius)
    print(f"Retrieved {len(xmm)} XMM sources (raw)")
    xmm = filter_xmm_quality(xmm)

    # Build master catalog cumulatively by matching and appending
    combined_df = gaia.copy()
    # Apply Gaia parameters to its error ellipses
    gp = catalog_params['gaia']
    combined_df['errMaj'] *= gp['factor']
    combined_df['errMin'] *= gp['factor']
    # Clip to minimum radius
    combined_df['errMaj'] = combined_df['errMaj'].clip(lower=gp['min_radius'])
    combined_df['errMin'] = combined_df['errMin'].clip(lower=gp['min_radius'])
    # Store minimal original catalogs for reference as dict of DataFrames
    catalogs: dict[str, object] = {}
    # Save Gaia catalog (indexed by gaia_id) as DataFrame with factor and min_radius
    # Build Gaia columns list defensively in case cache lacks PM fields
    gaia_cols_base = ['ra_deg','dec_deg','errMaj','errMin','errPA']
    gaia_cols_pm = ['pmra','pmdec','pmra_error','pmdec_error','ref_epoch']
    gaia_cols = gaia_cols_base + [c for c in gaia_cols_pm if c in gaia.columns]
    catalogs['gaia'] = {
        'data': gaia.set_index('gaia_id')[gaia_cols],
        **catalog_params['gaia']
    }

    for df_other, id_col in [
        (tmass,     '2MASS'),
        (wise,      'wise_id'),
        (chan,      'chandra_id'),
        (xmm,       'xmm_id'),
    ]:
        # Strip '_id' suffix for catalog key
        key = id_col.replace('_id', '')
        params = catalog_params[key]
        # Always include RA and Dec in each catalog
        cols: list[str] = ['ra_deg', 'dec_deg']
        for ellipse in ['errMaj', 'errMin', 'errPA']:
            if ellipse in df_other.columns:
                cols.append(ellipse)
        catalogs[key] = {
            'data': df_other.set_index(id_col)[cols],
            **params
        }
        print(f"Matching {id_col} with master catalog ...")
        # cross-match df_other to combined_df
        pdf_path = f"{key}_matches.pdf" if args.pdf else None
        matches = cross_match(
            combined_df, df_other, id_col,
            min_radius = params['min_radius'] * u.arcsec,
            pdf_path   = None,  # plot AFTER merge so panels reflect final combined_df
            scale_factor = params['factor'],
            all_catalogs = catalogs,
            plot_mode = 'match'
        )
        # Initialize this catalog's ID column with the same dtype as df_other[id_col].
        # Use a sensible sentinel: numeric dtypes get -1 and object/string dtypes get None.
        col_dtype = df_other[id_col].dtype
        sentinel: object
        # pandas dtypes have a ``kind`` attribute: ``i`` for signed int, ``u`` for unsigned,
        # ``f`` for float and ``O`` for object.  We treat numeric kinds uniformly.
        if getattr(col_dtype, "kind", None) in ("i", "u", "f"):
            sentinel = -1
        else:
            sentinel = None
        combined_df[id_col] = pd.Series([sentinel] * len(combined_df),
                                         index=combined_df.index,
                                         dtype=col_dtype)
        rows_to_append: list[dict] = []
        # Index df_other by its identifier for direct lookup
        df_other_indexed = df_other.set_index(id_col)
        # First, update existing rows for matches
        for i in combined_df.index:
            matched_ids = matches[i]
            if not matched_ids:
                continue
            # Handle first match in place
            first_oid = matched_ids[0]
            other_row = df_other_indexed.loc[first_oid]
            # Scale other ellipse and enforce minimum radius
            scaled_maj = max(other_row['errMaj'] * params['factor'],
                             params['min_radius'])
            scaled_min = max(other_row['errMin'] * params['factor'],
                             params['min_radius'])
            scaled_pa  = other_row['errPA']
            # Choose smaller ellipse
            if scaled_maj < combined_df.at[i, 'errMaj']:
                combined_df.at[i, 'ra_deg'] = other_row['ra_deg']
                combined_df.at[i, 'dec_deg'] = other_row['dec_deg']
                combined_df.at[i, 'errMaj']  = scaled_maj
                combined_df.at[i, 'errMin']  = scaled_min
                combined_df.at[i, 'errPA']   = scaled_pa
            combined_df.at[i, id_col] = first_oid
            # Duplicate for any additional matches
            for oid in matched_ids[1:]:
                # Convert lookup key to integer when possible
                try:
                    lookup_key = int(oid)
                except ValueError:
                    lookup_key = oid
                other_row = df_other_indexed.loc[lookup_key]
                # Scale other ellipse and enforce minimum radius
                scaled_maj = max(other_row['errMaj'] * params['factor'],
                                 params['min_radius'])
                scaled_min = max(other_row['errMin'] * params['factor'],
                                 params['min_radius'])
                scaled_pa  = other_row['errPA']
                # Determine which ellipse to use
                if scaled_maj < combined_df.at[i, 'errMaj']:
                    ra_sel, dec_sel = other_row['ra_deg'], other_row['dec_deg']
                    maj_sel, min_sel, pa_sel = scaled_maj, scaled_min, scaled_pa
                else:
                    ra_sel = combined_df.at[i, 'ra_deg']
                    dec_sel = combined_df.at[i, 'dec_deg']
                    maj_sel = combined_df.at[i, 'errMaj']
                    min_sel = combined_df.at[i, 'errMin']
                    pa_sel  = combined_df.at[i, 'errPA']
                # Create a new row dict from the current row
                base = combined_df.loc[i].to_dict()
                base.update({
                    'ra_deg': ra_sel,
                    'dec_deg': dec_sel,
                    'errMaj': maj_sel,
                    'errMin': min_sel,
                    'errPA':  pa_sel,
                    id_col:   oid
                })
                rows_to_append.append(base)
        # Then add sources in other_df with no matches
        rev_matches = cross_match(
            df_other, combined_df, id_col,
            min_radius = params['min_radius'] * u.arcsec,
            pdf_path   = None,
            scale_factor = params['factor'],
            all_catalogs = catalogs,
            plot_mode = 'match'
        )
        for j, ids in enumerate(rev_matches):
            if not ids:
                other_row = df_other.iloc[j]
                # Scale ellipse and enforce minimum radius
                maj   = max(other_row['errMaj'] * params['factor'],
                            params['min_radius'])
                min_  = max(other_row['errMin'] * params['factor'],
                            params['min_radius'])
                pa    = other_row['errPA']
                entry = {
                    'ra_deg': other_row['ra_deg'],
                    'dec_deg': other_row['dec_deg'],
                    'errMaj': maj,
                    'errMin': min_,
                    'errPA': pa,
                    id_col: other_row[id_col],
                    'gaia_id': -1
                }
                rows_to_append.append(entry)
        # Append all new rows to the existing combined_df
        if rows_to_append:
            combined_df = pd.concat([combined_df, pd.DataFrame(rows_to_append)], ignore_index=True)
        # For 2MASS, add a blend-aware association pass (no radius inflation)
        if id_col == '2MASS':
            combined_df = attach_2mass_blends(
                combined_df,
                df_other,
                r_group_arcsec=1.5,
                max_pair_sep_arcsec=2.0,
                r_centroid_arcsec=0.7,
                max_dG_mag=2.0,
                prefer_blflag=True
            )
            # Remove synthetic rows for this 2MASS id when a valid Gaia link exists
            combined_df = _dedupe_unmatched_for_id(combined_df, id_col)
        # Generate post-merge plots so panels reflect final combined_df
        if args.pdf:
            plot_after_merge(
                combined_df,
                df_other,
                id_col,
                catalogs,
                pdf_path=f"{key}_matches.pdf",
                plot_mode='match',
                ncols=GRID_COLS,
                nrows=GRID_ROWS,
                invert_cmap=INVERT_CMAP,
                draw_images=(not args.no_images),
                aladin_dir=os.environ.get('ALADIN_SCRIPTS_DIR', 'aladin_scripts'),
                samp_enabled=args.samp,
                samp_addr=args.samp_addr
            )
    # Keep only coordinate, error ellipse, and identification columns in the final catalog
    # id_columns = [
    #     'ra_deg', 'dec_deg',
    #     'errMaj', 'errMin', 'errPA',
    #     'gaia_id', '2MASS', 'wise_id', 'chandra_id', 'xmm_id'
    # ]
    # combined_df = combined_df[id_columns]
    combined_df.to_csv('ngc2264_combined.csv', index=False)
    print(f"Master catalog written to ngc2264_combined.csv ({len(combined_df)} rows)")

    # Write a master HTML index that links to the first existing page for each catalog
    try:
        scripts_dir = os.environ.get('ALADIN_SCRIPTS_DIR', 'aladin_scripts')
        first_pages = []
        html_dir = os.path.join(scripts_dir, 'html')
        if os.path.isdir(html_dir):
            try:
                import re
                first_map = {}
                for name in sorted(os.listdir(html_dir)):
                    m = re.match(r"^(.*)_page(\d+)\.html$", name)
                    if not m:
                        continue
                    prefix = m.group(1)
                    n = int(m.group(2))
                    if (prefix not in first_map) or (n < first_map[prefix][0]):
                        first_map[prefix] = (n, name)
                first_pages = [v[1] for k, v in sorted(first_map.items(), key=lambda kv: kv[0].lower())]
            except Exception:
                # Fallback to page 1 only
                for name in sorted(os.listdir(html_dir)):
                    if name.endswith('_page1.html'):
                        first_pages.append(name)
        if first_pages:
            lines = [
                '<!doctype html>',
                '<meta charset="utf-8">',
                '<title>Aladin Index</title>',
                '<style>body{font:14px -apple-system,BlinkMacSystemFont,Segoe UI,Arial} li{margin:6px 0}</style>',
                '<h2>Aladin Index</h2>',
                '<ul>'
            ]
            for p in first_pages:
                # label is the prefix before _page1.html
                lab = p.replace('_page1.html','')
                lines.append(f'<li><a href="{scripts_dir}/html/{p}">{lab}</a></li>')
            lines.append('</ul>')
            with open('aladin_index.html','w') as f:
                f.write('\n'.join(lines) + '\n')
            print('[index] Wrote aladin_index.html')
    except Exception as e:
        print(f"[index] Could not write master HTML index: {e}")
    return combined_df, catalogs


def build_data_for_web(refresh: bool = False,
                       chandra_csv_path: str = '/Users/ettoref/ASTRONOMY/DATA/N2264_XMM_alt/N2264_acis12.csv',
                       include_catalogs=None):
    """Prepare combined_df and catalogs for the dynamic web server without
    writing PDFs or HTML pages. Reuses the same query and matching logic as
    main(), but stops after producing the merged catalog in memory.

    Parameters
    ----------
    refresh : bool
        If True, ignore local caches when querying catalogs.
    chandra_csv_path : str
        Path to pre-exported Chandra CSV with columns ra_deg, dec_deg, epos/pub_n.

    Returns
    -------
    (combined_df, catalogs)
        The merged DataFrame and the per-catalog metadata dict used by plotting.
    """
    global REFRESH
    REFRESH = bool(refresh)
    # Per-catalog parameters: factor and min_radius (arcsec)
    catalog_params = {
        'gaia':    {'factor': 2.0, 'min_radius': 0.05, 'epoch': 2016.0},
        '2MASS':   {'factor': 2.0, 'min_radius': 0.10, 'epoch': 2000.0},
        'wise':    {'factor': 2.0, 'min_radius': 1.00, 'epoch': 2010.5},
        'chandra': {'factor': 2.0, 'min_radius': 5.00, 'epoch': 2002.0},
        'xmm':     {'factor': 2.0, 'min_radius': 5.00, 'epoch': 2003.0}
    }
    # Normalize include list
    include_set = None
    if include_catalogs is not None:
        try:
            if isinstance(include_catalogs, str):
                include_catalogs = [include_catalogs]
            include_set = set([str(x).strip().lower() for x in include_catalogs if str(x).strip()])
        except Exception:
            include_set = None
    center = SkyCoord(ra=100.25 * u.deg, dec=9.883333 * u.deg, frame='icrs')
    radius = 1.00 * u.deg

    gaia = query_gaia(center, radius)
    gaia = filter_gaia_quality(gaia)
    tmass = wise = chan = xmm = None
    if (include_set is None) or ('2mass' in include_set or '2mass_j' in include_set or '2mass-j' in include_set or '2MASS' in (include_catalogs if include_catalogs else [])):
        tmass = query_2mass(center, radius)
        tmass = filter_2mass_quality(tmass)
    if (include_set is None) or ('wise' in include_set):
        wise = query_wise(center, radius)
        wise = filter_wise_quality(wise)
    if (include_set is None) or ('chandra' in include_set):
        chan = query_chandra_csv(center, radius, csv_path=chandra_csv_path)
        chan = filter_chandra_quality(chan)
    if (include_set is None) or ('xmm' in include_set):
        xmm  = query_xmm(center, radius)
        xmm  = filter_xmm_quality(xmm)

    combined_df = gaia.copy()
    gp = catalog_params['gaia']
    combined_df['errMaj'] *= gp['factor']
    combined_df['errMin'] *= gp['factor']
    combined_df['errMaj'] = combined_df['errMaj'].clip(lower=gp['min_radius'])
    combined_df['errMin'] = combined_df['errMin'].clip(lower=gp['min_radius'])

    catalogs: dict[str, object] = {}
    gaia_cols_base = ['ra_deg','dec_deg','errMaj','errMin','errPA']
    gaia_cols_pm = ['pmra','pmdec','pmra_error','pmdec_error','ref_epoch']
    gaia_cols = gaia_cols_base + [c for c in gaia_cols_pm if c in gaia.columns]
    catalogs['gaia'] = {
        'data': gaia.set_index('gaia_id')[gaia_cols],
        'frame': gaia.copy(),
        **catalog_params['gaia']
    }

    for df_other, id_col in [
        (tmass,     '2MASS'),
        (wise,      'wise_id'),
        (chan,      'chandra_id'),
        (xmm,       'xmm_id'),
    ]:
        if df_other is None:
            continue
        key = id_col.replace('_id', '')
        params = catalog_params[key]
        cols: list[str] = ['ra_deg', 'dec_deg']
        for ellipse in ['errMaj', 'errMin', 'errPA']:
            if ellipse in df_other.columns:
                cols.append(ellipse)
        catalogs[key] = {
            'data': df_other.set_index(id_col)[cols],
            'frame': df_other.copy(),
            **params
        }
        # Cross-match and merge (same logic as main, minus plotting/PDF)
        matches = cross_match(
            combined_df, df_other, id_col,
            min_radius = params['min_radius'] * u.arcsec,
            pdf_path   = None,
            scale_factor = params['factor'],
            all_catalogs = catalogs,
            plot_mode = 'match'
        )
        # Initialize id column with the same dtype and sentinel
        col_dtype = df_other[id_col].dtype
        sentinel: object
        if getattr(col_dtype, "kind", None) in ("i", "u", "f"):
            sentinel = -1
        else:
            sentinel = None
        combined_df[id_col] = pd.Series([sentinel] * len(combined_df),
                                         index=combined_df.index,
                                         dtype=col_dtype)
        rows_to_append: list[dict] = []
        df_other_indexed = df_other.set_index(id_col)
        for i in combined_df.index:
            matched_ids = matches[i]
            if not matched_ids:
                continue
            first_oid = matched_ids[0]
            other_row = df_other_indexed.loc[first_oid]
            scaled_maj = max(other_row['errMaj'] * params['factor'], params['min_radius'])
            scaled_min = max(other_row['errMin'] * params['factor'], params['min_radius'])
            scaled_pa  = other_row['errPA']
            if scaled_maj < combined_df.at[i, 'errMaj']:
                combined_df.at[i, 'ra_deg'] = other_row['ra_deg']
                combined_df.at[i, 'dec_deg'] = other_row['dec_deg']
                combined_df.at[i, 'errMaj']  = scaled_maj
                combined_df.at[i, 'errMin']  = scaled_min
                combined_df.at[i, 'errPA']   = scaled_pa
            combined_df.at[i, id_col] = first_oid
            for oid in matched_ids[1:]:
                try:
                    lookup_key = int(oid)
                except ValueError:
                    lookup_key = oid
                other_row = df_other_indexed.loc[lookup_key]
                scaled_maj = max(other_row['errMaj'] * params['factor'], params['min_radius'])
                scaled_min = max(other_row['errMin'] * params['factor'], params['min_radius'])
                scaled_pa  = other_row['errPA']
                entry = {
                    'ra_deg': other_row['ra_deg'],
                    'dec_deg': other_row['dec_deg'],
                    'errMaj': min(combined_df.at[i, 'errMaj'], scaled_maj),
                    'errMin': min(combined_df.at[i, 'errMin'], scaled_min),
                    'errPA':  scaled_pa,
                }
                entry.update({c: combined_df.at[i, c] for c in ['gaia_id','2MASS','wise_id','chandra_id','xmm_id'] if c in combined_df.columns})
                entry[id_col] = lookup_key
                rows_to_append.append(entry)
        # Append synthetic rows for unmatched others
        rev_matches = cross_match(
            df_other, combined_df, id_col,
            min_radius = params['min_radius'] * u.arcsec,
            pdf_path   = None,
            scale_factor = params['factor'],
            all_catalogs = catalogs,
            plot_mode = 'match'
        )
        for j, ids in enumerate(rev_matches):
            if not ids:
                other_row = df_other.iloc[j]
                maj = max(other_row['errMaj'] * params['factor'], params['min_radius'])
                min_ = max(other_row['errMin'] * params['factor'], params['min_radius'])
                pa = other_row['errPA']
                entry = {
                    'ra_deg': other_row['ra_deg'],
                    'dec_deg': other_row['dec_deg'],
                    'errMaj': maj,
                    'errMin': min_,
                    'errPA': pa,
                    id_col: other_row[id_col],
                    'gaia_id': -1
                }
                rows_to_append.append(entry)
        if rows_to_append:
            combined_df = pd.concat([combined_df, pd.DataFrame(rows_to_append)], ignore_index=True)
        if id_col == '2MASS':
            combined_df = attach_2mass_blends(
                combined_df, df_other,
                r_group_arcsec=1.5,
                max_pair_sep_arcsec=2.0,
                r_centroid_arcsec=0.7,
                max_dG_mag=2.0,
                prefer_blflag=True
            )
            combined_df = _dedupe_unmatched_for_id(combined_df, id_col)

    return combined_df, catalogs


if __name__ == '__main__':
    combined_df, catalogs = main()
