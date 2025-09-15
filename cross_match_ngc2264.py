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

    matches = []
    for i in range(len(source_df)):
        # Offsets in degrees: RA scaled by cos(dec) for sky projection
        dec0 = np.deg2rad(source_df.iloc[i]['dec_deg'])
        ra0  = source_df.iloc[i]['ra_deg']
        dra  = (other_df['ra_deg'].values - ra0) * np.cos(dec0)
        ddec = other_df['dec_deg'].values - source_df.iloc[i]['dec_deg']

        # Apply proper-motion shift of the source position to the per-row target epoch (if available)
        dx_pm = 0.0
        dy_pm = 0.0
        epoch_j = _epoch_of_row(j)
        if (epoch_j is not None) and has_pm and has_epoch:
            try:
                ref_ep = float(source_df.iloc[i]['ref_epoch'])
            except Exception:
                ref_ep = None
            if ref_ep is not None and not np.isnan(ref_ep):
                dt = float(epoch_j - ref_ep)
                try:
                    pmra_val = float(source_df.iloc[i]['pmra'])
                except Exception:
                    pmra_val = 0.0
                try:
                    pmdec_val = float(source_df.iloc[i]['pmdec'])
                except Exception:
                    pmdec_val = 0.0
                # pmra is mu_alpha* (RA*cosDec) in mas/yr; convert to deg shift along RA* axis
                dx_pm = (pmra_val / 1000.0 / 3600.0) * dt
                dy_pm = (pmdec_val / 1000.0 / 3600.0) * dt

        ids = []
        for j in range(len(other_df)):
            # Combined covariance matrix
            C = src_cov[i] + oth_cov[j]
            # Add PM-uncertainty inflation for the baseline to this row's epoch
            if (epoch_j is not None) and has_pm_err and has_epoch:
                try:
                    ref_ep = float(source_df.iloc[i]['ref_epoch'])
                except Exception:
                    ref_ep = None
                if ref_ep is not None and not np.isnan(ref_ep):
                    dt = float(epoch_j - ref_ep)
                    try:
                        pmra_e = float(source_df.iloc[i]['pmra_error'])
                    except Exception:
                        pmra_e = 0.0
                    try:
                        pmdec_e = float(source_df.iloc[i]['pmdec_error'])
                    except Exception:
                        pmdec_e = 0.0
                    sig_ra_deg = abs(pmra_e) / 1000.0 / 3600.0 * abs(dt)
                    sig_de_deg = abs(pmdec_e) / 1000.0 / 3600.0 * abs(dt)
                    C = C + np.diag([sig_ra_deg**2, sig_de_deg**2])
            # Invert covariance, with fallback regularization if singular
            try:
                invC = la.inv(C)
            except la.LinAlgError:
                jitter = np.eye(2) * 1e-12
                try:
                    invC = la.inv(C + jitter)
                except la.LinAlgError:
                    # Skip this pair if still singular
                    continue
            # Subtract the source PM shift so comparison is at the catalog row's epoch
            vec = np.array([dra[j] - dx_pm, ddec[j] - dy_pm])
            # Mahalanobis squared distance
            D2 = vec.dot(invC).dot(vec)
            # Only match if within 1-sigma ellipse
            if D2 <= 1.0:
                # Preserve the original ID type: use column value if present, else index
                if other_id_col in other_df.columns:
                    ids.append(other_df.iloc[j][other_id_col])
                else:
                    ids.append(other_df.index[j])
        matches.append(ids)

    # --- Precompute reverse assignment (other -> master) with the same Mahalanobis test ---
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
                                    width=2*oth_maj_plot[k],
                                    height=2*oth_min_plot[k],
                                    angle=other_df.iloc[k]['errPA'],
                                    edgecolor='lightblue', linestyle='-', fill=False)
                        ax.add_patch(e)
                # Plot all master-catalog ellipses in light red
                for m in range(len(source_df)):
                    if abs(dx_src.iloc[m]) <= margin and abs(dy_src.iloc[m]) <= margin:
                        e = Ellipse((dx_src.iloc[m], dy_src.iloc[m]),
                                    width=2*src_maj_plot[m],
                                    height=2*src_min_plot[m],
                                    angle=source_df.iloc[m]['errPA'],
                                    edgecolor='lightcoral', linestyle='-', fill=False)
                        ax.add_patch(e)
                # Overlay and highlight the specific new source (blue solid)
                ell_new = Ellipse((0, 0),
                                  width=2*oth_maj_plot[j],
                                  height=2*oth_min_plot[j],
                                  angle=other_df.iloc[j]['errPA'],
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
                                            width=2*_maj_m,
                                            height=2*_min_m,
                                            angle=crow['errPA'],
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
      - 'auto'     : try SkyView, then fall back to HiPS on failure (default)
      - 'hips'     : use HiPS directly (skip SkyView)
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
    mode = os.environ.get("TMASS_FETCH", "auto").lower().strip()

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
                     draw_images: bool = True) -> None:
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
        pa = float((np.degrees(np.arctan2(v[0], v[1])) % 360.0))
        return maj, min_, pa

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

    # Precompute new-catalog ellipse axes for plotting
    if plot_mode == 'match':
        oth_maj_plot = np.maximum(other_df['errMaj'].values * factor, min_r)
        oth_min_plot = np.maximum(other_df['errMin'].values * factor, min_r)
    else:
        oth_maj_plot = other_df['errMaj'].values
        oth_min_plot = other_df['errMin'].values

    for start in range(0, total, per_page):
        end = min(start + per_page, total)
        indices = list(range(start, end))
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 8.5), constrained_layout=True)
        axes_flat = np.atleast_1d(axes).ravel()
        for ax in axes_flat[len(indices):]:
            fig.delaxes(ax)
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
                                width=2*oth_maj_plot[k],
                                height=2*oth_min_plot[k],
                                angle=other_df.iloc[k]['errPA'],
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
                         width=2*oth_maj_plot[j], height=2*oth_min_plot[j],
                         angle=other_df.iloc[j]['errPA'],
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
                            ell_m = Ellipse((dxm, dym), width=2*_maj_m, height=2*_min_m,
                                            angle=crow['errPA'], edgecolor=CAT_COLORS.get(kcat, 'red'), linestyle='-',
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
                                     angle=crow['errPA'], edgecolor=CAT_COLORS.get(kcat, 'red'), linestyle='-',
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

            # Build rows sorted by distance from center
            if 'gaia_id' in combined_df.columns:
                matched_set = {m for m in master_indices
                            if (not pd.isna(combined_df.iloc[m].get('gaia_id', np.nan))
                                and combined_df.iloc[m].get('gaia_id', -1) != -1)}
            else:
                matched_set = set(master_indices)
            cand_with_dist = [(int(m), float(np.hypot(dx_src.iloc[m], dy_src.iloc[m]))) for m in cand_indices]
            cand_with_dist.sort(key=lambda t: t[1])

            def _fmt(val):
                return str(val) if val is not None else '—'

            # Build table columns and headers; if GAIA is present, add a Gmag column right after it
            col_keys = list(display_cats)
            add_gmag = ('gaia' in display_cats) and ('Gmag' in combined_df.columns)
            if add_gmag:
                gi = col_keys.index('gaia')
                col_keys.insert(gi + 1, 'gaia_g')  # synthetic key for Gmag display
            headers = [
                ('G' if k == 'gaia_g' else k.upper())
                for k in col_keys
            ]
            raw_rows = []  # list of (cells_list, is_matched)
            for m, _dist in cand_with_dist:
                mrow = combined_df.iloc[m]
                cells = []
                for kcat in col_keys:
                    if kcat == 'gaia_g':
                        # Show Gaia G magnitude next to the GAIA id
                        if ('gaia_id' in combined_df.columns and not _is_missing(mrow.get('gaia_id', None))
                                and ('Gmag' in combined_df.columns) and pd.notna(mrow.get('Gmag', np.nan))):
                            try:
                                gval = float(mrow['Gmag'])
                                cells.append(f"{gval:.2f}")
                            except Exception:
                                cells.append(_fmt(None))
                        else:
                            cells.append(_fmt(None))
                        continue
                    colname = col_for[kcat]
                    if colname in combined_df.columns and not _is_missing(mrow.get(colname, None)):
                        cells.append(_fmt(_short_of(kcat, mrow[colname])))
                    else:
                        cells.append(_fmt(None))
                raw_rows.append((cells, (m in matched_set)))

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
                    key_color = CAT_COLORS.get('gaia' if key_i == 'gaia_g' else key_i, 'black')
                    ax.text(0.05, y0, head, transform=trans,
                            va='top', ha='left', family='monospace', fontsize=font_size,
                            fontweight='bold', color=key_color)

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
                        key_color = CAT_COLORS.get('gaia' if key_i == 'gaia_g' else key_i, 'black')
                        ax.text(0.05, y, text_cell.ljust(col_w[i]), transform=trans,
                                va='top', ha='left', family='monospace', fontsize=font_size,
                                fontweight=('bold' if is_matched else 'normal'),
                                color=key_color)

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
        if pdf:
            pdf.savefig(fig)
            plt.close(fig)
    if pdf:
        pdf.close()


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
    radius = 0.1 * u.deg

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
    catalogs['gaia'] = {
        'data': gaia.set_index('gaia_id')[['ra_deg','dec_deg','errMaj','errMin','errPA',
                                           'pmra','pmdec','pmra_error','pmdec_error','ref_epoch']],
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
                draw_images=(not args.no_images)
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
    return combined_df, catalogs


if __name__ == '__main__':
    combined_df, catalogs = main()
