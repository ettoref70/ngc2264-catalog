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
except ImportError as exc:
    print("Missing required packages: {}".format(exc))
    print("Install dependencies with: pip install astropy astroquery pandas")
    sys.exit(1)


import os
# Directory to cache downloaded catalogs
CACHE_DIR = "catalogs"
os.makedirs(CACHE_DIR, exist_ok=True)

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
           source_id, phot_g_mean_mag
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
        'source_id', 'phot_g_mean_mag'
    ]]
    # Rename and cast Gaia ID to integer
    df = df.rename(columns={
        'ra': 'ra_deg',
        'dec': 'dec_deg',
        'ra_error': 'ra_err_mas',
        'dec_error': 'dec_err_mas',
        'ra_dec_corr': 'ra_dec_corr',
        'source_id': 'gaia_id',
        'phot_g_mean_mag': 'Gmag'
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
    df = df.drop(columns=[
        'ra_err_mas', 'dec_err_mas',
        'ra_err_deg', 'dec_err_deg',
        'ra_dec_corr'
    ])
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
    v = Vizier(
        columns=[
            'RAJ2000', 'DEJ2000',
            'errMaj', 'errMin', 'errPA',
            'Jmag', 'Hmag', 'Kmag', '2MASS'
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
    df = df.drop(columns=['ra_err_deg', 'dec_err_deg', 'pos_err_maj_arcsec', 'pos_err_min_arcsec', 'pos_err_pa_deg'])
    df.to_csv(cache_file, index=False)
    return df[[
        'ra_deg', 'dec_deg',
        'errMaj', 'errMin', 'errPA',
        'Jmag', 'Hmag', 'Kmag', '2MASS'
    ]]


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
        columns=['RAJ2000', 'DEJ2000', 'W1mpro', 'W2mpro', 'W3mpro', 'W4mpro',
                 'AllWISE', 'eeMaj', 'eeMin', 'eePA'],
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
    # Drop per-axis uncertainty and raw positional error columns
    df = df.drop(columns=['ra_err_deg', 'dec_err_deg', 'pos_err_maj_arcsec', 'pos_err_min_arcsec', 'pos_err_pa_deg'])
    df.to_csv(cache_file, index=False)
    return df[final_cols]


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
    # Drop per-axis uncertainty and raw positional error columns
    tab = tab.drop(columns=['ra_err_deg', 'dec_err_deg', 'pos_err_arcsec'])
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
    # Drop per-axis uncertainty and raw positional error columns
    tab = tab.drop(columns=['ra_err_deg', 'dec_err_deg', 'pos_err_arcsec'])
    tab.to_csv(cache_file, index=False)
    return tab[[
        'ra_deg', 'dec_deg',
        'errMaj', 'errMin', 'errPA',
        'xmm_id'
    ]]


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

    # Precompute 2x2 covariance matrices for source objects
    # Apply scale factor to ellipse axes before converting to degrees
    min_deg = min_radius.to(u.deg).value
    src_maj = np.maximum((source_df['errMaj'] * scale_factor).values / 3600.0, min_deg)  # deg
    src_min = np.maximum((source_df['errMin'] * scale_factor).values / 3600.0, min_deg)  # deg
    src_pa  = np.deg2rad(source_df['errPA'].values)  # radians
    n_src = len(source_df)
    src_cov = np.zeros((n_src, 2, 2))
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
        src_cov[i] = a*a * np.outer(vmaj, vmaj) + b*b * np.outer(vmin, vmin)

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

        ids = []
        for j in range(len(other_df)):
            # Combined covariance matrix
            C = src_cov[i] + oth_cov[j]
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
            vec = np.array([dra[j], ddec[j]])
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
            axes_flat = axes.flatten()
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


# === Post-merge plotting so panels reflect final combined_df ===
def plot_after_merge(combined_df: pd.DataFrame,
                     other_df: pd.DataFrame,
                     id_col: str,
                     all_catalogs: dict,
                     pdf_path: str,
                     plot_mode: str = 'match') -> None:
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
    ncols, nrows = 2, 2
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
        axes_flat = axes.flatten()
        for ax in axes_flat[len(indices):]:
            fig.delaxes(ax)
        for ax_idx, j in enumerate(indices):
            ax = axes_flat[ax_idx]
            new_id = other_df.iloc[j][id_col]
            ra0 = other_df.iloc[j]['ra_deg']
            dec0 = other_df.iloc[j]['dec_deg']

            # Offsets for background
            dx_oth = (other_df['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
            dy_oth = (other_df['dec_deg'] - dec0) * 3600.0
            dx_src = (combined_df['ra_deg'] - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
            dy_src = (combined_df['dec_deg'] - dec0) * 3600.0

            # Determine plotting margin in arcsec
            margin = max(oth_maj_plot[j], oth_min_plot[j]) * 2.0
            margin = max(margin, 2.5)

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
                            ell_m = Ellipse((dxm, dym), width=2*_maj_m, height=2*_min_m,
                                            angle=crow['errPA'], edgecolor=CAT_COLORS.get(kcat, 'red'), linestyle='-',
                                            linewidth=1.2, fill=False, label=None)
                            ax.add_patch(ell_m)
                            drawn_components.add((kcat, lookup_key))

            # Background: draw components for UNMATCHED master rows in lighter/transparent colors
            # and draw a thin gray cross at every master position
            matched_set = set(master_indices)
            for m in cand_indices:
                # gray cross at master position
                ax.plot(dx_src.iloc[m], dy_src.iloc[m], marker='+', color='gray', markersize=5,
                        markeredgewidth=0.6, linestyle='None')
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
            matched_set = set(master_indices)
            cand_with_dist = [(int(m), float(np.hypot(dx_src.iloc[m], dy_src.iloc[m]))) for m in cand_indices]
            cand_with_dist.sort(key=lambda t: t[1])

            def _fmt(val):
                return str(val) if val is not None else '—'

            headers = [k.upper() for k in display_cats]
            raw_rows = []  # list of (cells_list, is_matched)
            for m, _dist in cand_with_dist:
                mrow = combined_df.iloc[m]
                cells = []
                for kcat in display_cats:
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
                # Layout parameters
                y0 = 0.95
                line_h = 0.055  # axes fraction per line
                font_size = 9
                char_w_pts = font_size * 0.6  # approximate monospace char width in points
                # Compute x-offsets in points for each column (include two spaces between columns)
                x_offsets = []
                acc = 0.0
                for i in range(len(headers)):
                    x_offsets.append(acc)
                    acc += (col_w[i] + 2) * char_w_pts
                # Draw header per column with catalog colors
                for i, head in enumerate(headers):
                    trans = offset_copy(ax.transAxes, fig=fig, x=x_offsets[i], y=0, units='points')
                    ax.text(0.05, y0, head, transform=trans,
                            va='top', ha='left', family='monospace', fontsize=font_size,
                            fontweight='bold', color=CAT_COLORS.get(display_cats[i], 'black'))
                # Draw separator as dashes per column (neutral color)
                for i in range(len(headers)):
                    trans = offset_copy(ax.transAxes, fig=fig, x=x_offsets[i], y=0, units='points')
                    ax.text(0.05, y0 - line_h, '-' * col_w[i], transform=trans,
                            va='top', ha='left', family='monospace', fontsize=font_size)
                # Draw rows per column, color-coded and bold if matched
                for r_idx, (cells, is_matched) in enumerate(raw_rows):
                    y = y0 - (r_idx + 2) * line_h
                    for i, cell in enumerate(cells):
                        trans = offset_copy(ax.transAxes, fig=fig, x=x_offsets[i], y=0, units='points')
                        ax.text(0.05, y, cell.ljust(col_w[i]), transform=trans,
                                va='top', ha='left', family='monospace', fontsize=font_size,
                                fontweight=('bold' if is_matched else 'normal'),
                                color=CAT_COLORS.get(display_cats[i], 'black'))

            # Title with short id for the current catalog
            short_new = _short_of(key, new_id)
            if short_new is not None:
                ax.set_title(f"{key} #{short_new}", pad=10)
            else:
                ax.set_title(f"{key}", pad=10)

            ax.set_xlim(-margin, margin)
            ax.set_ylim(-margin, margin)
            ax.set_aspect('equal', 'box')
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
    args, _ = parser.parse_known_args()
    global REFRESH
    REFRESH = args.refresh
    # Per-catalog parameters: factor and min_radius (arcsec)
    catalog_params = {
        'gaia':    {'factor': 2.0, 'min_radius': 0.05},
        '2MASS':   {'factor': 2.0, 'min_radius': 0.10},
        'wise':    {'factor': 2.0, 'min_radius': 1.00},
        'chandra': {'factor': 2.0, 'min_radius': 5.00},
        'xmm':     {'factor': 2.0, 'min_radius': 5.00}
    }
    # Define the central coordinate of NGC 2264 and search radius
    center = SkyCoord(ra=100.25 * u.deg, dec=9.883333 * u.deg, frame='icrs')
    radius = 0.1 * u.deg

    print("Querying Gaia DR3 ...")
    gaia = query_gaia(center, radius)
    print(f"Retrieved {len(gaia)} Gaia sources")

    print("Querying 2MASS ...")
    tmass = query_2mass(center, radius)
    print(f"Retrieved {len(tmass)} 2MASS sources")

    print("Querying WISE ...")
    wise = query_wise(center, radius)
    print(f"Retrieved {len(wise)} WISE sources")

    print("Querying Chandra ...")
    chan = query_chandra(center, radius)
    print(f"Retrieved {len(chan)} Chandra sources")

    print("Querying XMM‑Newton ...")
    xmm = query_xmm(center, radius)
    print(f"Retrieved {len(xmm)} XMM sources")

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
        'data': gaia.set_index('gaia_id')[['ra_deg','dec_deg','errMaj','errMin','errPA']],
        **catalog_params['gaia']
    }

    for df_other, id_col in [
        (tmass,     '2MASS'),
        (wise,      'wise_id'),
        # (chan,      'chandra_id'),
        # (xmm,       'xmm_id'),
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
        # Generate post-merge plots so panels reflect final combined_df
        if args.pdf:
            plot_after_merge(
                combined_df,
                df_other,
                id_col,
                catalogs,
                pdf_path=f"{key}_matches.pdf",
                plot_mode='match'
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