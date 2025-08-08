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

"""

from __future__ import annotations

import sys
from typing import List, Optional
import numpy as np
import numpy.linalg as la

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages

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
    if not REFRESH and os.path.exists(cache_file):
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
    if not REFRESH and os.path.exists(cache_file):
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
    if not REFRESH and os.path.exists(cache_file):
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
    if not REFRESH and os.path.exists(cache_file):
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
    if not REFRESH and os.path.exists(cache_file):
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
                scale_factor: float = 2.0) -> List[List]:
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
    max_sep : Quantity
        Maximum separation allowed between objects; matches outside this
        radius are ignored.

    Returns
    -------
    list of lists
        Each element is a list of identifiers (strings) matched to the
        corresponding row in `source_df`.  Empty lists indicate no match.
    """
    if other_df.empty:
        return [[] for _ in range(len(source_df))]

    pdf = PdfPages(pdf_path) if pdf_path else None

    src_coord = SkyCoord(source_df['ra_deg'].values * u.deg,
                         source_df['dec_deg'].values * u.deg)
    oth_coord = SkyCoord(other_df['ra_deg'].values * u.deg,
                         other_df['dec_deg'].values * u.deg)

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

    if pdf:
        # Grid settings: 2 columns × 2 rows per PDF page
        ncols, nrows = 2, 2
        per_page = ncols * nrows
        total = len(other_df)
        for start in range(0, total, per_page):
            end = min(start + per_page, total)
            indices = list(range(start, end))
            fig, axes = plt.subplots(nrows, ncols, figsize=(11, 8.5))
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
                # Determine plotting margin in arcsec
                margin = max(other_df.iloc[j]['errMaj'], other_df.iloc[j]['errMin']) * 2
                # Enforce at least ±1 arcsec half‐width (5x5" area)
                margin = max(margin, 2.5)
                # Plot all new-catalog ellipses in light blue
                for k in range(len(other_df)):
                    if abs(dx_oth.iloc[k]) <= margin and abs(dy_oth.iloc[k]) <= margin:
                        e = Ellipse((dx_oth.iloc[k], dy_oth.iloc[k]),
                                    width=2*other_df.iloc[k]['errMaj'],
                                    height=2*other_df.iloc[k]['errMin'],
                                    angle=other_df.iloc[k]['errPA'],
                                    edgecolor='lightblue', linestyle='-', fill=False)
                        ax.add_patch(e)
                # Plot all master-catalog ellipses in light red
                for m in range(len(source_df)):
                    if abs(dx_src.iloc[m]) <= margin and abs(dy_src.iloc[m]) <= margin:
                        e = Ellipse((dx_src.iloc[m], dy_src.iloc[m]),
                                    width=2*source_df.iloc[m]['errMaj'],
                                    height=2*source_df.iloc[m]['errMin'],
                                    angle=source_df.iloc[m]['errPA'],
                                    edgecolor='lightcoral', linestyle='-', fill=False)
                        ax.add_patch(e)
                # Overlay and highlight the specific new source (blue solid)
                ell_new = Ellipse((0, 0),
                                  width=2*other_df.iloc[j]['errMaj'],
                                  height=2*other_df.iloc[j]['errMin'],
                                  angle=other_df.iloc[j]['errPA'],
                                  edgecolor='blue', linestyle='-',
                                  fill=False, label=f"new: {other_id_col.replace('_id','')}")
                ax.add_patch(ell_new)
                # Overlay and highlight matched master ellipses (red solid)
                plotted_master = False
                for i, ids in enumerate(matches):
                    if str(other_df.iloc[j][other_id_col]) in ids:
                        dx1 = dx_src.iloc[i]
                        dy1 = dy_src.iloc[i]
                        ell_m = Ellipse((dx1, dy1),
                                        width=2*source_df.iloc[i]['errMaj'],
                                        height=2*source_df.iloc[i]['errMin'],
                                        angle=source_df.iloc[i]['errPA'],
                                        edgecolor='red', linestyle='-',
                                        fill=False,
                                        label='master' if not plotted_master else None)
                        ax.add_patch(ell_m)
                        plotted_master = True
                # Configure axes
                ax.set_xlim(-margin, margin)
                ax.set_ylim(-margin, margin)
                ax.set_aspect('equal', 'box')
                ax.set_xlabel('ΔRA [arcsec]')
                ax.set_ylabel('ΔDec [arcsec]')
                ax.legend(loc='lower right')
                # Annotate IDs
                new_id = other_df.iloc[j][other_id_col]
                master_idx = next((i for i, ids in enumerate(matches) if str(new_id) in ids), None)
                if master_idx is not None and not pd.isna(source_df.iloc[master_idx]['gaia_id']):
                    master_val = source_df['gaia_id'].iat[master_idx]
                    try:
                        master_id = str(int(master_val))
                    except (ValueError, TypeError):
                        master_id = 'None'
                else:
                    master_id = 'None'
                ax.text(0.05, 0.95,
                        f"{other_id_col.replace('_id','')}: {new_id}\nmaster_gaia: {master_id}",
                        transform=ax.transAxes, va='top', ha='left')
                ax.set_title(f"Source {new_id}")
            pdf.savefig(fig)
            plt.close(fig)
        pdf.close()

    return matches


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
        'wise':    {'factor': 2.0, 'min_radius': 0.50},
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
        'data': gaia.set_index('gaia_id')[[
            'ra_deg', 'dec_deg',
            'errMaj', 'errMin', 'errPA'
        ]],
        **catalog_params['gaia']
    }

    for df_other, id_col in [
        (tmass,     '2MASS'),
        # (wise,      'wise_id'),
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
            pdf_path   = pdf_path,
            scale_factor = params['factor']
        )
        # Initialize this catalog's ID column with the same dtype as df_other[id_col]
        col_dtype = df_other[id_col].dtype
        combined_df[id_col] = pd.Series([-1] * len(combined_df),
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
            scale_factor = params['factor']
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