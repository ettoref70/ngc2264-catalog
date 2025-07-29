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

try:
    from astroquery.gaia import Gaia
    from astroquery.vizier import Vizier
    from astroquery.irsa import Irsa
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
    # Use async cone search for efficiency
    job = Gaia.cone_search_async(center, radius)
    results = job.get_results()
    # Convert to DataFrame and keep relevant columns
    gaia_df = results.to_pandas()[['ra', 'dec', 'source_id', 'phot_g_mean_mag']]
    return gaia_df.rename(columns={
        'ra': 'ra_deg',
        'dec': 'dec_deg',
        'phot_g_mean_mag': 'Gmag'
    })


def query_2mass(center: SkyCoord, radius: u.Quantity) -> pd.DataFrame:
    """Query 2MASS Point Source Catalogue via Vizier.

    Returns a DataFrame with positions and photometry.
    """
    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'Jmag', 'Hmag', 'Kmag', '2MASS'],
               catalog='II/246', row_limit=-1)
    tables = v.query_region(center, radius=radius)
    if len(tables) == 0:
        return pd.DataFrame(columns=['ra_deg', 'dec_deg', 'Jmag', 'Hmag', 'Kmag', '2MASS'])
    tab = tables[0]
    df = tab.to_pandas()
    df = df.rename(columns={'RAJ2000': 'ra_deg', 'DEJ2000': 'dec_deg'})
    return df[['ra_deg', 'dec_deg', 'Jmag', 'Hmag', 'Kmag', '2MASS']]


def query_wise(center: SkyCoord, radius: u.Quantity) -> pd.DataFrame:
    """Query the AllWISE catalogue via IRSA.

    Returns a DataFrame with positions and photometry.
    """
    # IRSA expects radius in degrees; convert quantity to deg
    rad_deg = radius.to(u.deg).value
    wise_tab = Irsa.query_region(center, catalog='allwise_p3as_psc', spatial='Cone', radius=rad_deg)
    if wise_tab is None:
        return pd.DataFrame(columns=['ra_deg', 'dec_deg', 'w1mpro', 'w2mpro', 'w3mpro', 'w4mpro', 'allwise_id'])
    df = wise_tab.to_pandas()
    df = df.rename(columns={'ra': 'ra_deg', 'dec': 'dec_deg'})
    return df[['ra_deg', 'dec_deg', 'w1mpro', 'w2mpro', 'w3mpro', 'w4mpro', 'designation']]


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
    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'ACIS'], catalog='J/A+A/670/A37/tablea1', row_limit=-1)
    tables = v.query_region(center, radius=radius)
    if len(tables) == 0:
        return pd.DataFrame(columns=['ra_deg', 'dec_deg', 'chandra_id'])
    tab = tables[0]
    df = tab.to_pandas()
    # Rename columns to uniform names used in matching routines
    df = df.rename(columns={'RAJ2000': 'ra_deg', 'DEJ2000': 'dec_deg', 'ACIS': 'chandra_id'})
    return df[['ra_deg', 'dec_deg', 'chandra_id']]


def query_xmm(center: SkyCoord, radius: u.Quantity) -> pd.DataFrame:
    """Query the Flaccomio et al. (2023) XMM‑Newton source list via VizieR.

    The catalogue `J/A+A/670/A37/table2` provides 944 XMM‑Newton detections
    with positions in columns RAJ2000 and DEJ2000 and a source
    identification number in the `XMM` column【583771511037106†L107-L120】.  This function
    uses astroquery.vizier to perform a cone search and returns the
    results with uniform column names.  If the query yields no sources,
    an empty DataFrame with the expected columns is returned.
    """
    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'XMM'], catalog='J/A+A/670/A37/table2', row_limit=-1)
    tables = v.query_region(center, radius=radius)
    if len(tables) == 0:
        return pd.DataFrame(columns=['ra_deg', 'dec_deg', 'xmm_id'])
    tab = tables[0]
    df = tab.to_pandas()
    df = df.rename(columns={'RAJ2000': 'ra_deg', 'DEJ2000': 'dec_deg', 'XMM': 'xmm_id'})
    return df[['ra_deg', 'dec_deg', 'xmm_id']]


def cross_match(source_df: pd.DataFrame,
                other_df: pd.DataFrame,
                other_id_col: str,
                max_sep: u.Quantity) -> List[List[Optional[str]]]:
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
    src_coord = SkyCoord(source_df['ra_deg'].values * u.deg,
                          source_df['dec_deg'].values * u.deg)
    oth_coord = SkyCoord(other_df['ra_deg'].values * u.deg,
                          other_df['dec_deg'].values * u.deg)
    # match_to_catalog_sky finds the nearest neighbour for each source
    idx, sep2d, _ = src_coord.match_to_catalog_sky(oth_coord)
    matches = []
    for i, sep in enumerate(sep2d):
        if sep < max_sep:
            # There may be multiple potential matches within max_sep; find all
            close = oth_coord.separation(src_coord[i]) < max_sep
            ids = other_df.loc[close, other_id_col].astype(str).unique().tolist()
            matches.append(ids)
        else:
            matches.append([])
    return matches


def main() -> None:
    # Define the central coordinate of NGC 2264 and search radius
    center = SkyCoord(ra=100.25 * u.deg, dec=9.883333 * u.deg, frame='icrs')
    radius = 1.0 * u.deg

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

    # Prepare output DataFrame based on Gaia
    out_df = gaia.copy()
    # Add placeholder columns for matches
    out_df['tmass_id'] = None
    out_df['wise_id'] = None
    out_df['chandra_id'] = None
    out_df['xmm_id'] = None

    # Cross‑match with 2MASS (1 arcsec radius)
    print("Matching 2MASS with Gaia ...")
    tmass_matches = cross_match(out_df, tmass, '2MASS', 1.0 * u.arcsec)
    out_df['tmass_id'] = [",".join(ids) if ids else None for ids in tmass_matches]

    # Cross‑match with WISE (1 arcsec radius)
    print("Matching WISE with Gaia ...")
    wise_matches = cross_match(out_df, wise.rename(columns={'designation': 'wise_id'}), 'wise_id', 1.0 * u.arcsec)
    out_df['wise_id'] = [",".join(ids) if ids else None for ids in wise_matches]

    # Cross‑match with Chandra (5 arcsec radius)
    print("Matching Chandra with Gaia ...")
    # The column name in the Chandra table is now 'chandra_id' as returned by
    # query_chandra().  Ensure the identifier is correctly propagated into the
    # output table.
    chandra_matches = cross_match(out_df, chan, 'chandra_id', 5.0 * u.arcsec)
    out_df['chandra_id'] = [",".join(ids) if ids else None for ids in chandra_matches]

    # Cross‑match with XMM‑Newton (5 arcsec radius)
    print("Matching XMM‑Newton with Gaia ...")
    # The column name in the XMM table is now 'xmm_id' as returned by
    # query_xmm().
    xmm_matches = cross_match(out_df, xmm, 'xmm_id', 5.0 * u.arcsec)
    out_df['xmm_id'] = [",".join(ids) if ids else None for ids in xmm_matches]

    # Write to CSV
    out_df.to_csv('ngc2264_combined.csv', index=False)
    print("Combined catalogue written to ngc2264_combined.csv")


if __name__ == '__main__':
    main()