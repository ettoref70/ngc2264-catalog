# NGC 2264 Multi‑wavelength Catalog

This repository contains a Python script to assemble a multi‑band catalog of the star‑forming region **NGC 2264**.  The goal of the script is to download catalog data around NGC 2264 from several public archives, cross‑match the sources on the sky and produce a unified table of unique objects.

## Catalogs used

The following surveys/missions are queried:

1. **Gaia** (DR3) – provides high‑precision optical positions, proper motions and photometry for billions of stars.
2. **2MASS** (II/246) – near‑infrared photometry (J, H, K bands).
3. **WISE** (AllWISE PSC) – mid‑infrared photometry (W1–W4 bands).
4. **Chandra** – X‑ray detections from the *Chandra* ACIS observations of NGC 2264 compiled in the
   **Flaccomio et al. (2023)** study (table A.1 of A&A 670 A37).  The catalogue is accessed through
   VizieR (`J/A+A/670/A37/tablea1`) and provides RA/Dec and an ACIS identifier for each source.
5. **XMM‑Newton** – X‑ray detections from the *XMM‑Newton* observations of the same region, published
   in **Flaccomio et al. (2023)** (table 2 of A&A 670 A37).  This table is accessed via VizieR
   (`J/A+A/670/A37/table2`) and includes source positions and an XMM identifier.

## How it works

The script performs the following steps:

1. **Define the region of interest:** NGC 2264 is centred near RA ≈ 06h 41m (100.25°) and Dec ≈ +09° 53′【211272000973635†L2-L5】.  A search radius of 1 degree is used (corresponding to a 2×2 deg field).
2. **Query each catalogue:** Using the `astroquery` package, the script performs cone searches against each archive.  Gaia DR3 is queried via `Gaia.cone_search_async`; 2MASS via a VizieR cone search; AllWISE via the IRSA API; and the X‑ray lists are retrieved from VizieR tables `J/A+A/670/A37/tablea1` (Chandra) and `J/A+A/670/A37/table2` (XMM‑Newton) published by Flaccomio et al. (2023).
3. **Cross‑match sources:** The Gaia list is taken as the base catalogue because of its superior astrometric accuracy.  Each of the other tables is cross‑matched against the Gaia catalogue using `astropy.coordinates.SkyCoord` and a configurable matching radius (1 arcsec for infrared/optical, 5 arcsec for X‑ray).  When multiple matches occur within the tolerance, *all* candidate IDs are retained in the final table (stored as comma‑separated strings) instead of forcing a single match.
4. **Save results:** The script writes the unified catalogue to a CSV file (`ngc2264_combined.csv`) containing Gaia identifiers and positions along with the matched identifiers and photometry from the other surveys.

## Requirements

This script requires the following Python packages:

```bash
pip install astropy astroquery pandas
```

Network connectivity is required to query the online archives.  The script has been designed for clarity rather than speed; large radius searches may take several minutes.

## Usage

Run the script from the command line:

```bash
python cross_match_ngc2264.py
```

It will print progress messages to the console and write the combined catalogue to `ngc2264_combined.csv` in the current directory.

## Disclaimer

The cross‑matching radii and catalogues were chosen to provide a broad multi‑wavelength overview of NGC 2264.  Users may wish to adjust the search radius, matching thresholds or selected columns based on the specifics of their scientific application.