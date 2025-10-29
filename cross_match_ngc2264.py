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
import logging
import math
from typing import List, Optional
from collections import defaultdict
from contextlib import contextmanager, nullcontext
from time import perf_counter
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
import json
from typing import Dict, Any
from html import escape

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

# Configure logging for problem detection
logger = logging.getLogger(__name__)
# Check if we should enable debug logging via environment variable
if _os.environ.get('PROBLEM_DETECTION_DEBUG', '0') in ('1', 'true', 'yes', 'on'):
    logging.basicConfig(
        level=logging.DEBUG,
        format='[%(levelname)s] %(message)s'
    )
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)


import os
import base64
from pathlib import Path
# Directory to cache downloaded catalogs

CACHE_DIR = "catalogs"
os.makedirs(CACHE_DIR, exist_ok=True)

# Subdirectory for cached sky images (e.g., 2MASS J cutouts)
IMAGES_DIR = os.path.join(CACHE_DIR, "2MASS_J")
os.makedirs(IMAGES_DIR, exist_ok=True)

import argparse


# -------------------- Edit-link helpers --------------------
CURRENT_EDITS_RADIUS: float | None = None
NAV_KEYS_VERSION = "20251023d"

NAV_KEYS_TEMPLATE_PATH = Path(__file__).resolve().parent / 'aladin_scripts' / 'html' / 'nav_keys.js'
NAV_KEYS_MINIMAL = (
    "// nav_keys.js fallback: keyboard shortcuts unavailable\n"
    "(function(){ console.warn('nav_keys.js fallback loaded; keyboard shortcuts disabled.'); })();\n"
)
NAV_KEYS_EMBEDDED_B64 = (
"Ci8vIEtleWJvYXJkIG5hdmlnYXRpb246IOKGkC/ihpIgcGFnZXM7ICdhJyB0byB0cmlnZ2VyIEFs"
"YWRpbi4gQm91bmRzLWF3YXJlLgooZnVuY3Rpb24oKXsKICB2YXIgS0VZX0RFQlVHX1VSTCA9ICcv"
"YXBpL2tleV9kZWJ1Zyc7CiAgdmFyIF9teCA9IG51bGw7CiAgdmFyIF9teSA9IG51bGw7CgogIGRv"
"Y3VtZW50LmFkZEV2ZW50TGlzdGVuZXIoJ21vdXNlbW92ZScsIGZ1bmN0aW9uKGV2KXsKICAgIF9t"
"eCA9IGV2LmNsaWVudFg7CiAgICBfbXkgPSBldi5jbGllbnRZOwogIH0sIHtwYXNzaXZlOnRydWV9"
"KTsKCiAgZnVuY3Rpb24gaXNFZGl0YWJsZShlbCl7CiAgICBpZighZWwpIHJldHVybiBmYWxzZTsK"
"ICAgIHZhciB0YWcgPSAoZWwudGFnTmFtZSB8fCAnJykudG9VcHBlckNhc2UoKTsKICAgIHJldHVy"
"biBlbC5pc0NvbnRlbnRFZGl0YWJsZSB8fCB0YWcgPT09ICdJTlBVVCcgfHwgdGFnID09PSAnVEVY"
"VEFSRUEnIHx8IHRhZyA9PT0gJ1NFTEVDVCc7CiAgfQoKICBmdW5jdGlvbiBjbG9zZXN0KGVsLCBz"
"ZWwpewogICAgd2hpbGUoZWwgJiYgZWwubm9kZVR5cGUgPT09IDEpewogICAgICBpZihlbC5tYXRj"
"aGVzICYmIGVsLm1hdGNoZXMoc2VsKSkgcmV0dXJuIGVsOwogICAgICBlbCA9IGVsLnBhcmVudEVs"
"ZW1lbnQ7CiAgICB9CiAgICByZXR1cm4gbnVsbDsKICB9CgogIGZ1bmN0aW9uIF9hdXRvS2V5KCl7"
"CiAgICB0cnl7CiAgICAgIHZhciBiYXNlID0gKHR5cGVvZiBQQUdFX0tFWSA9PT0gJ3N0cmluZycg"
"JiYgUEFHRV9LRVkpID8gUEFHRV9LRVkgOiAnZ2xvYmFsJzsKICAgICAgcmV0dXJuICdBVVRPX0FM"
"QURJTl8nICsgYmFzZTsKICAgIH1jYXRjaChlKXsKICAgICAgcmV0dXJuICdBVVRPX0FMQURJTl9n"
"bG9iYWwnOwogICAgfQogIH0KCiAgZnVuY3Rpb24gX2NvbnN1bWVBdXRvU3VwcHJlc3MoKXsKICAg"
"IHZhciBzdXBwcmVzc2VkID0gZmFsc2U7CiAgICB0cnl7CiAgICAgIHN1cHByZXNzZWQgPSBzZXNz"
"aW9uU3RvcmFnZS5nZXRJdGVtKCdBVVRPX0FMQURJTl9TVVBQUkVTUycpID09PSAnMSc7CiAgICAg"
"IGlmKHN1cHByZXNzZWQpewogICAgICAgIHNlc3Npb25TdG9yYWdlLnJlbW92ZUl0ZW0oJ0FVVE9f"
"QUxBRElOX1NVUFBSRVNTJyk7CiAgICAgIH0KICAgIH1jYXRjaChlKXsKICAgICAgc3VwcHJlc3Nl"
"ZCA9IGZhbHNlOwogICAgfQogICAgcmV0dXJuIHN1cHByZXNzZWQ7CiAgfQoKICBmdW5jdGlvbiBy"
"ZWFkQXV0bygpewogICAgdHJ5ewogICAgICB2YXIga2V5ID0gX2F1dG9LZXkoKTsKICAgICAgdmFy"
"IHZhbCA9IG51bGw7CiAgICAgIGlmKHdpbmRvdy5sb2NhbFN0b3JhZ2UpewogICAgICAgIHZhbCA9"
"IGxvY2FsU3RvcmFnZS5nZXRJdGVtKGtleSk7CiAgICAgICAgaWYodmFsID09PSBudWxsKXsKICAg"
"ICAgICAgIHZhbCA9IGxvY2FsU3RvcmFnZS5nZXRJdGVtKCdBVVRPX0FMQURJTicpOwogICAgICAg"
"IH0KICAgICAgfQogICAgICByZXR1cm4gdmFsID09PSAnMSc7CiAgICB9Y2F0Y2goZSl7CiAgICAg"
"IHJldHVybiBmYWxzZTsKICAgIH0KICB9CgogIGZ1bmN0aW9uIHdyaXRlQXV0byhvbil7CiAgICB0"
"cnl7CiAgICAgIHZhciBrZXkgPSBfYXV0b0tleSgpOwogICAgICBpZih3aW5kb3cubG9jYWxTdG9y"
"YWdlKXsKICAgICAgICBsb2NhbFN0b3JhZ2Uuc2V0SXRlbShrZXksIG9uID8gJzEnIDogJzAnKTsK"
"ICAgICAgICBsb2NhbFN0b3JhZ2Uuc2V0SXRlbSgnQVVUT19BTEFESU4nLCBvbiA/ICcxJyA6ICcw"
"Jyk7CiAgICAgIH0KICAgIH1jYXRjaChlKXt9CiAgICBpZighb24pewogICAgICB0cnl7IHNlc3Np"
"b25TdG9yYWdlLnJlbW92ZUl0ZW0oJ0FVVE9fQUxBRElOX1NVUFBSRVNTJyk7IH1jYXRjaChfZSl7"
"fQogICAgfQogIH0KCiAgZnVuY3Rpb24gdHJpZ2dlckFsYWRpbigpewogICAgdmFyIGNhbmQgPSBu"
"dWxsOwogICAgdHJ5ewogICAgICBpZih0eXBlb2YgX214ID09PSAnbnVtYmVyJyAmJiB0eXBlb2Yg"
"X215ID09PSAnbnVtYmVyJyl7CiAgICAgICAgdmFyIGVsID0gZG9jdW1lbnQuZWxlbWVudEZyb21Q"
"b2ludChfbXgsIF9teSk7CiAgICAgICAgY2FuZCA9IGNsb3Nlc3QoZWwsICdhLmJ0bi5hbFtkYXRh"
"LWFqc10nKTsKICAgICAgfQogICAgfWNhdGNoKGUpe30KICAgIGlmKCFjYW5kKXsKICAgICAgdHJ5"
"eyBjYW5kID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcignYS5idG4uYWxbZGF0YS1hanNdJyk7IH1j"
"YXRjaChlKXt9CiAgICB9CiAgICBpZighY2FuZCl7CiAgICAgIHRyeXsgY2FuZCA9IGRvY3VtZW50"
"LnF1ZXJ5U2VsZWN0b3IoJ1tkYXRhLWF1dG8tYWxhZGluPSIxIl0nKTsgfWNhdGNoKGUpe30KICAg"
"IH0KICAgIGlmKGNhbmQpewogICAgICB0cnl7CiAgICAgICAgY2FuZC5jbGljaygpOwogICAgICAg"
"IHJldHVybiB0cnVlOwogICAgICB9Y2F0Y2goZSl7fQogICAgfQogICAgdHJ5ewogICAgICBpZihB"
"cnJheS5pc0FycmF5KHdpbmRvdy5BTEFESU5fSVRFTVMpICYmIHdpbmRvdy5BTEFESU5fSVRFTVMu"
"bGVuZ3RoKXsKICAgICAgICB2YXIgaWR4ID0gKHR5cGVvZiB3aW5kb3cuX3NyY0lkeCA9PT0gJ251"
"bWJlcicpID8gd2luZG93Ll9zcmNJZHggOiAwOwogICAgICAgIGlmKGlkeCA8IDAgfHwgaWR4ID49"
"IHdpbmRvdy5BTEFESU5fSVRFTVMubGVuZ3RoKXsKICAgICAgICAgIGlkeCA9IDA7CiAgICAgICAg"
"fQogICAgICAgIHZhciBhanMgPSB3aW5kb3cuQUxBRElOX0lURU1TW2lkeF0gfHwgd2luZG93LkFM"
"QURJTl9JVEVNU1swXTsKICAgICAgICB2YXIgYmFzZSA9ICh0eXBlb2Ygd2luZG93LlNBTVBfVVJM"
"ID09PSAnc3RyaW5nJyAmJiB3aW5kb3cuU0FNUF9VUkwpID8gd2luZG93LlNBTVBfVVJMIDogJ2h0"
"dHA6Ly8xMjcuMC4wLjE6ODc2NSc7CiAgICAgICAgaWYoYmFzZSAmJiBiYXNlLmNoYXJBdChiYXNl"
"Lmxlbmd0aCAtIDEpID09PSAnLycpewogICAgICAgICAgYmFzZSA9IGJhc2Uuc2xpY2UoMCwgLTEp"
"OwogICAgICAgIH0KICAgICAgICB2YXIgdXJsID0gYmFzZSArICcvcnVuX3NhbXA/ZmlsZT0nICsg"
"ZW5jb2RlVVJJQ29tcG9uZW50KGFqcykgKyAnJl90PScgKyBEYXRlLm5vdygpOwogICAgICAgIHZh"
"ciBzaW5rID0gZG9jdW1lbnQuZ2V0RWxlbWVudHNCeU5hbWUoJ2FsYWRpbl9zaW5rJylbMF07CiAg"
"ICAgICAgaWYoc2luayl7CiAgICAgICAgICB0cnl7IHNpbmsuc3JjID0gdXJsOyByZXR1cm4gdHJ1"
"ZTsgfWNhdGNoKGUpe30KICAgICAgICB9CiAgICAgICAgd2luZG93Lm9wZW4odXJsLCAnX2JsYW5r"
"Jyk7CiAgICAgICAgcmV0dXJuIHRydWU7CiAgICAgIH0KICAgIH1jYXRjaChlKXt9CiAgICByZXR1"
"cm4gZmFsc2U7CiAgfQoKICBmdW5jdGlvbiBmaW5kTmF2KGRpcil7CiAgICBkaXIgPSAoZGlyID09"
"PSAncHJldicpID8gJ3ByZXYnIDogJ25leHQnOwogICAgdHJ5ewogICAgICB2YXIgYW5jaG9ycyA9"
"IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3JBbGwoJy5uYXYgYVtocmVmXScpOwogICAgICBmb3IodmFy"
"IGk9MDtpPGFuY2hvcnMubGVuZ3RoO2krKyl7CiAgICAgICAgdmFyIGEgPSBhbmNob3JzW2ldOwog"
"ICAgICAgIGlmKCFhKSBjb250aW51ZTsKICAgICAgICB2YXIgbGFiZWwgPSAoYS50ZXh0Q29udGVu"
"dCB8fCAnJykudG9Mb3dlckNhc2UoKTsKICAgICAgICBpZihhLmRhdGFzZXQgJiYgYS5kYXRhc2V0"
"Lm5hdil7CiAgICAgICAgICB2YXIgZHMgPSBTdHJpbmcoYS5kYXRhc2V0Lm5hdiB8fCAnJykudG9M"
"b3dlckNhc2UoKTsKICAgICAgICAgIGlmKGRzID09PSBkaXIpIHJldHVybiBhOwogICAgICAgIH0K"
"ICAgICAgICBpZihkaXIgPT09ICdwcmV2Jyl7CiAgICAgICAgICBpZihsYWJlbC5pbmRleE9mKCdw"
"cmV2IHBhZ2UnKSAhPT0gLTEgfHwgbGFiZWwuaW5kZXhPZigncHJldmlvdXMgcGFnZScpICE9PSAt"
"MSkgcmV0dXJuIGE7CiAgICAgICAgICBpZihsYWJlbC5pbmRleE9mKCfCqyBwcmV2JykgIT09IC0x"
"IHx8IGxhYmVsLmluZGV4T2YoJ3ByZXYgwqsnKSAhPT0gLTEpIHJldHVybiBhOwogICAgICAgIH1l"
"bHNlIGlmKGRpciA9PT0gJ25leHQnKXsKICAgICAgICAgIGlmKGxhYmVsLmluZGV4T2YoJ25leHQg"
"cGFnZScpICE9PSAtMSkgcmV0dXJuIGE7CiAgICAgICAgICBpZihsYWJlbC5pbmRleE9mKCfCuyBu"
"ZXh0JykgIT09IC0xIHx8IGxhYmVsLmluZGV4T2YoJ25leHQgwrsnKSAhPT0gLTEpIHJldHVybiBh"
"OwogICAgICAgIH0KICAgICAgfQogICAgfWNhdGNoKGUpe30KICAgIHRyeXsKICAgICAgdmFyIGJ0"
"biA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoJy5uYXYgW2RhdGEtbmF2PSInICsgZGlyICsgJyJd"
"Jyk7CiAgICAgIGlmKGJ0bil7CiAgICAgICAgdmFyIGhyZWYgPSBidG4uZ2V0QXR0cmlidXRlKCdo"
"cmVmJykgfHwgYnRuLmdldEF0dHJpYnV0ZSgnZGF0YS1ocmVmJyk7CiAgICAgICAgaWYoaHJlZil7"
"CiAgICAgICAgICB2YXIgYXV0byA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnQoJ2EnKTsKICAgICAg"
"ICAgIGF1dG8uc2V0QXR0cmlidXRlKCdocmVmJywgaHJlZik7CiAgICAgICAgICByZXR1cm4gYXV0"
"bzsKICAgICAgICB9CiAgICAgIH0KICAgIH1jYXRjaChlKXt9CiAgICB0cnl7CiAgICAgIGlmKHR5"
"cGVvZiBQQUdFX0tFWSA9PT0gJ3N0cmluZycgJiYgUEFHRV9LRVkpewogICAgICAgIHZhciBjdXJy"
"ZW50ID0gTnVtYmVyKFBBR0VfTlVNIHx8IDApOwogICAgICAgIGlmKCFOdW1iZXIuaXNGaW5pdGUo"
"Y3VycmVudCkgfHwgY3VycmVudCA8IDEpIGN1cnJlbnQgPSAxOwogICAgICAgIHZhciB0b3RhbCA9"
"IF90b3RhbFBhZ2VzKCk7CiAgICAgICAgdmFyIHRhcmdldCA9IGRpciA9PT0gJ3ByZXYnID8gKGN1"
"cnJlbnQgLSAxKSA6IChjdXJyZW50ICsgMSk7CiAgICAgICAgaWYoTnVtYmVyLmlzRmluaXRlKHRv"
"dGFsKSAmJiB0b3RhbCA+IDApewogICAgICAgICAgaWYodGFyZ2V0IDwgMSB8fCB0YXJnZXQgPiB0"
"b3RhbCkgcmV0dXJuIG51bGw7CiAgICAgICAgfWVsc2UgaWYodGFyZ2V0IDwgMSl7CiAgICAgICAg"
"ICByZXR1cm4gbnVsbDsKICAgICAgICB9CiAgICAgICAgaWYodGFyZ2V0ID09PSBjdXJyZW50KSBy"
"ZXR1cm4gbnVsbDsKICAgICAgICB2YXIgYSA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnQoJ2EnKTsK"
"ICAgICAgICB2YXIgaGFzaCA9ICcjc2Vlaz0nICsgZGlyOwogICAgICAgIGEuc2V0QXR0cmlidXRl"
"KCdocmVmJywgUEFHRV9LRVkgKyAnX3BhZ2UnICsgdGFyZ2V0ICsgJy5odG1sJyArIGhhc2gpOwog"
"ICAgICAgIHJldHVybiBhOwogICAgICB9CiAgICB9Y2F0Y2goZSl7fQogICAgdHJ5ewogICAgICB2"
"YXIgcGF0aCA9IFN0cmluZyhsb2NhdGlvbi5wYXRobmFtZSB8fCAnJyk7CiAgICAgIHZhciBtID0g"
"cGF0aC5tYXRjaCgvKFteXC9dKylfcGFnZShcZCspXC5odG1sJC8pOwogICAgICBpZihtKXsKICAg"
"ICAgICB2YXIgYmFzZSA9IG1bMV07CiAgICAgICAgdmFyIG51bSA9IHBhcnNlSW50KG1bMl0sIDEw"
"KSB8fCAwOwogICAgICAgIHZhciB0YXJnZXQyID0gZGlyID09PSAncHJldicgPyAobnVtIC0gMSkg"
"OiAobnVtICsgMSk7CiAgICAgICAgaWYodGFyZ2V0MiA+PSAxKXsKICAgICAgICAgIHZhciBzdHVi"
"ID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudCgnYScpOwogICAgICAgICAgc3R1Yi5zZXRBdHRyaWJ1"
"dGUoJ2hyZWYnLCBiYXNlICsgJ19wYWdlJyArIHRhcmdldDIgKyAnLmh0bWwjc2Vlaz0nICsgZGly"
"KTsKICAgICAgICAgIHJldHVybiBzdHViOwogICAgICAgIH0KICAgICAgfQogICAgfWNhdGNoKGUp"
"e30KICAgIHJldHVybiBudWxsOwogIH0KCiAgZnVuY3Rpb24gX3RvdGFsUGFnZXMoKXsKICAgIHRy"
"eXsKICAgICAgcmV0dXJuICh0eXBlb2YgVE9UQUxfUEFHRVMgPT09ICdudW1iZXInICYmIFRPVEFM"
"X1BBR0VTID4gMCkgPyBUT1RBTF9QQUdFUyA6IDE7CiAgICB9Y2F0Y2goZSl7CiAgICAgIHJldHVy"
"biAxOwogICAgfQogIH0KCiAgZnVuY3Rpb24gbG9nS2V5KGV2dCl7CiAgICB0cnl7CiAgICAgIHZh"
"ciBwYXJhbXMgPSBbCiAgICAgICAgJ2tleT0nICsgZW5jb2RlVVJJQ29tcG9uZW50KGV2dC5rZXkg"
"fHwgJycpLAogICAgICAgICdjb2RlPScgKyBlbmNvZGVVUklDb21wb25lbnQoZXZ0LmNvZGUgfHwg"
"JycpLAogICAgICAgICdjdHJsPScgKyAoZXZ0LmN0cmxLZXkgPyAnMScgOiAnMCcpLAogICAgICAg"
"ICdhbHQ9JyArIChldnQuYWx0S2V5ID8gJzEnIDogJzAnKSwKICAgICAgICAnc2hpZnQ9JyArIChl"
"dnQuc2hpZnRLZXkgPyAnMScgOiAnMCcpLAogICAgICAgICdtZXRhPScgKyAoZXZ0Lm1ldGFLZXkg"
"PyAnMScgOiAnMCcpLAogICAgICAgICd0cz0nICsgRGF0ZS5ub3coKQogICAgICBdLmpvaW4oJyYn"
"KTsKICAgICAgdmFyIGltZyA9IG5ldyBJbWFnZSgpOwogICAgICBpbWcuc3JjID0gS0VZX0RFQlVH"
"X1VSTCArICc/JyArIHBhcmFtczsKICAgIH1jYXRjaChlKXt9CiAgfQoKICBmdW5jdGlvbiBfaXNQ"
"cm9ibGVtUGFnZShsaXN0LCBwYWdlKXsKICAgIGlmKCFBcnJheS5pc0FycmF5KGxpc3QpKSByZXR1"
"cm4gZmFsc2U7CiAgICB2YXIgdGFyZ2V0ID0gTnVtYmVyKHBhZ2UpOwogICAgaWYoIU51bWJlci5p"
"c0Zpbml0ZSh0YXJnZXQpKSByZXR1cm4gZmFsc2U7CiAgICBmb3IodmFyIGk9MDtpPGxpc3QubGVu"
"Z3RoO2krKyl7CiAgICAgIGlmKE51bWJlcihsaXN0W2ldKSA9PT0gdGFyZ2V0KSByZXR1cm4gdHJ1"
"ZTsKICAgIH0KICAgIHJldHVybiBmYWxzZTsKICB9CgogIGZ1bmN0aW9uIF9oYXNVbnNraXBwZWRQ"
"cm9ibGVtcygpewogICAgdHJ5ewogICAgICB2YXIgcGIgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9y"
"QWxsKCcucGFuZWxib3gnKTsKICAgICAgZm9yKHZhciBpPTA7aTxwYi5sZW5ndGg7aSsrKXsKICAg"
"ICAgICB2YXIgcCA9IHBiW2ldOwogICAgICAgIGlmKHAuZ2V0QXR0cmlidXRlKCdkYXRhLXByb2Js"
"ZW0nKSA9PT0gJzEnICYmIHAuZ2V0QXR0cmlidXRlKCdkYXRhLXNraXAnKSAhPT0gJzEnKXsKICAg"
"ICAgICAgIHJldHVybiB0cnVlOwogICAgICAgIH0KICAgICAgfQogICAgfWNhdGNoKGUpe30KICAg"
"IHJldHVybiBmYWxzZTsKICB9CgogIGZ1bmN0aW9uIF9tYXliZUF1dG9BZHZhbmNlKCl7CiAgICB0"
"cnl7CiAgICAgIHZhciBvbmx5T24gPSAoZnVuY3Rpb24oKXsKICAgICAgICB0cnl7CiAgICAgICAg"
"ICB2YXIga2V5ID0gJ09OTFlfUFJPQl8nICsgKFBBR0VfS0VZIHx8ICcnKTsKICAgICAgICAgIHJl"
"dHVybiBsb2NhbFN0b3JhZ2UuZ2V0SXRlbShrZXkpID09PSAnMSc7CiAgICAgICAgfWNhdGNoKGUp"
"eyByZXR1cm4gZmFsc2U7IH0KICAgICAgfSkoKTsKICAgICAgaWYoIW9ubHlPbikgcmV0dXJuIGZh"
"bHNlOwogICAgICBpZihfaGFzVW5za2lwcGVkUHJvYmxlbXMoKSkgcmV0dXJuIGZhbHNlOwogICAg"
"ICB2YXIgZGlyID0gKGZ1bmN0aW9uKCl7CiAgICAgICAgdHJ5ewogICAgICAgICAgdmFyIGggPSBT"
"dHJpbmcobG9jYXRpb24uaGFzaCB8fCAnJyk7CiAgICAgICAgICB2YXIgbSA9IGgubWF0Y2goL3Nl"
"ZWs9KG5leHR8cHJldikvKTsKICAgICAgICAgIHJldHVybiBtID8gbVsxXSA6ICduZXh0JzsKICAg"
"ICAgICB9Y2F0Y2goZSl7IHJldHVybiAnbmV4dCc7IH0KICAgICAgfSkoKTsKICAgICAgdmFyIHRv"
"dGFsID0gX3RvdGFsUGFnZXMoKTsKICAgICAgdmFyIGN1cnJlbnQgPSBOdW1iZXIoUEFHRV9OVU0g"
"fHwgMCk7CiAgICAgIGlmKCFOdW1iZXIuaXNGaW5pdGUoY3VycmVudCkgfHwgY3VycmVudCA8IDEp"
"eyBjdXJyZW50ID0gMTsgfQogICAgICB2YXIgdGFyZ2V0ID0gKGRpciA9PT0gJ3ByZXYnKSA/IChj"
"dXJyZW50IC0gMSkgOiAoY3VycmVudCArIDEpOwogICAgICBpZih0YXJnZXQgPCAxIHx8IHRhcmdl"
"dCA+IHRvdGFsKXsKICAgICAgICByZXR1cm4gZmFsc2U7CiAgICAgIH0KICAgICAgaWYodGFyZ2V0"
"ICE9PSBjdXJyZW50KXsKICAgICAgICB3aW5kb3cubG9jYXRpb24ucmVwbGFjZShQQUdFX0tFWSAr"
"ICdfcGFnZScgKyB0YXJnZXQgKyAnLmh0bWwjc2Vlaz0nICsgZGlyKTsKICAgICAgICByZXR1cm4g"
"dHJ1ZTsKICAgICAgfQogICAgfWNhdGNoKGUpe30KICAgIHJldHVybiBmYWxzZTsKICB9CgogIGZ1"
"bmN0aW9uIG1heWJlVHJpZ2dlckF1dG8oKXsKICAgIHRyeXsKICAgICAgdmFyIGNiID0gZG9jdW1l"
"bnQuZ2V0RWxlbWVudEJ5SWQoJ2F1dG9fYWxhZGluJyk7CiAgICAgIGlmKGNiICYmIGNiLmNoZWNr"
"ZWQpewogICAgICAgIHNldFRpbWVvdXQoZnVuY3Rpb24oKXsgdHJpZ2dlckFsYWRpbigpOyB9LCAy"
"MDApOwogICAgICB9CiAgICB9Y2F0Y2goZSl7fQogIH0KCiAgZnVuY3Rpb24gY29sbGVjdFRvZ2ds"
"ZXMoKXsKICAgIHZhciBidXR0b25zID0gQXJyYXkucHJvdG90eXBlLnNsaWNlLmNhbGwoZG9jdW1l"
"bnQucXVlcnlTZWxlY3RvckFsbCgnLmJ0bi50Z2wnKSk7CiAgICBpZighYnV0dG9ucy5sZW5ndGgp"
"IHJldHVybiBidXR0b25zOwogICAgdmFyIHBhbmVsID0gbnVsbDsKICAgIGlmKHR5cGVvZiBfbXgg"
"PT09ICdudW1iZXInICYmIHR5cGVvZiBfbXkgPT09ICdudW1iZXInKXsKICAgICAgdmFyIGVsID0g"
"ZG9jdW1lbnQuZWxlbWVudEZyb21Qb2ludChfbXgsIF9teSk7CiAgICAgIHBhbmVsID0gY2xvc2Vz"
"dChlbCwgJy5wYW5lbGJveCcpOwogICAgICBpZighcGFuZWwpewogICAgICAgIHZhciBidG4gPSBj"
"bG9zZXN0KGVsLCAnLmJ0bi50Z2wnKTsKICAgICAgICBpZihidG4peyBwYW5lbCA9IGNsb3Nlc3Qo"
"YnRuLCAnLnBhbmVsYm94Jyk7IH0KICAgICAgfQogICAgfQogICAgaWYocGFuZWwpewogICAgICB2"
"YXIgYmFzZSA9IHBhbmVsLmdldEF0dHJpYnV0ZSgnZGF0YS1iYXNlJyk7CiAgICAgIHZhciBmaWx0"
"ZXJlZCA9IGJ1dHRvbnMuZmlsdGVyKGZ1bmN0aW9uKGJ0bil7CiAgICAgICAgaWYoYmFzZSl7CiAg"
"ICAgICAgICB2YXIgaWQgPSBidG4uZ2V0QXR0cmlidXRlKCdkYXRhLXBhbmVsJyk7CiAgICAgICAg"
"ICBpZihpZCkgcmV0dXJuIGlkID09PSBiYXNlOwogICAgICAgIH0KICAgICAgICByZXR1cm4gcGFu"
"ZWwuY29udGFpbnMoYnRuKTsKICAgICAgfSk7CiAgICAgIGlmKGZpbHRlcmVkLmxlbmd0aCkgYnV0"
"dG9ucyA9IGZpbHRlcmVkOwogICAgfQogICAgYnV0dG9ucy5zb3J0KGZ1bmN0aW9uKGEsIGIpewog"
"ICAgICB2YXIgcmEgPSBhLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpOwogICAgICB2YXIgcmIgPSBi"
"LmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpOwogICAgICBpZihNYXRoLmFicyhyYS50b3AgLSByYi50"
"b3ApID4gNCkgcmV0dXJuIHJhLnRvcCAtIHJiLnRvcDsKICAgICAgcmV0dXJuIHJhLmxlZnQgLSBy"
"Yi5sZWZ0OwogICAgfSk7CiAgICByZXR1cm4gYnV0dG9uczsKICB9CgogIGZ1bmN0aW9uIHRyaWdn"
"ZXJUb2dnbGUocmFuayl7CiAgICB2YXIgYnV0dG9ucyA9IGNvbGxlY3RUb2dnbGVzKCk7CiAgICBp"
"ZighYnV0dG9ucy5sZW5ndGgpIHJldHVybiBmYWxzZTsKICAgIHZhciBpZHggPSBNYXRoLm1pbihN"
"YXRoLm1heChyYW5rIC0gMSwgMCksIGJ1dHRvbnMubGVuZ3RoIC0gMSk7CiAgICB2YXIgYnRuID0g"
"YnV0dG9uc1tpZHhdOwogICAgaWYoIWJ0bikgcmV0dXJuIGZhbHNlOwogICAgdmFyIGhyZWYgPSBi"
"dG4uZ2V0QXR0cmlidXRlKCdkYXRhLWhyZWYtZm9yY2UnKSB8fCBidG4uZ2V0QXR0cmlidXRlKCdo"
"cmVmJyk7CiAgICBpZighaHJlZil7CiAgICAgIHRyeXsgYnRuLmNsaWNrKCk7IHJldHVybiB0cnVl"
"OyB9Y2F0Y2goZSl7IHJldHVybiBmYWxzZTsgfQogICAgfQogICAgd2luZG93LmxvY2F0aW9uLmhy"
"ZWYgPSBocmVmOwogICAgcmV0dXJuIHRydWU7CiAgfQoKICBmdW5jdGlvbiBzZXR1cEF1dG8oc2tp"
"cEltbWVkaWF0ZSl7CiAgICB0cnl7CiAgICAgIHZhciBjYiA9IGRvY3VtZW50LmdldEVsZW1lbnRC"
"eUlkKCdhdXRvX2FsYWRpbicpOwogICAgICBpZighY2IpIHJldHVybjsKICAgICAgdmFyIG9uID0g"
"cmVhZEF1dG8oKTsKICAgICAgdmFyIHN1cHByZXNzZWQgPSBfY29uc3VtZUF1dG9TdXBwcmVzcygp"
"OwogICAgICBjYi5jaGVja2VkID0gb247CiAgICAgIGNiLmFkZEV2ZW50TGlzdGVuZXIoJ2NoYW5n"
"ZScsIGZ1bmN0aW9uKCl7CiAgICAgICAgd3JpdGVBdXRvKGNiLmNoZWNrZWQpOwogICAgICAgIGlm"
"KGNiLmNoZWNrZWQpewogICAgICAgICAgdHJ5eyBzZXNzaW9uU3RvcmFnZS5yZW1vdmVJdGVtKCdB"
"VVRPX0FMQURJTl9TVVBQUkVTUycpOyB9Y2F0Y2goX2Upe30KICAgICAgICAgIHZhciBob3BwZWQg"
"PSBfbWF5YmVBdXRvQWR2YW5jZSgpOwogICAgICAgICAgaWYoIWhvcHBlZCl7CiAgICAgICAgICAg"
"IG1heWJlVHJpZ2dlckF1dG8oKTsKICAgICAgICAgIH0KICAgICAgICB9CiAgICAgIH0pOwogICAg"
"ICBpZihvbiAmJiAhc3VwcHJlc3NlZCAmJiAhc2tpcEltbWVkaWF0ZSl7CiAgICAgICAgbWF5YmVU"
"cmlnZ2VyQXV0bygpOwogICAgICB9CiAgICB9Y2F0Y2goZSl7fQogIH0KCiAgZnVuY3Rpb24gaXNQ"
"cm9ibGVtUGFuZWwoZWwpewogICAgaWYoIWVsKSByZXR1cm4gZmFsc2U7CiAgICB0cnl7IHJldHVy"
"biBlbC5nZXRBdHRyaWJ1dGUoJ2RhdGEtcHJvYmxlbScpID09PSAnMSc7IH1jYXRjaChlKXsgcmV0"
"dXJuIGZhbHNlOyB9CiAgfQoKICBpZighQXJyYXkuaXNBcnJheSh3aW5kb3cuUFJPQkxFTV9QQUdF"
"UykpewogICAgd2luZG93LlBST0JMRU1fUEFHRVMgPSBbXTsKICB9CgogIGRvY3VtZW50LmFkZEV2"
"ZW50TGlzdGVuZXIoJ2tleWRvd24nLCBmdW5jdGlvbihlKXsKICAgIGlmKGlzRWRpdGFibGUoZS50"
"YXJnZXQpKSByZXR1cm47CiAgICBsb2dLZXkoZSk7CiAgICB2YXIgayA9IChlLmtleSB8fCAnJyk7"
"CiAgICBpZighZS5jdHJsS2V5ICYmICFlLm1ldGFLZXkgJiYgIWUuYWx0S2V5KXsKICAgICAgaWYo"
"ayA9PT0gJzEnIHx8IGsgPT09ICcyJyB8fCBrID09PSAnMycpewogICAgICAgIGlmKHRyaWdnZXJU"
"b2dnbGUocGFyc2VJbnQoaywgMTApKSl7CiAgICAgICAgICBlLnByZXZlbnREZWZhdWx0KCk7CiAg"
"ICAgICAgICByZXR1cm47CiAgICAgICAgfQogICAgICB9CiAgICB9CiAgICB2YXIgb25seSA9IGRv"
"Y3VtZW50LmdldEVsZW1lbnRCeUlkKCdvbmx5X3Byb2InKTsKICAgIHZhciBvbmx5T24gPSAhIShv"
"bmx5ICYmIG9ubHkuY2hlY2tlZCk7CiAgICBpZihrID09PSAnQXJyb3dMZWZ0Jyl7CiAgICAgIGlm"
"KG9ubHlPbil7CiAgICAgICAgZS5wcmV2ZW50RGVmYXVsdCgpOwogICAgICAgIHZhciBjdXJyZW50"
"TCA9IE51bWJlcihQQUdFX05VTSB8fCAwKTsKICAgICAgICBpZighTnVtYmVyLmlzRmluaXRlKGN1"
"cnJlbnRMKSB8fCBjdXJyZW50TCA8IDEpeyBjdXJyZW50TCA9IDE7IH0KICAgICAgICBpZihjdXJy"
"ZW50TCA+IDEpewogICAgICAgICAgd2luZG93LmxvY2F0aW9uLmhyZWYgPSBQQUdFX0tFWSArICdf"
"cGFnZScgKyAoY3VycmVudEwgLSAxKSArICcuaHRtbCNzZWVrPXByZXYnOwogICAgICAgIH0KICAg"
"ICAgICByZXR1cm47CiAgICAgIH0KICAgICAgdmFyIHByZXYgPSBmaW5kTmF2KCdwcmV2Jyk7CiAg"
"ICAgIGlmKHByZXYgJiYgcHJldi5ocmVmKXsKICAgICAgICBlLnByZXZlbnREZWZhdWx0KCk7CiAg"
"ICAgICAgd2luZG93LmxvY2F0aW9uLmhyZWYgPSBwcmV2LmhyZWY7CiAgICAgIH0KICAgIH0KICAg"
"IGVsc2UgaWYoayA9PT0gJ0Fycm93UmlnaHQnKXsKICAgICAgaWYob25seU9uKXsKICAgICAgICBl"
"LnByZXZlbnREZWZhdWx0KCk7CiAgICAgICAgdmFyIGN1cnJlbnRSID0gTnVtYmVyKFBBR0VfTlVN"
"IHx8IDApOwogICAgICAgIGlmKCFOdW1iZXIuaXNGaW5pdGUoY3VycmVudFIpIHx8IGN1cnJlbnRS"
"IDwgMSl7IGN1cnJlbnRSID0gMTsgfQogICAgICAgIHZhciB0b3RhbFBhZ2VzID0gX3RvdGFsUGFn"
"ZXMoKTsKICAgICAgICBpZihjdXJyZW50UiA8IHRvdGFsUGFnZXMpewogICAgICAgICAgd2luZG93"
"LmxvY2F0aW9uLmhyZWYgPSBQQUdFX0tFWSArICdfcGFnZScgKyAoY3VycmVudFIgKyAxKSArICcu"
"aHRtbCNzZWVrPW5leHQnOwogICAgICAgIH0KICAgICAgICByZXR1cm47CiAgICAgIH0KICAgICAg"
"dmFyIG5leHQgPSBmaW5kTmF2KCduZXh0Jyk7CiAgICAgIGlmKG5leHQgJiYgbmV4dC5ocmVmKXsK"
"ICAgICAgICBlLnByZXZlbnREZWZhdWx0KCk7CiAgICAgICAgd2luZG93LmxvY2F0aW9uLmhyZWYg"
"PSBuZXh0LmhyZWY7CiAgICAgIH0KICAgIH0KICAgIGVsc2UgaWYoay50b0xvd2VyQ2FzZSgpID09"
"PSAnYScgJiYgIWUuY3RybEtleSAmJiAhZS5tZXRhS2V5ICYmICFlLmFsdEtleSl7CiAgICAgIGlm"
"KHRyaWdnZXJBbGFkaW4oKSl7IGUucHJldmVudERlZmF1bHQoKTsgfQogICAgfQogICAgZWxzZSBp"
"ZihrLnRvTG93ZXJDYXNlKCkgPT09ICdnJyAmJiAhZS5jdHJsS2V5ICYmICFlLm1ldGFLZXkgJiYg"
"IWUuYWx0S2V5KXsKICAgICAgdHJ5ewogICAgICAgIHZhciBpbnAgPSBkb2N1bWVudC5xdWVyeVNl"
"bGVjdG9yKCcubmF2IGZvcm0uanVtcCBpbnB1dFt0eXBlPW51bWJlcl0nKTsKICAgICAgICBpZihp"
"bnApeyBpbnAuZm9jdXMoKTsgaW5wLnNlbGVjdCgpOyBlLnByZXZlbnREZWZhdWx0KCk7IHJldHVy"
"bjsgfQogICAgICB9Y2F0Y2goX2Upe30KICAgICAgdHJ5ewogICAgICAgIHZhciB0b3RhbCA9IF90"
"b3RhbFBhZ2VzKCk7CiAgICAgICAgdmFyIGRlZnYgPSBTdHJpbmcoKHR5cGVvZiBQQUdFX05VTSA9"
"PT0gJ251bWJlcicpID8gUEFHRV9OVU0gOiAxKTsKICAgICAgICB2YXIgdiA9IHByb21wdCgnR28g"
"dG8gcGFnZSAoMS4uJyArIHRvdGFsICsgJyk6JywgZGVmdik7CiAgICAgICAgaWYodiAhPT0gbnVs"
"bCl7IGUucHJldmVudERlZmF1bHQoKTsgZ290b1BhZ2Uodik7IHJldHVybjsgfQogICAgICB9Y2F0"
"Y2goX2Upe30KICAgIH0KICB9LCB0cnVlKTsKCiAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcign"
"RE9NQ29udGVudExvYWRlZCcsIGZ1bmN0aW9uKCl7CiAgICB0cnl7CiAgICAgIGNvbnN0IHNrPV9y"
"ZWFkU2tpcCgpOwogICAgICBkb2N1bWVudC5xdWVyeVNlbGVjdG9yQWxsKCcucGFuZWxib3gnKS5m"
"b3JFYWNoKGZ1bmN0aW9uKHApewogICAgICAgIGNvbnN0IGI9cC5nZXRBdHRyaWJ1dGUoJ2RhdGEt"
"YmFzZScpOwogICAgICAgIGlmKHNrW2JdKXsKICAgICAgICAgIHAuc2V0QXR0cmlidXRlKCdkYXRh"
"LXNraXAnLCcxJyk7CiAgICAgICAgICBjb25zdCBjYj1wLnF1ZXJ5U2VsZWN0b3IoJy5za2lwYm94"
"Jyk7CiAgICAgICAgICBpZihjYikgY2IuY2hlY2tlZD10cnVlOwogICAgICAgIH0KICAgICAgfSk7"
"CiAgICB9Y2F0Y2goZSl7fQogICAgY29uc3Qgb25seT1kb2N1bWVudC5nZXRFbGVtZW50QnlJZCgn"
"b25seV9wcm9iJyk7CiAgICBpZihvbmx5KXsKICAgICAgb25seS5jaGVja2VkID0gX3JlYWRPbmx5"
"KCk7CiAgICAgIG9ubHkuYWRkRXZlbnRMaXN0ZW5lcignY2hhbmdlJywgZnVuY3Rpb24oKXsKICAg"
"ICAgICBfd3JpdGVPbmx5KG9ubHkuY2hlY2tlZCk7CiAgICAgICAgdXBkYXRlUGFuZWxzKCk7CiAg"
"ICAgICAgaWYob25seS5jaGVja2VkKXsKICAgICAgICAgIHZhciBob3BwZWQgPSBfbWF5YmVBdXRv"
"QWR2YW5jZSgpOwogICAgICAgICAgaWYoIWhvcHBlZCl7CiAgICAgICAgICAgIG1heWJlVHJpZ2dl"
"ckF1dG8oKTsKICAgICAgICAgIH0KICAgICAgICB9CiAgICAgIH0pOwogICAgfQogICAgdXBkYXRl"
"UGFuZWxzKCk7CiAgICB2YXIgYXV0b0hvcHBlZCA9IF9tYXliZUF1dG9BZHZhbmNlKCk7CiAgICBz"
"ZXR1cEF1dG8oYXV0b0hvcHBlZCk7CiAgICBpZighYXV0b0hvcHBlZCl7CiAgICAgIG1heWJlVHJp"
"Z2dlckF1dG8oKTsKICAgIH0KICB9LCB0cnVlKTsKfSkoKTsK"
)

def _embedded_nav_keys() -> str:
    try:
        return base64.b64decode(NAV_KEYS_EMBEDDED_B64).decode('utf-8')
    except Exception:
        return NAV_KEYS_MINIMAL


def _load_nav_keys_template() -> str:
    try:
        if NAV_KEYS_TEMPLATE_PATH.exists():
            return NAV_KEYS_TEMPLATE_PATH.read_text()
    except Exception:
        pass
    return _embedded_nav_keys()


NAV_JS_TEMPLATE = _load_nav_keys_template()

# -------------------- Blend detection constants --------------------
# Thresholds for identifying 2MASS blends (multiple Gaia sources in unresolved 2MASS source)
BLEND_GROUP_RADIUS_ARCSEC = 1.5  # Max separation from 2MASS position for group members
BLEND_MAX_PAIR_SEP_ARCSEC = 2.0  # Max pairwise distance between Gaia members
BLEND_CENTROID_TOLERANCE_ARCSEC = 0.5  # Max distance between 2MASS and Gaia centroid
BLEND_MAX_MAGNITUDE_DIFF = 1.0  # Max G-band magnitude range in group


def _effective_radius(radius_deg: float | None = None) -> float | None:
    if radius_deg is not None:
        return radius_deg
    try:
        if CURRENT_EDITS_RADIUS is not None:
            return float(CURRENT_EDITS_RADIUS)
    except Exception:
        pass
    for var in ('DYN_RADIUS_DEG', 'RADIUS_DEG', 'CURRENT_EDITS_RADIUS'):
        val = os.environ.get(var)
        if val is None:
            continue
        try:
            return float(val)
        except Exception:
            continue
    return None


class _PlotProfiler:
    """Lightweight optional profiler for plot_after_merge."""

    def __init__(self, enabled: bool, prefix: str = "plot_after_merge"):
        self.enabled = bool(enabled)
        if not self.enabled:
            return
        self.prefix = prefix
        self._totals: dict[str, float] = defaultdict(float)
        self._counts: dict[str, int] = defaultdict(int)
        self._stack: list[tuple[str, float]] = []

    def tic(self, name: str) -> None:
        if not self.enabled:
            return
        self._stack.append((name, perf_counter()))

    def toc(self, name: str) -> float:
        if not self.enabled:
            return 0.0
        for idx in range(len(self._stack) - 1, -1, -1):
            tag, start = self._stack[idx]
            if tag == name:
                self._stack.pop(idx)
                dt = perf_counter() - start
                self._totals[name] += dt
                self._counts[name] += 1
                return dt
        return 0.0

    @contextmanager
    def section(self, name: str):
        self.tic(name)
        try:
            yield
        finally:
            self.toc(name)

    def add(self, name: str, duration: float, *, count: int = 1) -> None:
        if not self.enabled:
            return
        self._totals[name] += float(duration)
        self._counts[name] += int(count)

    def emit(self, label: str = "") -> None:
        if not self.enabled:
            return
        header = f"[profile] {self.prefix}" if not label else f"[profile] {self.prefix} {label}"
        print(header)
        for key, total in sorted(self._totals.items(), key=lambda kv: kv[1], reverse=True):
            calls = self._counts.get(key, 0)
            avg = (total / calls) if calls else total
            print(f"[profile]  {key:24s} total={total*1000:.1f} ms  calls={calls:3d}  avg={avg*1000:.1f} ms")


@contextmanager
def _diag_timer(profiler: _PlotProfiler, log: list[tuple[str, float]] | None, label: str):
    start = perf_counter()
    try:
        yield
    finally:
        duration = perf_counter() - start
        profiler.add(f'diag_{label}', duration)
        if log is not None:
            log.append((label, duration))

def _radius_suffix(radius_deg: float | None) -> str:
    radius = _effective_radius(radius_deg)
    if radius is None:
        return ''
    radius_str = f"{radius:.3f}".rstrip('0').rstrip('.')
    if not radius_str:
        radius_str = '0'
    return f"_r{radius_str}"


def _edits_path(radius_deg: float | None = None) -> str:
    suffix = _radius_suffix(radius_deg)
    base = os.environ.get('EDIT_LINKS_PATH')
    if base:
        base_path = Path(base)
        if base_path.suffix:
            if suffix:
                return str(base_path.with_name(base_path.stem + suffix + base_path.suffix))
            return str(base_path)
        # Treat as directory (even if it does not yet exist)
        return str((base_path / f"links{suffix or ''}.json"))
    base_dir = Path('edits')
    try:
        base_dir.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    filename = f"links{suffix or ''}.json"
    return str(base_dir / filename)


def _load_edits(radius_deg: float | None = None) -> Dict[str, Any]:
    """Load edit links from JSON file with proper error handling."""
    primary_path = Path(_edits_path(radius_deg))
    fallback_path: Path | None = None
    if radius_deg is not None and not primary_path.exists():
        fallback_path = Path(_edits_path(None))
        if fallback_path == primary_path:
            fallback_path = None

    paths_to_try = [p for p in [primary_path, fallback_path] if p is not None]

    for path in paths_to_try:
        if not path.exists():
            continue
        try:
            with open(path, 'r') as f:
                try:
                    return json.load(f)
                except json.JSONDecodeError as e:
                    print(f"Warning: Invalid JSON in {path}: {e}", file=sys.stderr)
                    continue
        except Exception as e:
            print(f"Warning: Error loading {path}: {e}", file=sys.stderr)
            continue

    return {}


def _parse_master_key(master_key: str) -> tuple[str, Any] | None:
    try:
        if ':' not in master_key:
            return None
        kcat, kid = master_key.split(':', 1)
        kcat = kcat.strip()
        if kcat.lower() == 'gaia':
            try:
                return ('gaia_id', int(kid))
            except Exception:
                return ('gaia_id', kid)
        if kcat == '2MASS' or kcat.lower() == '2mass':
            return ('2MASS', kid)
        if kcat.lower() == 'wise':
            return ('wise_id', kid)
        if kcat.lower() == 'chandra':
            return ('chandra_id', kid)
        if kcat.lower() == 'xmm':
            return ('xmm_id', kid)
    except Exception:
        return None
    return None


def _apply_edits_to_matches(combined_df: 'pd.DataFrame',
                            matches: list[list],
                            id_col: str,
                            edits_for_cat: Dict[str, Any]) -> list[list]:
    """Apply remove/force edits to the raw matches list before merging.

    - remove: delete (other, master) assignments if present
    - force: ensure (other -> target master) is present and unique
    """
    if not edits_for_cat:
        return matches

    delete_ids: set[str] = set()
    for entry in (edits_for_cat.get('delete') or []):
        if entry is None:
            continue
        if isinstance(entry, dict):
            val = entry.get('id', entry.get('other'))
        else:
            val = entry
        if val is not None:
            delete_ids.add(str(val))

    # Build list of force assignments allowing multiple masters for the same 'other'
    force_list: list[tuple[Any, int]] = []
    remove_set: set[tuple[int, Any]] = set()

    # Helper to find row indices for a master key
    def _find_master_idx(col: str, val: Any) -> list[int]:
        try:
            if col not in combined_df.columns:
                return []
            if col == 'gaia_id':
                try:
                    return list(combined_df.index[combined_df[col].astype('Int64') == int(val)])
                except Exception:
                    return list(combined_df.index[combined_df[col].astype(str) == str(val)])
            else:
                return list(combined_df.index[combined_df[col].astype(str) == str(val)])
        except Exception:
            return []

    if delete_ids:
        for i in range(len(matches)):
            if matches[i]:
                matches[i] = [oid for oid in matches[i] if str(oid) not in delete_ids]

    # Process removes (support optional row hint to target a specific duplicate master row)
    for rem in edits_for_cat.get('remove', []) or []:
        other = rem.get('other')
        mkey = rem.get('master')
        row_hint = rem.get('row') if isinstance(rem, dict) else None
        if not other or not mkey:
            continue
        spec = None
        # master may be dict {cat_id: value} or string 'cat:value'
        if isinstance(mkey, dict):
            try:
                col, val = next(((k, v) for k, v in mkey.items()))
            except Exception:
                continue
            spec = (col, val)
        elif isinstance(mkey, str):
            spec = _parse_master_key(mkey)
        if not spec:
            continue
        col, val = spec
        # If a specific row is hinted, prefer it
        used = False
        try:
            if row_hint is not None:
                irow = int(row_hint)
                if (irow >= 0) and (irow < len(combined_df)):
                    try:
                        if str(combined_df.at[irow, col]) == str(val):
                            remove_set.add((int(irow), other))
                            used = True
                    except Exception:
                        pass
        except Exception:
            pass
        if not used:
            idxs = _find_master_idx(col, val)
            for i in idxs:
                remove_set.add((int(i), other))

    # Process forces
    for forc in edits_for_cat.get('force', []) or []:
        other = forc.get('other')
        mkey = forc.get('master')
        row_hint = forc.get('row') if isinstance(forc, dict) else None
        if not other or not mkey:
            continue
        spec = None
        if isinstance(mkey, dict):
            try:
                col, val = next(((k, v) for k, v in mkey.items()))
            except Exception:
                continue
            spec = (col, val)
        elif isinstance(mkey, str):
            spec = _parse_master_key(mkey)
        if not spec:
            continue
        col, val = spec
        # If row hint is provided, prefer it when valid and consistent
        idxs: list[int] = []
        try:
            if row_hint is not None:
                irow = int(row_hint)
                if (irow >= 0) and (irow < len(combined_df)):
                    # If the hinted row matches the same master key, use it
                    try:
                        if str(combined_df.at[irow, col]) == str(val):
                            idxs = [irow]
                    except Exception:
                        pass
        except Exception:
            pass
        if not idxs:
            idxs = _find_master_idx(col, val)
        if not idxs:
            continue
        # Pick the first occurrence deterministically (but allow multiple entries overall)
        force_list.append((other, int(idxs[0])))

    # Apply removes
    if remove_set:
        for i in range(len(matches)):
            if not matches[i]:
                continue
            matches[i] = [oid for oid in matches[i] if (i, oid) not in remove_set]

    # Apply forces (additive: do not remove other placements; allow multiple masters)
    if force_list:
        for other, idx in force_list:
            while idx >= len(matches):
                matches.append([])
            row_list = matches[idx] or []
            if other not in row_list:
                row_list.append(other)
            matches[idx] = row_list

    return matches


def apply_link_edits_to_combined(base_df: 'pd.DataFrame',
                                 catalogs: Dict[str, Any],
                                 edits: Dict[str, Any] | None) -> 'pd.DataFrame':
    """Apply force/remove link edits onto a pre-matched combined catalogue.

    Parameters
    ----------
    base_df : DataFrame
        Combined catalogue *before* applying manual link edits. This DataFrame
        is not modified.
    catalogs : dict
        Metadata dictionary returned by the matcher. Only used to keep the
        return signature uniform; the current implementation does not mutate
        it but downstream callers rely on the structure remaining unchanged.
    edits : dict or None
        Mapping loaded from ``links.json`` describing manual edits.

    Returns
    -------
    DataFrame
        A deep copy of ``base_df`` with edits applied. The copy includes a
        ``'_baseline_unforced'`` entry in ``.attrs`` pointing back to
        ``base_df`` so callers can reapply edits without recomputing matches.
    """
    import pandas as _pd
    base_df = base_df if base_df is not None else _pd.DataFrame()
    combined_df = base_df.copy(deep=True)
    combined_df.attrs['_baseline_unforced'] = base_df
    if not isinstance(edits, dict) or combined_df.empty:
        return combined_df

    def _sentinel_for(series: _pd.Series):
        dtype = series.dtype
        try:
            from pandas.api import types as _types
        except Exception:
            _types = None  # type: ignore
        if _types is not None and _types.is_numeric_dtype(dtype):
            return -1
        if getattr(dtype, 'kind', None) in ('i', 'u', 'f'):
            return -1
        return None

    def _coerce_other(value, series: _pd.Series):
        dtype = series.dtype
        try:
            from pandas.api import types as _types
        except Exception:
            _types = None  # type: ignore
        try:
            if _types is not None and _types.is_integer_dtype(dtype):
                return int(value)
            if _types is not None and _types.is_float_dtype(dtype):
                return float(value)
        except Exception:
            pass
        try:
            if getattr(dtype, 'kind', None) in ('i', 'u'):
                return int(value)
            if getattr(dtype, 'kind', None) == 'f':
                return float(value)
        except Exception:
            pass
        return value

    def _resolve_master_rows(df: _pd.DataFrame, column: str, value, row_hint) -> list[int]:
        idxs: list[int] = []
        if column not in df.columns:
            return idxs
        # Honour explicit row hint if it is consistent with the requested master
        if row_hint is not None:
            try:
                irow = int(row_hint)
                if 0 <= irow < len(df):
                    if str(df.at[irow, column]) == str(value):
                        idxs.append(int(irow))
                        return idxs
            except Exception:
                idxs = []
        try:
            mask = df[column].astype(str) == str(value)
            idxs = [int(i) for i in df.index[mask]]
        except Exception:
            idxs = []
        return idxs

    def _section_for(edits_dict: Dict[str, Any], key: str) -> Dict[str, Any]:
        # Accept a few common variants (case-insensitive, wise_id, etc.)
        if key in edits_dict:
            return edits_dict.get(key) or {}
        low = key.lower()
        for cand in edits_dict.keys():
            if str(cand).lower() == low:
                return edits_dict.get(cand) or {}
        aliases = {
            'wise_id': 'wise',
            '2mass': '2MASS',
            '2mass_j': '2MASS',
            '2massj': '2MASS',
            '2mass-j': '2MASS',
        }
        alt = aliases.get(key.lower()) if isinstance(key, str) else None
        if alt and alt in edits_dict:
            return edits_dict.get(alt) or {}
        return {}

    deleted_map: dict[str, set[str]] = {}

    for key in ('2MASS', 'wise', 'chandra', 'xmm'):
        sect = _section_for(edits, key)
        if not sect:
            deleted_map.pop(key, None)
            continue
        id_col = '2MASS' if key == '2MASS' else f"{key}_id"
        if id_col not in combined_df.columns:
            continue
        series = combined_df[id_col]
        sentinel = _sentinel_for(series)
        delete_entries = sect.get('delete', []) or []
        delete_ids: set[str] = set()
        for entry in delete_entries:
            if entry is None:
                continue
            if isinstance(entry, dict):
                val = entry.get('id', entry.get('other'))
            else:
                val = entry
            if val is not None:
                delete_ids.add(str(val))
        if delete_ids:
            deleted_map[key] = set(delete_ids)
        else:
            deleted_map.pop(key, None)
        if delete_ids:
            try:
                mask_del = series.astype(str).isin(delete_ids)
                if mask_del.any():
                    combined_df.loc[mask_del, id_col] = sentinel
            except Exception:
                pass

        # Removals: clear id from the specified master rows
        for rem in sect.get('remove', []) or []:
            if not isinstance(rem, dict):
                continue
            other = rem.get('other')
            master = rem.get('master')
            if not other or not master:
                continue
            spec = _parse_master_key(master) if isinstance(master, str) else None
            if not spec:
                continue
            col, val = spec
            idxs = _resolve_master_rows(combined_df, col, val, rem.get('row'))
            for idx in idxs:
                try:
                    if str(combined_df.at[idx, id_col]) == str(other):
                        combined_df.at[idx, id_col] = sentinel
                except Exception:
                    pass

        # Forces: ensure only the requested master keeps the link
        for forc in sect.get('force', []) or []:
            if not isinstance(forc, dict):
                continue
            other = forc.get('other')
            master = forc.get('master')
            if not other or not master:
                continue
            spec = _parse_master_key(master) if isinstance(master, str) else None
            if not spec:
                continue
            col, val = spec
            idxs = _resolve_master_rows(combined_df, col, val, forc.get('row'))
            if not idxs:
                continue
            target_idx = idxs[0]
            assign_val = _coerce_other(other, series)
            try:
                combined_df.at[target_idx, id_col] = assign_val
            except Exception:
                combined_df.at[target_idx, id_col] = other

        combined_df = _dedupe_unmatched_for_id(combined_df, id_col)
        if delete_ids:
            try:
                mask_drop = (combined_df.get('gaia_id', _pd.Series(dtype=float)) == -1)
                if sentinel is None:
                    mask_drop = mask_drop & combined_df[id_col].isna()
                else:
                    mask_drop = mask_drop & (combined_df[id_col] == sentinel)
                if mask_drop.any():
                    combined_df = combined_df.loc[~mask_drop].copy()
            except Exception:
                pass
    combined_df.reset_index(drop=True, inplace=True)
    if deleted_map:
        combined_df.attrs['_deleted_ids'] = {k: sorted(v) for k, v in deleted_map.items()}
    else:
        combined_df.attrs.pop('_deleted_ids', None)
    return combined_df

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
    # Robustness: ensure finite uncertainties; fall back to small floors if needed
    for _c, _floor in (('errMaj', 0.05), ('errMin', 0.05)):
        try:
            df[_c] = pd.to_numeric(df[_c], errors='coerce')
            df.loc[~np.isfinite(df[_c]), _c] = _floor
            df.loc[df[_c] <= 0, _c] = _floor
        except Exception:
            pass
    try:
        df['errPA'] = pd.to_numeric(df['errPA'], errors='coerce').fillna(0.0)
    except Exception:
        pass
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
    # RUWE filtering disabled per request: keep all Gaia rows regardless of RUWE
    # (previously: RUWE <= 1.4 with NaN pass). If needed, re-enable by applying
    # a mask on df['ruwe'] here.
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
                        r_group_arcsec: float = BLEND_GROUP_RADIUS_ARCSEC,
                        max_pair_sep_arcsec: float = BLEND_MAX_PAIR_SEP_ARCSEC,
                        r_centroid_arcsec: float = BLEND_CENTROID_TOLERANCE_ARCSEC,
                        max_dG_mag: float = BLEND_MAX_MAGNITUDE_DIFF,
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


def chain_other_to_unmatched_master(combined_df: pd.DataFrame,
                                    other_df: pd.DataFrame,
                                    *,
                                    id_col: str,
                                    factor: float,
                                    min_radius_arcsec: float) -> pd.DataFrame:
    """Link 'other' catalog sources directly to unmatched master rows (gaia_id == -1)
    when their positional ellipses overlap (Mahalanobis D^2 <= 1), so that catalogs
    without Gaia counterparts can still be chained together via the master.

    This runs after the normal Gaia-anchored match has been applied. It fills
    ``combined_df[id_col]`` for rows with ``gaia_id == -1`` where a best-match
    from ``other_df`` exists. Conflicts (multiple masters competing for the same
    other id) are resolved by choosing the smallest D^2. Entries already matched
    in the main pass are not modified.
    """
    if combined_df.empty or other_df is None or other_df.empty:
        return combined_df

    # Candidate master rows: synthetic (gaia_id == -1) and missing this id_col
    if id_col not in combined_df.columns:
        combined_df[id_col] = pd.Series([None] * len(combined_df))
    try:
        miss_mask = combined_df[id_col].isna() | (combined_df[id_col] == -1)
    except Exception:
        miss_mask = combined_df[id_col].isna()
    master_mask = miss_mask & (combined_df.get('gaia_id', -1) == -1)
    master_idx = list(np.nonzero(master_mask.values if hasattr(master_mask, 'values') else master_mask)[0])
    if not master_idx:
        return combined_df

    # Precompute other catalog ellipse (scaled and min-clipped) in arcsec
    oth_maj_raw = pd.to_numeric(other_df['errMaj'], errors='coerce').astype(float)
    oth_min_raw = pd.to_numeric(other_df['errMin'], errors='coerce').astype(float)
    oth_pa_deg = pd.to_numeric(other_df['errPA'], errors='coerce').astype(float)
    default_factor = float(factor)
    default_min = float(min_radius_arcsec)
    match_factors = np.full(len(other_df), default_factor, dtype=float)
    if 'match_factor' in other_df.columns:
        mf = pd.to_numeric(other_df['match_factor'], errors='coerce').to_numpy()
        match_factors = np.where(np.isfinite(mf), mf, match_factors)
    match_min = np.full(len(other_df), default_min, dtype=float)
    if 'match_min_radius' in other_df.columns:
        mm = pd.to_numeric(other_df['match_min_radius'], errors='coerce').to_numpy()
        match_min = np.where(np.isfinite(mm), np.maximum(mm, default_min), match_min)
    oth_maj_as = np.where(np.isfinite(oth_maj_raw), oth_maj_raw, default_min)
    oth_min_as = np.where(np.isfinite(oth_min_raw), oth_min_raw, default_min)
    oth_maj_as = np.maximum(oth_maj_as * match_factors, match_min)
    oth_min_as = np.maximum(oth_min_as * match_factors, match_min)

    def _cov_from_axes_arcsec(maj_as, min_as, pa_deg):
        phi = np.deg2rad(pa_deg)
        sphi, cphi = np.sin(phi), np.cos(phi)
        vmaj = np.array([sphi, cphi])
        vmin = np.array([-cphi, sphi])
        return (maj_as*maj_as) * np.outer(vmaj, vmaj) + (min_as*min_as) * np.outer(vmin, vmin)

    # Index other_df by id for quick lookup
    df_other_indexed = other_df.set_index(id_col, drop=False)

    # Build assignment map: other_id -> (best_D2, master_i)
    best_for_other: dict = {}

    for i in master_idx:
        mrow = combined_df.iloc[i]
        ra0 = float(mrow['ra_deg']); dec0 = float(mrow['dec_deg'])
        # Master covariance (already in arcsec in combined_df)
        try:
            m_maj = float(mrow['errMaj']); m_min = float(mrow['errMin']); m_pa = float(mrow['errPA'])
        except Exception:
            continue
        Cmas = _cov_from_axes_arcsec(m_maj, m_min, m_pa)
        # Vectorized residuals to all other sources in arcsec
        dx = (other_df['ra_deg'].values - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
        dy = (other_df['dec_deg'].values - dec0) * 3600.0
        # Quick circular prefilter
        r_pref = (np.maximum(m_maj, m_min) + float(np.nanmax(np.maximum(oth_maj_as, oth_min_as))))
        sel = np.where((dx*dx + dy*dy) <= (r_pref*r_pref))[0]
        for j in sel:
            Coth = _cov_from_axes_arcsec(oth_maj_as[j], oth_min_as[j], float(oth_pa_deg.iloc[j]))
            Csum = Cmas + Coth
            try:
                invC = la.inv(Csum)
            except la.LinAlgError:
                invC = la.inv(Csum + np.eye(2) * 1e-12)
            v = np.array([dx[j], dy[j]])
            d2 = float(v.dot(invC).dot(v))
            if d2 <= 1.0:
                oid = other_df.iloc[j][id_col]
                prev = best_for_other.get(oid)
                if (prev is None) or (d2 < prev[0]):
                    best_for_other[oid] = (d2, i)

    # Apply assignments (only to still-missing rows)
    for oid, (d2, i) in best_for_other.items():
        if pd.isna(combined_df.at[i, id_col]) or (combined_df.at[i, id_col] == -1):
            combined_df.at[i, id_col] = oid
            # Update master ellipse to the smaller-AREA ellipse between the master and this 'other' row
            try:
                other_row = df_other_indexed.loc[oid]
            except Exception:
                other_row = None
            if other_row is not None:
                try:
                    row_factor = float(other_row.get('match_factor', default_factor))
                except Exception:
                    row_factor = default_factor
                try:
                    row_min_radius = float(other_row.get('match_min_radius', default_min))
                except Exception:
                    row_min_radius = default_min
                row_min_radius = max(row_min_radius, default_min)
                try:
                    # Scale other ellipse and enforce min radius
                    maj_o = max(float(other_row['errMaj']) * row_factor, row_min_radius)
                    min_o = max(float(other_row['errMin']) * row_factor, row_min_radius)
                    pa_o  = float(other_row['errPA'])
                except Exception:
                    maj_o = row_min_radius; min_o = row_min_radius; pa_o = 0.0
                try:
                    maj_m = float(combined_df.at[i, 'errMaj']); min_m = float(combined_df.at[i, 'errMin'])
                except Exception:
                    maj_m = float('nan'); min_m = float('nan')
                area_m = maj_m * min_m if (np.isfinite(maj_m) and np.isfinite(min_m)) else float('inf')
                area_o = maj_o * min_o
                if (not np.isfinite(area_m)) or (area_o < area_m):
                    combined_df.at[i, 'ra_deg'] = float(other_row['ra_deg'])
                    combined_df.at[i, 'dec_deg'] = float(other_row['dec_deg'])
                    combined_df.at[i, 'errMaj']  = maj_o
                    combined_df.at[i, 'errMin']  = min_o
                    combined_df.at[i, 'errPA']   = pa_o

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
    n_other = len(other_df)
    factor_default = float(scale_factor)
    min_radius_arcsec_default = float(min_radius.to(u.arcsec).value)
    min_deg_default = min_radius_arcsec_default / 3600.0
    if 'match_factor' in other_df.columns:
        row_factors = pd.to_numeric(other_df['match_factor'], errors='coerce').to_numpy()
        row_factors = np.where(np.isfinite(row_factors), row_factors, factor_default)
    else:
        row_factors = np.full(n_other, factor_default, dtype=float)
    if 'match_min_radius' in other_df.columns:
        mm = pd.to_numeric(other_df['match_min_radius'], errors='coerce').to_numpy()
        row_min_arcsec = np.where(np.isfinite(mm), np.maximum(mm, min_radius_arcsec_default), min_radius_arcsec_default)
    else:
        row_min_arcsec = np.full(n_other, min_radius_arcsec_default, dtype=float)
    row_min_deg = row_min_arcsec / 3600.0

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
    # IMPORTANT: source_df ellipses are already scaled/clipped upstream according to their own
    # catalog policy (e.g., Gaia ×2 and min 0.05"). Do NOT rescale them here with the
    # secondary catalog's factor/min; just convert to degrees and ensure a tiny positive floor.
    eps_deg = 1e-6 / 3600.0
    src_maj = np.maximum((source_df['errMaj']).values / 3600.0, eps_deg)  # deg
    src_min = np.maximum((source_df['errMin']).values / 3600.0, eps_deg)  # deg
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
    oth_maj = np.maximum((other_df['errMaj'].values * row_factors) / 3600.0, row_min_deg)  # deg
    oth_min = np.maximum((other_df['errMin'].values * row_factors) / 3600.0, row_min_deg)  # deg
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
                    # Fall back to using the DataFrame index when no identifier column exists
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
            # Scale only the secondary catalog's ellipse for plotting; the master/source ellipses
            # are already pre-scaled/clipped upstream (e.g., Gaia ×2 and min 0.05").
            oth_maj_plot = np.maximum(other_df['errMaj'].values * row_factors, row_min_arcsec)
            oth_min_plot = np.maximum(other_df['errMin'].values * row_factors, row_min_arcsec)
            src_maj_plot = source_df['errMaj'].values
            src_min_plot = source_df['errMin'].values
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
                                    width=2*oth_maj_plot[k],
                                    height=2*oth_min_plot[k],
                                    angle=(90.0 - float(other_df.iloc[k]['errPA'])),
                                    edgecolor='lightblue', linestyle='-', fill=False)
                        ax.add_patch(e)
                # Plot all master-catalog ellipses in light red
                for m in range(len(source_df)):
                    if abs(dx_src.iloc[m]) <= margin and abs(dy_src.iloc[m]) <= margin:
                        e = Ellipse((dx_src.iloc[m], dy_src.iloc[m]),
                                    width=2*src_maj_plot[m],
                                    height=2*src_min_plot[m],
                                    angle=(90.0 - float(source_df.iloc[m]['errPA'])),
                                    edgecolor='lightcoral', linestyle='-', fill=False)
                        ax.add_patch(e)
                # Overlay and highlight the specific new source (blue solid)
                ell_new = Ellipse((0, 0),
                                  width=2*oth_maj_plot[j],
                                  height=2*oth_min_plot[j],
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
                                try:
                                    _f = float(crow.get('match_factor', all_catalogs[key].get('factor', 1.0)))
                                except Exception:
                                    _f = float(all_catalogs[key].get('factor', 1.0))
                                try:
                                    _mr = float(crow.get('match_min_radius', all_catalogs[key].get('min_radius', 0.0)))
                                except Exception:
                                    _mr = float(all_catalogs[key].get('min_radius', 0.0))
                                _mr = max(_mr, float(all_catalogs[key].get('min_radius', 0.0)))
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
                            # Draw a small colored cross at each overlaid component center for visibility
                            try:
                                _ccol = CAT_COLORS.get(kcat, 'black') if 'CAT_COLORS' in globals() else 'black'
                                ax.plot(dxm, dym, marker='+', color=_ccol, markersize=6,
                                        markeredgewidth=0.9, linestyle='None', alpha=0.9)
                            except Exception:
                                pass
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
        """Final fallback: use CDS HiPS → FITS service.

        Resolution controls (via environment variables):
          - TMASS_ARCSEC_PER_PIX: arcsec/pixel (default 1.0; lower = higher res)
          - TMASS_MAX_SIZE:       max width/height in pixels (default 4096)
        """
        try:
            print("[2MASS J] attempting hips2fits …")
            fov = max(width_deg, height_deg)
            # Controls for output pixel scale and image size
            try:
                arcsec_per_pix = float(str(os.environ.get('TMASS_ARCSEC_PER_PIX', '1.0')))
            except Exception:
                arcsec_per_pix = 1.0
            arcsec_per_pix = float(np.clip(arcsec_per_pix, 0.2, 10.0))
            try:
                max_size = int(str(os.environ.get('TMASS_MAX_SIZE', '4096')))
            except Exception:
                max_size = 4096
            max_size = int(np.clip(max_size, 256, 8192))
            # Compute requested linear size from field-of-view and pixel scale
            size = int(np.clip(round((fov * 3600.0) / arcsec_per_pix), 128, max_size))
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
    factor_default = float(params.get('factor', 1.0))
    min_radius_arcsec_default = float(params.get('min_radius', 0.0))
    row_factors = np.full(len(other_df), factor_default, dtype=float)
    if 'match_factor' in other_df.columns:
        mf = pd.to_numeric(other_df['match_factor'], errors='coerce').to_numpy()
        row_factors = np.where(np.isfinite(mf), mf, row_factors)
    row_min_arcsec = np.full(len(other_df), min_radius_arcsec_default, dtype=float)
    if 'match_min_radius' in other_df.columns:
        mm = pd.to_numeric(other_df['match_min_radius'], errors='coerce').to_numpy()
        row_min_arcsec = np.where(np.isfinite(mm), np.maximum(mm, min_radius_arcsec_default), row_min_arcsec)
    # Color palette per catalog (consistent across plot and table)
    CAT_COLORS = {
        'gaia':    '#1f77b4',  # blue
        '2MASS':   '#2ca02c',  # green
        'wise':    '#ff7f0e',  # orange
        'chandra': '#d62728',  # red
        'xmm':     '#9467bd',  # purple
    }
    # Optional debug overlay for D^2 internals (PDF_DEBUG_D2=1|true|yes)
    try:
        PDF_DEBUG_D2 = str(os.environ.get('PDF_DEBUG_D2', '0')).strip().lower() in ('1', 'true', 'yes', 'on')
    except Exception:
        PDF_DEBUG_D2 = False
    debug_env = str(os.environ.get('ALADIN_SHOW_DEBUG', '1' if PDF_DEBUG_D2 else '0')).strip().lower() in ('1', 'true', 'yes', 'on')
    cur_color = CAT_COLORS.get(key, '#17becf')  # color for the current (new) catalog

    profile_enabled = str(os.environ.get('ALADIN_PROFILE', '0')).strip().lower() in ('1', 'true', 'yes', 'on')
    pdf = PdfPages(pdf_path) if pdf_path else None
    profiler = _PlotProfiler(profile_enabled, prefix=f"plot_after_merge[{key}]")

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
    # Problem pages tracking (persisted across incremental generation)
    import json as _json
    scripts_dir = os.environ.get('ALADIN_SCRIPTS_DIR', 'aladin_scripts')
    html_dir = os.path.join(scripts_dir, 'html')
    try:
        os.makedirs(html_dir, exist_ok=True)
    except Exception:
        pass
    idx_path = os.path.join(html_dir, f"{key}_index.json")
    persisted_problem_pages: set[int] = set()
    if os.path.exists(idx_path):
        try:
            with open(idx_path, 'r') as _f:
                idx_data = _json.load(_f) or {}
            prev_list = idx_data.get('problem_pages')
            if isinstance(prev_list, list):
                for p in prev_list:
                    try:
                        persisted_problem_pages.add(int(p))
                    except Exception:
                        continue
        except Exception as exc:
            logger.warning(f"Could not read existing index {idx_path}: {exc}")
    if not persisted_problem_pages:
        fallback_candidates = [
            os.path.join(html_dir, f"{key}_index_TEST.json"),
            os.path.join(html_dir, f"{key}_index_seed.json"),
        ]
        radius_suffix_env = os.environ.get('ALADIN_PAGE_SUFFIX', '')
        if radius_suffix_env:
            fallback_candidates.extend([
                os.path.join(html_dir, f"{key}{radius_suffix_env}_index_TEST.json"),
                os.path.join(html_dir, f"{key}{radius_suffix_env}_index_seed.json"),
            ])
        for cand in fallback_candidates:
            if not cand or not os.path.exists(cand):
                continue
            try:
                with open(cand, 'r') as _f:
                    idx_data = _json.load(_f) or {}
                fallback_list = idx_data.get('problem_pages')
                if isinstance(fallback_list, list):
                    added = 0
                    for p in fallback_list:
                        try:
                            persisted_problem_pages.add(int(p))
                            added += 1
                        except Exception:
                            continue
                    if added:
                        logger.info(f"Seeded {added} problem pages from fallback index {cand}")
            except Exception as exc:
                logger.warning(f"Failed to load fallback index {cand}: {exc}")
    problem_page_set: set[int] = set(persisted_problem_pages)
    pages_marked_problem: set[int] = set()
    pages_cleared_problem: set[int] = set()

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

    # Strict ID-based lookup of the current catalog's ellipse for the panel row.
    # Returns (maj_scaled_as, min_scaled_as, pa_deg, raw_maj_as, raw_min_as).
    def _scaled_ellipse_by_id(_key: str, _id_val, _id_col: str):
        if not isinstance(all_catalogs, dict):
            raise ValueError("catalog metadata not available")
        _meta = all_catalogs.get(_key)
        if not isinstance(_meta, dict):
            raise ValueError(f"metadata missing for catalog '{_key}'")
        _frame = _meta.get('frame')
        if _frame is None or _id_col not in _frame.columns:
            raise ValueError(f"catalog frame missing or id column '{_id_col}' not present")

        # Robust lookup: accept int/float/str ids; normalize common variants (e.g., '439.0' → '439').
        def _find_row(_frame, _id_col, _id_val, _key):
            trials = []
            # 1) Native value
            trials.append(_id_val)
            # 2) int(value) for numeric strings/floats
            try:
                trials.append(int(float(_id_val)))
            except Exception:
                pass
            # 3) string representation(s)
            s = str(_id_val).strip()
            trials.append(s)
            # Trim trailing '.0' if present (e.g., '439.0' → '439')
            if s.endswith('.0'):
                trials.append(s[:-2])
            # Special-case normalization for 2MASS IDs
            if _key == '2MASS' and (' ' in s):
                trials.append(s.replace(' ', '+'))
                trials.append(s.replace(' ', '-'))
            # De-duplicate while preserving order
            seen = set()
            _uniq = []
            for t in trials:
                if t not in seen:
                    seen.add(t)
                    _uniq.append(t)
            # Try a few matching strategies per trial
            import pandas as _pd
            for t in _uniq:
                # Exact dtype match
                try:
                    rows = _frame.loc[_frame[_id_col] == t]
                    if not rows.empty:
                        return rows.iloc[0]
                except Exception:
                    pass
                # String-equality match
                try:
                    rows = _frame.loc[_frame[_id_col].astype(str) == str(t)]
                    if not rows.empty:
                        return rows.iloc[0]
                except Exception:
                    pass
                # Numeric-equality match
                try:
                    col_num = _pd.to_numeric(_frame[_id_col], errors='coerce')
                    t_num = _pd.to_numeric(_pd.Series([t]), errors='coerce').iloc[0]
                    if _pd.notna(t_num):
                        rows = _frame.loc[col_num == t_num]
                        if not rows.empty:
                            return rows.iloc[0]
                except Exception:
                    pass
            return None

        _r = _find_row(_frame, _id_col, _id_val, _key)
        if _r is None:
            raise ValueError(f"id {_id_val!r} not found in catalog '{_key}' frame")
        try:
            raw_maj = float(_r['errMaj']); raw_min = float(_r['errMin']); pa = float(_r['errPA'])
        except Exception as e:
            raise ValueError(f"invalid raw ellipse for id {_id_val!r}: {e}")
        try:
            _factor = float(_r.get('match_factor', _meta.get('factor', 1.0)))
        except Exception:
            _factor = float(_meta.get('factor', 1.0))
        try:
            _minr = float(_r.get('match_min_radius', _meta.get('min_radius', 0.0)))
        except Exception:
            _minr = float(_meta.get('min_radius', 0.0))
        maj_sc = max(raw_maj * _factor, _minr)
        min_sc = max(raw_min * _factor, _minr)
        return maj_sc, min_sc, pa, raw_maj, raw_min

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

    profiler.tic('precompute')
    per_page = ncols * nrows
    total = len(other_df)
    import math as _math
    total_pages = max(1, _math.ceil(total / max(1, per_page)))

    # Precompute new-catalog ellipse axes for plotting
    if plot_mode == 'match':
        oth_maj_plot = np.maximum(other_df['errMaj'].values * row_factors, row_min_arcsec)
        oth_min_plot = np.maximum(other_df['errMin'].values * row_factors, row_min_arcsec)
    else:
        oth_maj_plot = other_df['errMaj'].values
        oth_min_plot = other_df['errMin'].values
    profiler.toc('precompute')

    # Precompute tangent-plane coordinates for efficient per-panel subsets.
    if len(combined_df):
        ra_ref = float(np.nanmedian(combined_df['ra_deg']))
        dec_ref = float(np.nanmedian(combined_df['dec_deg']))
    elif len(other_df):
        ra_ref = float(np.nanmedian(other_df['ra_deg']))
        dec_ref = float(np.nanmedian(other_df['dec_deg']))
    else:
        ra_ref = 0.0
        dec_ref = 0.0
    cos_ref = np.cos(np.deg2rad(dec_ref))

    combined_ra_arr = combined_df['ra_deg'].to_numpy(dtype=float, copy=False)
    combined_dec_arr = combined_df['dec_deg'].to_numpy(dtype=float, copy=False)
    combined_x = (combined_ra_arr - ra_ref) * cos_ref * 3600.0
    combined_y = (combined_dec_arr - dec_ref) * 3600.0

    other_ra_arr = other_df['ra_deg'].to_numpy(dtype=float, copy=False)
    other_dec_arr = other_df['dec_deg'].to_numpy(dtype=float, copy=False)
    other_x = (other_ra_arr - ra_ref) * cos_ref * 3600.0
    other_y = (other_dec_arr - dec_ref) * 3600.0

    try:
        cell_arcsec = float(os.environ.get('ALADIN_PANEL_GRID_ARCSEC', '20'))
    except Exception:
        cell_arcsec = 20.0
    if not np.isfinite(cell_arcsec) or cell_arcsec <= 0.0:
        cell_arcsec = 20.0
    inv_cell = 1.0 / cell_arcsec

    def _build_grid(x_arr: np.ndarray, y_arr: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
        grid: dict[tuple[int, int], list[int]] = defaultdict(list)
        if x_arr.size == 0:
            return {}
        cell_x = np.floor(x_arr * inv_cell).astype(np.int64, copy=False)
        cell_y = np.floor(y_arr * inv_cell).astype(np.int64, copy=False)
        for idx, (cx, cy) in enumerate(zip(cell_x, cell_y)):
            grid[(int(cx), int(cy))].append(idx)
        return {key: np.asarray(values, dtype=np.int64) for key, values in grid.items()}

    def _collect_indices(grid: dict[tuple[int, int], np.ndarray],
                         center_x: float,
                         center_y: float,
                         radius_arcsec: float) -> np.ndarray:
        if not grid or radius_arcsec <= 0.0:
            return np.empty(0, dtype=np.int64)
        cell_radius = int(np.ceil(radius_arcsec / cell_arcsec)) + 1
        base_cx = int(np.floor(center_x * inv_cell))
        base_cy = int(np.floor(center_y * inv_cell))
        buckets: list[np.ndarray] = []
        for dx_cell in range(-cell_radius, cell_radius + 1):
            cx = base_cx + dx_cell
            for dy_cell in range(-cell_radius, cell_radius + 1):
                arr = grid.get((cx, base_cy + dy_cell))
                if arr is not None and arr.size:
                    buckets.append(arr)
        if not buckets:
            return np.empty(0, dtype=np.int64)
        return np.unique(np.concatenate(buckets))

    combined_grid = _build_grid(combined_x, combined_y)
    other_grid = _build_grid(other_x, other_y)

    if 'gaia_id' in combined_df.columns:
        valid_master_mask_global = (~combined_df['gaia_id'].isna()) & (combined_df['gaia_id'] != -1)
    else:
        valid_master_mask_global = np.ones(len(combined_df), dtype=bool)
    if id_col in combined_df.columns:
        id_series_global = combined_df[id_col]
    elif len(combined_df.columns):
        id_series_global = combined_df.iloc[:, 0]
    else:
        id_series_global = pd.Series([None] * len(combined_df))
    try:
        id_series_global_str = id_series_global.astype(str).to_numpy(copy=False)
    except Exception:
        id_series_global_str = np.asarray([str(v) for v in id_series_global])

    # Allow restricting generation to a single page via env: ALADIN_PAGE_ONLY
    try:
        _page_only = int(str(os.environ.get('ALADIN_PAGE_ONLY', '0')))
    except Exception:
        _page_only = 0

    logger.info(f"Starting page generation for catalog '{key}': total={total}, per_page={per_page}, "
                f"num_pages={math.ceil(total/per_page)}, aladin_dir={'provided' if aladin_dir else 'None'}")

    deleted_map_attr = combined_df.attrs.get('_deleted_ids', {})
    if isinstance(deleted_map_attr, dict):
        deleted_ids_global = {k: {str(v) for v in (vals or [])} for k, vals in deleted_map_attr.items()}
    else:
        deleted_ids_global = {}

    for start in range(0, total, per_page):
        profiler.tic('page_total')
        page_num = (start // per_page) + 1
        if _page_only > 0 and page_num != _page_only:
            logger.debug(f"Skipping page {page_num} (ALADIN_PAGE_ONLY={_page_only})")
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
            deleted_ids = deleted_ids_global.get(key, set())
            is_deleted = str(new_id) in deleted_ids
            ra0 = other_df.iloc[j]['ra_deg']
            dec0 = other_df.iloc[j]['dec_deg']
            epoch_row = _epoch_from_row(other_df.iloc[j])

            profiler.tic('panel_offsets')
            center_x = (ra0 - ra_ref) * cos_ref * 3600.0
            center_y = (dec0 - dec_ref) * 3600.0

            margin = max(oth_maj_plot[j], oth_min_plot[j]) * 2.0
            margin = max(margin, 5.0)  # at least ±2.5" half-width
            search_radius_src = margin + cell_arcsec
            search_radius_oth = margin + cell_arcsec

            oth_idx_arr = _collect_indices(other_grid, center_x, center_y, search_radius_oth)
            src_idx_arr = _collect_indices(combined_grid, center_x, center_y, search_radius_src)

            def _build_offset_dict(idx_array: np.ndarray,
                                   x_arr: np.ndarray,
                                   y_arr: np.ndarray) -> dict[int, tuple[float, float]]:
                offsets: dict[int, tuple[float, float]] = {}
                if idx_array.size:
                    dx_vals = x_arr[idx_array] - center_x
                    dy_vals = y_arr[idx_array] - center_y
                    for idx_val, dx_val, dy_val in zip(idx_array, dx_vals, dy_vals):
                        offsets[int(idx_val)] = (float(dx_val), float(dy_val))
                return offsets

            oth_offsets = _build_offset_dict(oth_idx_arr, other_x, other_y)
            if j not in oth_offsets:
                oth_offsets[int(j)] = (float(other_x[j] - center_x), float(other_y[j] - center_y))
            src_offsets = _build_offset_dict(src_idx_arr, combined_x, combined_y)

            def _oth_offset(idx: int) -> tuple[float, float]:
                idx = int(idx)
                val = oth_offsets.get(idx)
                if val is None:
                    val = (float(other_x[idx] - center_x), float(other_y[idx] - center_y))
                    oth_offsets[idx] = val
                return val

            def _src_offset(idx: int) -> tuple[float, float]:
                idx = int(idx)
                val = src_offsets.get(idx)
                if val is None:
                    val = (float(combined_x[idx] - center_x), float(combined_y[idx] - center_y))
                    src_offsets[idx] = val
                return val

            def _src_dx(idx: int) -> float:
                return _src_offset(idx)[0]

            def _src_dy(idx: int) -> float:
                return _src_offset(idx)[1]

            oth_indices = list(oth_offsets.keys())
            new_id_str = str(new_id)

            candidate_indices_raw = list(src_offsets.keys())
            _nearest_idx = None
            nearest_sep_arcsec = None
            if candidate_indices_raw:
                cand_arr = np.asarray(candidate_indices_raw, dtype=int)
                mask_valid = valid_master_mask_global[cand_arr]
                mask_not_self = (id_series_global_str[cand_arr] != new_id_str)
                nearby_candidates = cand_arr[mask_valid & mask_not_self]
                if nearby_candidates.size:
                    try:
                        gaia_df = all_catalogs.get('gaia', {}).get('data') if isinstance(all_catalogs, dict) else None
                    except Exception:
                        gaia_df = None
                    best = None
                    for ii in nearby_candidates:
                        if (gaia_df is not None) and ('gaia_id' in combined_df.columns):
                            gid = combined_df.iloc[ii].get('gaia_id', None)
                            if (gid is not None) and (gid in gaia_df.index):
                                grow = gaia_df.loc[gid]
                                dxm = (float(grow['ra_deg']) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                                dym = (float(grow['dec_deg']) - dec0) * 3600.0
                                dist = float(np.hypot(dxm, dym))
                            else:
                                dxm, dym = _src_offset(ii)
                                dist = float(np.hypot(dxm, dym))
                        else:
                            dxm, dym = _src_offset(ii)
                            dist = float(np.hypot(dxm, dym))
                        if (best is None) or (dist < best[1]):
                            best = (ii, dist)
                    if best is not None:
                        _nearest_idx = int(best[0])
                        nearest_sep_arcsec = float(best[1])

            cand_set = set(candidate_indices_raw)
            if _nearest_idx is not None:
                cand_set.add(int(_nearest_idx))
            cand_indices = [idx for idx in cand_set
                            if abs(_src_dx(idx)) <= margin and abs(_src_dy(idx)) <= margin]
            profiler.toc('panel_offsets')
            # PM-corrected closest separation (computed later when we know Δt); used for title display
            closest_sep_pm_arcsec = None

            # If this is a 2MASS panel and allowed, draw a 2MASS J-band background aligned via WCS
            profiler.tic('panel_tmass_fetch')
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
            profiler.toc('panel_tmass_fetch')

            profiler.toc('panel_offsets')
            profiler.tic('panel_overlay')
            # Plot new-catalog background (light blue)
            for k in oth_indices:
                dx_k, dy_k = _oth_offset(k)
                if abs(dx_k) <= margin and abs(dy_k) <= margin:
                    ls_bg = 'dotted' if (is_deleted and k == j) else '-'
                    e = Ellipse((dx_k, dy_k),
                                width=2*oth_maj_plot[k],
                                height=2*oth_min_plot[k],
                                angle=(90.0 - float(other_df.iloc[k]['errPA'])),
                                edgecolor='lightblue', linestyle=ls_bg, fill=False)
                    ax.add_patch(e)

            # Do not hide the synthetic self-row anymore in master-centric mode.
            # Previously we suppressed rows with (same id and gaia_id == -1). That prevented
            # legitimate matches in catalogs without Gaia counterpart from appearing.
            # Keep all candidate rows; they will be marked as matched via master_indices.

            # Overlay the specific new source using strict ID-based ellipse lookup
            try:
                _maj_use, _min_use, _pa_use, _raw_maj_use, _raw_min_use = _scaled_ellipse_by_id(key, other_df.iloc[j][id_col], id_col)
                _w_new = 2 * _maj_use
                _h_new = 2 * _min_use
                _ang_new = (90.0 - _pa_use)
            except Exception:
                # If lookup by ID fails, do not silently substitute; use current row values
                row_factor_fallback = float(row_factors[j]) if j < len(row_factors) else float(all_catalogs[key]['factor'])
                row_min_fallback = float(row_min_arcsec[j]) if j < len(row_min_arcsec) else float(all_catalogs[key]['min_radius'])
                _w_new = 2 * max(float(other_df.iloc[j]['errMaj']) * row_factor_fallback, row_min_fallback)
                _h_new = 2 * max(float(other_df.iloc[j]['errMin']) * row_factor_fallback, row_min_fallback)
                _ang_new = (90.0 - float(other_df.iloc[j]['errPA']))
            ls_new = 'dotted' if is_deleted else '-'
            e0 = Ellipse((0, 0), width=_w_new, height=_h_new,
                         angle=_ang_new, edgecolor=cur_color, linestyle=ls_new, fill=False,
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
                    if master_idx not in cand_set:
                        _src_offset(master_idx)
                        cand_set.add(master_idx)
                        cand_indices.append(master_idx)
                    mrow = combined_df.iloc[master_idx]
                    # Skip synthetic self-row for the current catalog (unmatched "self")
                    try:
                        _same_id = (str(mrow.get(id_col, '')) == str(new_id))
                    except Exception:
                        _same_id = False
                    try:
                        _is_synth = ('gaia_id' in combined_df.columns) and (
                            pd.isna(mrow.get('gaia_id', np.nan)) or mrow.get('gaia_id', -1) == -1
                        )
                    except Exception:
                        _is_synth = False
                    if _same_id and _is_synth:
                        continue
                    # Draw the master position (diamond) only if it has at least one preceding-catalog id
                    try:
                        mx, my = _src_offset(master_idx)
                        def _has_prev(mr, cur_key: str) -> bool:
                            order = ['gaia', '2MASS', 'wise', 'chandra', 'xmm']
                            try:
                                pos = order.index(cur_key)
                            except ValueError:
                                pos = 0
                            prev = set(order[:pos])
                            for kcat, col in [('gaia','gaia_id'), ('2MASS','2MASS'), ('wise','wise_id'), ('chandra','chandra_id'), ('xmm','xmm_id')]:
                                if (kcat in prev) and (col in combined_df.columns) and (not _is_missing(mr.get(col, None))):
                                    return True
                            return False
                        if _has_prev(mrow, key):
                            ax.plot(mx, my, marker='D', linestyle='None', markersize=5,
                                    markeredgewidth=1.0, markerfacecolor='none', color='#c0c0c0',
                                    alpha=0.95, zorder=12)
                            has_master_diamond = True
                    except Exception:
                        pass
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
                                try:
                                    _f = float(crow.get('match_factor', all_catalogs[kcat].get('factor', 1.0)))
                                except Exception:
                                    _f = float(all_catalogs[kcat].get('factor', 1.0))
                                try:
                                    _mr = float(crow.get('match_min_radius', all_catalogs[kcat].get('min_radius', 0.0)))
                                except Exception:
                                    _mr = float(all_catalogs[kcat].get('min_radius', 0.0))
                                _mr = max(_mr, float(all_catalogs[kcat].get('min_radius', 0.0)))
                                _maj_m = max(float(crow['errMaj']) * _f, _mr)
                                _min_m = max(float(crow['errMin']) * _f, _mr)
                            else:
                                _maj_m = float(crow['errMaj'])
                                _min_m = float(crow['errMin'])
                            lw = 1.8 if kcat == 'gaia' else 1.2
                            ell_m = Ellipse((dxm, dym), width=2*_maj_m, height=2*_min_m,
                                            angle=(90.0 - float(crow['errPA'])),
                                             edgecolor=CAT_COLORS.get(kcat, 'red'), linestyle='-',
                                             linewidth=lw, fill=False, label=None)
                            ax.add_patch(ell_m)
                            # Draw center cross for every matched component
                            try:
                                ax.plot(dxm, dym, marker='+', color=CAT_COLORS.get(kcat, 'black'),
                                        markersize=6, markeredgewidth=0.9, linestyle='None', alpha=0.9)
                            except Exception:
                                pass
                            # Emphasize matched GAIA position with a colored plus + PM arrow/inflated ellipse
                            if kcat == 'gaia':
                                ax.plot(dxm, dym, marker='+', color=CAT_COLORS['gaia'],
                                        markersize=7, markeredgewidth=1.2, linestyle='None')
                                # Always draw a PM arrow when PM information is available
                                try:
                                    if (epoch_row is not None
                                        and ('pmra' in crow.index) and ('pmdec' in crow.index)
                                        and ('pmra_error' in crow.index) and ('pmdec_error' in crow.index)
                                        and ('ref_epoch' in crow.index)):
                                        dt = float(epoch_row - float(crow['ref_epoch']))
                                        pm_dx = float(crow.get('pmra', 0.0)) * 0.001 * dt   # arcsec along RA*
                                        pm_dy = float(crow.get('pmdec', 0.0)) * 0.001 * dt  # arcsec Dec
                                        if np.isfinite(pm_dx) and np.isfinite(pm_dy):
                                            ax.arrow(dxm, dym, pm_dx, pm_dy,
                                                     width=0.0, head_width=0.03,
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
            # and draw a thin gray cross at every master position. A row is considered
            # matched if it belongs to master_indices, regardless of whether it has a Gaia id.
            matched_set = set(master_indices)
            for m in cand_indices:
                if m not in matched_set:
                    # Only draw diamonds for unmatched rows when there is at least one matched row
                    # (i.e., to provide context). Additionally, never draw the diamond for the
                    # synthetic self-row of the current catalog source (unmatched "self").
                    if master_indices:
                        try:
                            same_id = (str(combined_df.iloc[m][id_col]) == str(new_id))
                        except Exception:
                            same_id = False
                        try:
                            is_synth = ('gaia_id' in combined_df.columns) and (
                                pd.isna(combined_df.iloc[m].get('gaia_id', np.nan)) or
                                combined_df.iloc[m].get('gaia_id', -1) == -1
                            )
                        except Exception:
                            is_synth = False
                        if not (same_id and is_synth):
                            # Only draw if the unmatched row has at least one preceding-catalog id
                            order = ['gaia', '2MASS', 'wise', 'chandra', 'xmm']
                            try:
                                pos = order.index(key)
                            except ValueError:
                                pos = 0
                            prev = set(order[:pos])
                            mrow2 = combined_df.iloc[m]
                            has_prev = False
                            for kcat, col in [('gaia','gaia_id'), ('2MASS','2MASS'), ('wise','wise_id'), ('chandra','chandra_id'), ('xmm','xmm_id')]:
                                if (kcat in prev) and (col in combined_df.columns) and (not _is_missing(mrow2.get(col, None))):
                                    has_prev = True
                                    break
                            if has_prev:
                                mx_u, my_u = _src_offset(m)
                                ax.plot(
                                    mx_u, my_u,
                                    marker='D', linestyle='None',
                                    markersize=5, markeredgewidth=1.0,
                                    markerfacecolor='none', color='#c0c0c0', alpha=0.9
                                )
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
                        try:
                            _f = float(crow.get('match_factor', all_catalogs[kcat].get('factor', 1.0)))
                        except Exception:
                            _f = float(all_catalogs[kcat].get('factor', 1.0))
                        try:
                            _mr = float(crow.get('match_min_radius', all_catalogs[kcat].get('min_radius', 0.0)))
                        except Exception:
                            _mr = float(all_catalogs[kcat].get('min_radius', 0.0))
                        _mr = max(_mr, float(all_catalogs[kcat].get('min_radius', 0.0)))
                        _maj_m = max(float(crow['errMaj']) * _f, _mr)
                        _min_m = max(float(crow['errMin']) * _f, _mr)
                    else:
                        _maj_m = float(crow['errMaj'])
                        _min_m = float(crow['errMin'])
                    ell_bg = Ellipse((dxm, dym), width=2*_maj_m, height=2*_min_m,
                                     angle=(90.0 - float(crow['errPA'])), edgecolor=CAT_COLORS.get(kcat, 'red'), linestyle='-',
                                     linewidth=0.8, fill=False, alpha=0.28)
                    ax.add_patch(ell_bg)
                    # Add a small colored cross at the center of each background component
                    try:
                        ax.plot(dxm, dym, marker='+', color=CAT_COLORS.get(kcat, 'black'),
                                markersize=6, markeredgewidth=0.9, linestyle='None', alpha=0.9)
                    except Exception:
                        pass

            profiler.toc('panel_overlay')
            profiler.tic('panel_diagnostics')
            # --- Compute compact diagnostics (best separation and Mahalanobis D²) for title ---
            def _cov_from_axes(maj_as, min_as, pa_deg):
                # Clamp semi-axes to a small positive value and build covariance in arcsec^2
                eps = 1e-6
                try:
                    a = float(maj_as)
                    b = float(min_as)
                    pa = float(pa_deg)
                except Exception:
                    a = np.nan; b = np.nan; pa = 0.0
                if not np.isfinite(a):
                    a = eps
                if not np.isfinite(b):
                    b = eps
                if not np.isfinite(pa):
                    pa = 0.0
                a = max(a, eps)
                b = max(b, eps)
                phi = np.deg2rad(pa)
                sphi, cphi = np.sin(phi), np.cos(phi)
                vmaj = np.array([sphi, cphi])
                vmin = np.array([-cphi, sphi])
                C = (a*a) * np.outer(vmaj, vmaj) + (b*b) * np.outer(vmin, vmin)
                # Ensure strictly positive-definite by adding tiny jitter
                C = C + np.eye(2) * (eps * eps)
                # Replace any remaining non-finite numbers
                if not np.all(np.isfinite(C)):
                    C = np.nan_to_num(C, nan=eps, posinf=1e6, neginf=1e6)
                return C

            def _safe_inv(C):
                M = np.array(C, dtype=float)
                if not np.all(np.isfinite(M)):
                    M = np.nan_to_num(M, nan=0.0, posinf=1e6, neginf=1e6)
                # Symmetrize and add jitter
                M = 0.5 * (M + M.T) + np.eye(2) * 1e-9
                try:
                    return la.inv(M)
                except Exception:
                    try:
                        return la.inv(M + np.eye(2) * 1e-6)
                    except Exception:
                        try:
                            return la.pinv(M)
                        except Exception:
                            return None

            # Pretty-printer for 2x2 matrices that also sanitizes any stray NaNs/Infs
            def _fmt_mat(M):
                try:
                    A = np.array(M, dtype=float)
                except Exception:
                    A = np.asarray(M)
                A = np.nan_to_num(A, nan=0.0, posinf=1e6, neginf=-1e6)
                return np.array2string(A, formatter={'float_kind': lambda x: f"{x:.6f}"})

            # Compute D^2 using exactly the same math as in the matched case
            # (vector from combined_df offsets, same ellipse scaling/flooring and PM inflation)
            def _compute_d2_with_master_idx(m_idx: int):
                try:
                    # 2MASS ellipse used for matching (scaled + min_radius)
                    maj_use, min_use, pa_use, _, _ = _scaled_ellipse_by_id(key, other_df.iloc[j][id_col], id_col)
                    Cnew_loc = _cov_from_axes(maj_use, min_use, pa_use)
                    # Master ellipse from combined_df with 0.05" floor
                    mrow_loc = combined_df.iloc[int(m_idx)]
                    try:
                        maj_m = float(mrow_loc['errMaj']); min_m = float(mrow_loc['errMin']); pa_m = float(mrow_loc['errPA'])
                    except Exception:
                        maj_m = 0.05; min_m = 0.05; pa_m = 0.0
                    if not np.isfinite(maj_m): maj_m = 0.05
                    if not np.isfinite(min_m): min_m = 0.05
                    maj_m = max(maj_m, 0.05); min_m = max(min_m, 0.05)
                    Cmas_loc = _cov_from_axes(maj_m, min_m, pa_m)
                    # PM-uncertainty inflation (same as matched block)
                    if (epoch_row is not None) and ('ref_epoch' in combined_df.columns):
                        try:
                            ref_ep_local = float(mrow_loc.get('ref_epoch', np.nan))
                            if np.isfinite(ref_ep_local) and ('pmra_error' in combined_df.columns) and ('pmdec_error' in combined_df.columns):
                                dt_tmp_local = float(epoch_row - ref_ep_local)
                                sra = float(mrow_loc.get('pmra_error', 0.0))
                                sde = float(mrow_loc.get('pmdec_error', 0.0))
                                if np.isfinite(sra) and np.isfinite(sde):
                                    sig_ra = abs(sra) * 0.001 * abs(dt_tmp_local)
                                    sig_de = abs(sde) * 0.001 * abs(dt_tmp_local)
                                    Cmas_loc = Cmas_loc + np.diag([sig_ra**2, sig_de**2])
                        except Exception:
                            pass
                    Csum_loc = Cnew_loc + Cmas_loc
                    invC_loc = _safe_inv(Csum_loc)
                    if invC_loc is None:
                        return None
                    # PM-corrected residual vector using combined_df offsets (like the matched case)
                    pm_dx = pm_dy = 0.0
                    if (epoch_row is not None) and ('ref_epoch' in combined_df.columns):
                        try:
                            ref_ep = float(mrow_loc.get('ref_epoch', np.nan))
                            if np.isfinite(ref_ep) and ('pmra' in combined_df.columns) and ('pmdec' in combined_df.columns):
                                dt_tmp2 = float(epoch_row - ref_ep)
                                _pmra2 = float(mrow_loc.get('pmra', 0.0)); _pmde2 = float(mrow_loc.get('pmdec', 0.0))
                                if not np.isfinite(_pmra2): _pmra2 = 0.0
                                if not np.isfinite(_pmde2): _pmde2 = 0.0
                                pm_dx = _pmra2 * 0.001 * dt_tmp2
                                pm_dy = _pmde2 * 0.001 * dt_tmp2
                        except Exception:
                            pm_dx = pm_dy = 0.0
                    dx_loc, dy_loc = _src_offset(int(m_idx))
                    vec_loc = np.array([dx_loc + pm_dx, dy_loc + pm_dy])
                    d2_loc = float(vec_loc.dot(invC_loc).dot(vec_loc))
                    sep_loc = float(np.hypot(*vec_loc))
                    try:
                        dt_loc = float(epoch_row - float(mrow_loc['ref_epoch'])) if (epoch_row is not None) else None
                    except Exception:
                        dt_loc = None
                    return {
                        'D2': d2_loc,
                        'sep': sep_loc,
                        'dt': dt_loc,
                        'Cnew': Cnew_loc,
                        'Cmas': Cmas_loc,
                        'Csum': Csum_loc,
                        'vec': vec_loc,
                        'maj_m': maj_m, 'min_m': min_m, 'pa_m': pa_m,
                        'maj_use': maj_use, 'min_use': min_use, 'pa_use': pa_use,
                    }
                except Exception:
                    return None

            # Ensure the sum covariance is not unrealistically small by enforcing
            # a lower bound on its smallest eigenvalue. Physically, Cnew + Cmas is
            # positive semidefinite with λ_min >= min(λ_min(Cnew), λ_min(Cmas)).
            # Numerical mishaps can still yield a too‑small λ_min; we guard against
            # that by adding a tiny diagonal so that λ_min >= 0.5 * min(sigma_min^2).
            def _enforce_min_eig(Csum, min_eig_bound_as2):
                try:
                    w, _ = np.linalg.eigh(Csum)
                    bound = float(min_eig_bound_as2)
                    if not np.isfinite(bound):
                        return Csum
                    if w[0] < bound:
                        delta = float(bound - w[0])
                        if np.isfinite(delta) and (delta > 0):
                            Csum = Csum + np.eye(2) * delta
                    return Csum
                except Exception:
                    return Csum

            best_sep = None
            best_d2 = None
            best_dt = None
            best_master_maj = None
            best_master_min = None
            best_master_pa  = None
            used_master_idx = None
            nearest_d2 = None
            nearest_dt = None
            debug_info_lines = [] if debug_env else None
            diag_log: list[tuple[str, float]] | None = [] if debug_env else None
            diag_total_start = perf_counter()

            if master_indices:
                # Consider only real master rows (exclude synthetic self-rows with gaia_id == -1)
                try:
                    master_valid = [int(m) for m in master_indices
                                    if ('gaia_id' in combined_df.columns)
                                    and (not pd.isna(combined_df.iloc[int(m)].get('gaia_id', np.nan)))
                                    and (combined_df.iloc[int(m)].get('gaia_id', -1) != -1)]
                except Exception:
                    master_valid = []
                if master_valid:
                    with _diag_timer(profiler, diag_log, 'matched_branch'):
                        _maj_use, _min_use, _pa_use, _raw_maj_use, _raw_min_use = _scaled_ellipse_by_id(key, other_df.iloc[j][id_col], id_col)
                        Cnew = _cov_from_axes(_maj_use, _min_use, _pa_use)
                        master_floor_as = 0.05

                        def _area_at(idx):
                            try:
                                a = float(combined_df.iloc[idx]['errMaj'])
                                b = float(combined_df.iloc[idx]['errMin'])
                            except Exception:
                                a = np.nan; b = np.nan
                            if not np.isfinite(a):
                                a = master_floor_as
                            if not np.isfinite(b):
                                b = master_floor_as
                            a = max(a, master_floor_as)
                            b = max(b, master_floor_as)
                            return a * b

                        try:
                            used_master_idx = min(master_valid, key=_area_at)
                        except Exception:
                            used_master_idx = master_valid[0]
                        m = used_master_idx
                        pm_dx = 0.0
                        pm_dy = 0.0
                        if (epoch_row is not None) and ('pmra' in combined_df.columns) and ('pmdec' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                            try:
                                ref_ep = float(combined_df.iloc[m]['ref_epoch'])
                                if np.isfinite(ref_ep):
                                    dt_tmp = float(epoch_row - ref_ep)
                                    _pmra = float(combined_df.iloc[m]['pmra'])
                                    _pmde = float(combined_df.iloc[m]['pmdec'])
                                    if not np.isfinite(_pmra):
                                        _pmra = 0.0
                                    if not np.isfinite(_pmde):
                                        _pmde = 0.0
                                    pm_dx = _pmra * 0.001 * dt_tmp
                                    pm_dy = _pmde * 0.001 * dt_tmp
                            except Exception:
                                pm_dx = 0.0
                                pm_dy = 0.0
                        dx_m, dy_m = _src_offset(m)
                        vec = np.array([dx_m + pm_dx, dy_m + pm_dy])
                        mrow = combined_df.iloc[m]
                        try:
                            maj_m = float(mrow['errMaj'])
                            min_m = float(mrow['errMin'])
                            pa_m = float(mrow['errPA'])
                        except Exception:
                            maj_m = master_floor_as
                            min_m = master_floor_as
                            pa_m = 0.0
                        if not np.isfinite(maj_m):
                            maj_m = master_floor_as
                        if not np.isfinite(min_m):
                            min_m = master_floor_as
                        maj_m = max(maj_m, master_floor_as)
                        min_m = max(min_m, master_floor_as)
                        Cmas = _cov_from_axes(maj_m, min_m, pa_m)
                        try:
                            if (not np.all(np.isfinite(Cmas))) or (Cmas[0, 0] < 1e-6) or (Cmas[1, 1] < 1e-6):
                                phi = np.deg2rad(float(pa_m)) if np.isfinite(pa_m) else 0.0
                                sphi, cphi = np.sin(phi), np.cos(phi)
                                a2 = float(maj_m) * float(maj_m)
                                b2 = float(min_m) * float(min_m)
                                C11 = a2 * sphi * sphi + b2 * cphi * cphi
                                C22 = a2 * cphi * cphi + b2 * sphi * sphi
                                C12 = (a2 - b2) * sphi * cphi
                                Cmas = np.array([[C11, C12], [C12, C22]], dtype=float) + np.eye(2) * 1e-9
                        except Exception:
                            pass
                        if (epoch_row is not None) and ('ref_epoch' in combined_df.columns):
                            try:
                                ref_ep_local = float(mrow['ref_epoch'])
                                if np.isfinite(ref_ep_local) and ('pmra_error' in combined_df.columns) and ('pmdec_error' in combined_df.columns):
                                    dt_tmp_local = float(epoch_row - ref_ep_local)
                                    sra = float(mrow['pmra_error'])
                                    sde = float(mrow['pmdec_error'])
                                    if np.isfinite(sra) and np.isfinite(sde):
                                        sig_ra = abs(sra) * 0.001 * abs(dt_tmp_local)
                                        sig_de = abs(sde) * 0.001 * abs(dt_tmp_local)
                                        Cmas = Cmas + np.diag([sig_ra**2, sig_de**2])
                            except Exception:
                                pass
                        Csum = Cnew + Cmas
                        invC = _safe_inv(Csum)
                        if invC is None:
                            best_d2 = None
                            best_sep = None
                            if debug_env and debug_info_lines is not None:
                                debug_info_lines.append("D2: invC None (singular)")
                        else:
                            d2 = float(vec.dot(invC).dot(vec))
                            sep = float(np.hypot(*vec))
                            best_d2, best_sep = d2, sep
                            if epoch_row is not None and ('ref_epoch' in combined_df.columns):
                                try:
                                    best_dt = float(epoch_row - float(mrow['ref_epoch']))
                                except Exception:
                                    best_dt = None
                            best_master_maj = maj_m
                            best_master_min = min_m
                            best_master_pa = pa_m
                            if debug_env and debug_info_lines is not None:
                                try:
                                    evals, _ = np.linalg.eigh(Csum)
                                    evals = [float(x) for x in evals]
                                except Exception:
                                    evals = []
                                debug_info_lines.append(
                                    f"D2 master={used_master_idx} vec=({vec[0]:.4f},{vec[1]:.4f}) sep={sep:.4f}"
                                )
                                debug_info_lines.append(
                                    f"  Cnew: (maj={oth_maj_plot[j]:.2f}, min={oth_min_plot[j]:.2f}, pa={float(other_df.iloc[j]['errPA']):.1f})"
                                )
                                try:
                                    debug_info_lines.append("    Cnew_mat=" + _fmt_mat(Cnew))
                                except Exception:
                                    pass
                                debug_info_lines.append(
                                    f"  Cmas: (maj={maj_m:.2f}, min={min_m:.2f}, pa={pa_m:.1f}) eig={evals} D2={d2:.6f}"
                                )
                                try:
                                    debug_info_lines.append("    Cmas_mat=" + _fmt_mat(Cmas))
                                except Exception:
                                    pass
                else:
                    # No valid master rows despite a synthetic self-row: compute nearest D² like the unmatched branch
                    try:
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
                            try:
                                gaia_df2 = all_catalogs.get('gaia', {}).get('data') if isinstance(all_catalogs, dict) else None
                            except Exception:
                                gaia_df2 = None
                            try:
                                nearest_mask_arr = nearest_mask.to_numpy(dtype=bool, copy=False)  # type: ignore[attr-defined]
                            except Exception:
                                nearest_mask_arr = np.asarray(nearest_mask, dtype=bool)
                            idxs = np.nonzero(nearest_mask_arr)[0]
                            if idxs.size:
                                dx_vals = combined_x[idxs] - center_x
                                dy_vals = combined_y[idxs] - center_y
                                dist2_vals = dx_vals * dx_vals + dy_vals * dy_vals
                                best_pos = int(np.argmin(dist2_vals))
                                best_idx = int(idxs[best_pos])
                                _nearest_idx = best_idx
                                nearest_sep_arcsec = float(np.sqrt(dist2_vals[best_pos]))
                                # Cache the offset for downstream consumers that expect it to exist.
                                src_offsets.setdefault(best_idx, (float(dx_vals[best_pos]), float(dy_vals[best_pos])))
                    except Exception:
                        pass
                    # There is a synthetic self-row but no valid Gaia master: compute D² to the nearest valid master
                    mnear = None
                    if _nearest_idx is not None:
                        mnear = int(_nearest_idx)
                    else:
                        cand_for_nearest = np.array(list(cand_set), dtype=int) if cand_set else np.empty(0, dtype=int)
                        if cand_for_nearest.size:
                            mask_valid = valid_master_mask_global[cand_for_nearest]
                            mask_not_self = (id_series_global_str[cand_for_nearest] != new_id_str)
                            filtered = cand_for_nearest[mask_valid & mask_not_self]
                        else:
                            filtered = np.empty(0, dtype=int)
                        with _diag_timer(profiler, diag_log, 'nearest_fallback'):
                            if filtered.size:
                                dist_all = np.array([np.hypot(*_src_offset(int(idx))) for idx in filtered], dtype=float)
                                k = int(np.argmin(dist_all))
                                mnear = int(filtered[k])
                            elif cand_indices:
                                dist_all = np.array([np.hypot(*_src_offset(int(idx))) for idx in cand_indices], dtype=float)
                                k = int(np.argmin(dist_all))
                                mnear = int(cand_indices[k])
                            else:
                                raise RuntimeError('nearest index not found')
                            _maj_use, _min_use, _pa_use, _raw_maj_use, _raw_min_use = _scaled_ellipse_by_id(key, other_df.iloc[j][id_col], id_col)
                            Cnew = _cov_from_axes(_maj_use, _min_use, _pa_use)
                            mrow = combined_df.iloc[mnear]
                            try:
                                maj_m = float(mrow['errMaj'])
                                min_m = float(mrow['errMin'])
                                pa_m  = float(mrow['errPA'])
                            except Exception:
                                maj_m = master_floor_as; min_m = master_floor_as; pa_m = 0.0
                            if not np.isfinite(maj_m): maj_m = master_floor_as
                            if not np.isfinite(min_m): min_m = master_floor_as
                            maj_m = max(maj_m, master_floor_as)
                            min_m = max(min_m, master_floor_as)
                            Cmas = _cov_from_axes(maj_m, min_m, pa_m)
                            if (epoch_row is not None) and ('ref_epoch' in combined_df.columns):
                                dt_tmp_local = None
                                try:
                                    ref_ep_local = float(mrow['ref_epoch'])
                                    if np.isfinite(ref_ep_local):
                                        dt_tmp_local = float(epoch_row - ref_ep_local)
                                except Exception:
                                    dt_tmp_local = None
                                if (dt_tmp_local is not None) and ('pmra_error' in combined_df.columns) and ('pmdec_error' in combined_df.columns):
                                    try:
                                        sra = float(mrow['pmra_error'])
                                        sde = float(mrow['pmdec_error'])
                                        if np.isfinite(sra) and np.isfinite(sde):
                                            sig_ra = abs(sra) * 0.001 * abs(dt_tmp_local)
                                            sig_de = abs(sde) * 0.001 * abs(dt_tmp_local)
                                            Cmas = Cmas + np.diag([sig_ra**2, sig_de**2])
                                    except Exception:
                                        pass
                            Csum = Cnew + Cmas
                            # Guard: λ_min(Cnew+Cmas) ≥ λ_min(Cnew) + λ_min(Cmas) = min(maj_new,min_new)^2 + min(maj_m,min_m)^2
                            try:
                                _min_eig_bound = float(min(oth_maj_plot[j], oth_min_plot[j])**2 + min(maj_m, min_m)**2)
                                Csum = _enforce_min_eig(Csum, _min_eig_bound)
                            except Exception:
                                pass
                        invC = _safe_inv(Csum)
                        if invC is None:
                            nearest_d2 = None
                            closest_sep_pm_arcsec = None
                            continue
                        # PM-corrected residual vector to this nearest master
                        pm_dx = 0.0; pm_dy = 0.0
                        if (epoch_row is not None) and ('pmra' in combined_df.columns) and ('pmdec' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                            try:
                                ref_ep = float(mrow['ref_epoch'])
                                if np.isfinite(ref_ep):
                                    dt_tmp = float(epoch_row - ref_ep)
                                    _pmra = float(mrow['pmra'])
                                    _pmde = float(mrow['pmdec'])
                                    if not np.isfinite(_pmra):
                                        _pmra = 0.0
                                    if not np.isfinite(_pmde):
                                        _pmde = 0.0
                                    pm_dx = _pmra * 0.001 * dt_tmp
                                    pm_dy = _pmde * 0.001 * dt_tmp
                            except Exception:
                                pm_dx = 0.0; pm_dy = 0.0
                        # Residual vector w.r.t. original Gaia position when available
                        try:
                            gaia_df2 = all_catalogs.get('gaia', {}).get('data') if isinstance(all_catalogs, dict) else None
                        except Exception:
                            gaia_df2 = None
                        if (gaia_df2 is not None) and ('gaia_id' in combined_df.columns):
                            gid2 = combined_df.iloc[mnear].get('gaia_id', None)
                            if (gid2 is not None) and (gid2 in gaia_df2.index):
                                g2 = gaia_df2.loc[gid2]
                                dxm = (float(g2['ra_deg']) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                                dym = (float(g2['dec_deg']) - dec0) * 3600.0
                                vec = np.array([dxm + pm_dx, dym + pm_dy])
                            else:
                                dx_near, dy_near = _src_offset(mnear)
                                vec = np.array([dx_near + pm_dx, dy_near + pm_dy])
                        else:
                            dx_near2, dy_near2 = _src_offset(mnear)
                            vec = np.array([dx_near2 + pm_dx, dy_near2 + pm_dy])
                        nearest_d2 = float(vec.dot(invC).dot(vec))
                        closest_sep_pm_arcsec = float(np.hypot(*vec))
                        try:
                            nearest_dt = float(epoch_row - float(mrow['ref_epoch'])) if (epoch_row is not None) else None
                        except Exception:
                            nearest_dt = None
                        best_master_maj = maj_m
                        best_master_min = min_m
                        best_master_pa  = pa_m
                        if debug_env and debug_info_lines is not None:
                            try:
                                evals, _ = np.linalg.eigh(Csum)
                                evals = [float(x) for x in evals]
                            except Exception:
                                evals = []
                            try:
                                _gid_dbg = mrow.get('gaia_id', None)
                            except Exception:
                                _gid_dbg = None
                            debug_info_lines.append(
                                f"D2 nearest master={mnear} gaia_id={_gid_dbg} vec=({vec[0]:.4f},{vec[1]:.4f}) sep={closest_sep_pm_arcsec:.4f}"
                            )
                            debug_info_lines.append(
                                f"  Cnew: (maj={oth_maj_plot[j]:.2f}, min={oth_min_plot[j]:.2f}, pa={float(other_df.iloc[j]['errPA']):.1f})"
                            )
                            try:
                                debug_info_lines.append("    Cnew_mat=" + _fmt_mat(Cnew))
                            except Exception:
                                pass
                            debug_info_lines.append(
                                f"  Cmas: (maj={maj_m:.2f}, min={min_m:.2f}, pa={pa_m:.1f}) eig={evals} D2={nearest_d2 if nearest_d2 is not None else float('nan'):.6f}"
                            )
                            try:
                                debug_info_lines.append("    Cmas_mat=" + _fmt_mat(Cmas))
                            except Exception:
                                pass
            else:
                # No master rows with the same id at all: compute D² to the nearest valid master by angular separation
                if _nearest_idx is not None:
                    mnear = int(_nearest_idx)
                else:
                    cand_for_nearest = np.array(list(cand_set), dtype=int) if cand_set else np.empty(0, dtype=int)
                    if cand_for_nearest.size:
                        mask_valid = valid_master_mask_global[cand_for_nearest]
                        mask_not_self = (id_series_global_str[cand_for_nearest] != new_id_str)
                        filtered = cand_for_nearest[mask_valid & mask_not_self]
                    else:
                        filtered = np.empty(0, dtype=int)
                    if filtered.size:
                        dist_all = np.array([np.hypot(*_src_offset(int(idx))) for idx in filtered], dtype=float)
                        k = int(np.argmin(dist_all))
                        mnear = int(filtered[k])
                    elif cand_indices:
                        dist_all = np.array([np.hypot(*_src_offset(int(idx))) for idx in cand_indices], dtype=float)
                        k = int(np.argmin(dist_all))
                        mnear = int(cand_indices[k])
                    else:
                        raise RuntimeError('nearest index not found')
                    _maj_use, _min_use, _pa_use, _raw_maj_use, _raw_min_use = _scaled_ellipse_by_id(key, other_df.iloc[j][id_col], id_col)
                    Cnew = _cov_from_axes(_maj_use, _min_use, _pa_use)
                    mrow = combined_df.iloc[mnear]
                    try:
                        maj_m = float(mrow['errMaj'])
                        min_m = float(mrow['errMin'])
                        pa_m  = float(mrow['errPA'])
                    except Exception:
                        maj_m = 0.05; min_m = 0.05; pa_m = 0.0
                    if not np.isfinite(maj_m): maj_m = 0.05
                    if not np.isfinite(min_m): min_m = 0.05
                    maj_m = max(maj_m, 0.05); min_m = max(min_m, 0.05)
                    Cmas = _cov_from_axes(maj_m, min_m, pa_m)
                    try:
                        if (not np.all(np.isfinite(Cmas))) or (Cmas[0,0] < 1e-6) or (Cmas[1,1] < 1e-6):
                            phi = np.deg2rad(float(pa_m)) if np.isfinite(pa_m) else 0.0
                            sphi, cphi = np.sin(phi), np.cos(phi)
                            a2 = float(maj_m) * float(maj_m)
                            b2 = float(min_m) * float(min_m)
                            C11 = a2 * sphi*sphi + b2 * cphi*cphi
                            C22 = a2 * cphi*cphi + b2 * sphi*sphi
                            C12 = (a2 - b2) * sphi * cphi
                            Cmas = np.array([[C11, C12], [C12, C22]], dtype=float) + np.eye(2) * 1e-9
                    except Exception:
                        pass
                    if (epoch_row is not None) and ('pmra_error' in combined_df.columns) and ('pmdec_error' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                        try:
                            ref_ep = float(mrow['ref_epoch']); dt_tmp = float(epoch_row - ref_ep)
                            sig_ra = abs(float(mrow['pmra_error'])) * 0.001 * abs(dt_tmp)
                            sig_de = abs(float(mrow['pmdec_error'])) * 0.001 * abs(dt_tmp)
                            Cmas = Cmas + np.diag([sig_ra**2, sig_de**2])
                        except Exception:
                            pass
                    Csum = Cnew + Cmas
                    # Guard: λ_min(Cnew+Cmas) ≥ λ_min(Cnew) + λ_min(Cmas)
                    try:
                        _min_eig_bound = float(min(_maj_use, _min_use)**2 + min(maj_m, min_m)**2)
                        Csum = _enforce_min_eig(Csum, _min_eig_bound)
                    except Exception:
                        pass
                    invC = _safe_inv(Csum)
                    # PM-corrected residual vector to this nearest master
                    pm_dx = 0.0; pm_dy = 0.0
                    if (epoch_row is not None) and ('pmra' in combined_df.columns) and ('pmdec' in combined_df.columns) and ('ref_epoch' in combined_df.columns):
                        try:
                            ref_ep = float(mrow['ref_epoch'])
                            if np.isfinite(ref_ep):
                                dt_tmp = float(epoch_row - ref_ep)
                                _pmra = float(mrow['pmra'])
                                _pmde = float(mrow['pmdec'])
                                if not np.isfinite(_pmra):
                                    _pmra = 0.0
                                if not np.isfinite(_pmde):
                                    _pmde = 0.0
                                pm_dx = _pmra * 0.001 * dt_tmp
                                pm_dy = _pmde * 0.001 * dt_tmp
                        except Exception:
                            pm_dx = 0.0; pm_dy = 0.0
                    # Residual vector using Gaia position if available
                    try:
                        gaia_df3 = all_catalogs.get('gaia', {}).get('data') if isinstance(all_catalogs, dict) else None
                    except Exception:
                        gaia_df3 = None
                    if (gaia_df3 is not None) and ('gaia_id' in combined_df.columns):
                        gid3 = combined_df.iloc[mnear].get('gaia_id', None)
                        if (gid3 is not None) and (gid3 in gaia_df3.index):
                            g3 = gaia_df3.loc[gid3]
                            dxm3 = (float(g3['ra_deg']) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                            dym3 = (float(g3['dec_deg']) - dec0) * 3600.0
                            vec = np.array([dxm3 + pm_dx, dym3 + pm_dy])
                        else:
                            dx_gaia_fallback, dy_gaia_fallback = _src_offset(mnear)
                            vec = np.array([dx_gaia_fallback + pm_dx, dy_gaia_fallback + pm_dy])
                    else:
                        dx_gaia_else, dy_gaia_else = _src_offset(mnear)
                        vec = np.array([dx_gaia_else + pm_dx, dy_gaia_else + pm_dy])
                    nearest_d2 = float(vec.dot(invC).dot(vec))
                    closest_sep_pm_arcsec = float(np.hypot(*vec))
                    try:
                        nearest_dt = float(epoch_row - float(mrow['ref_epoch'])) if (epoch_row is not None) else None
                    except Exception:
                        nearest_dt = None
                    best_master_maj = maj_m
                    best_master_min = min_m
                    best_master_pa  = pa_m

            if debug_env:
                total_diag = (perf_counter() - diag_total_start) if diag_total_start is not None else 0.0
                print(f"[diag] {key} panel {page_num}:{j} total={total_diag*1000:.1f} ms")
                if diag_log:
                    for name, dt in diag_log:
                        print(f"[diag]    {name:20s} {dt*1000:.1f} ms")
            profiler.toc('panel_diagnostics')
            cand_indices = list(dict.fromkeys(cand_indices))
            profiler.tic('panel_table')
            # === Dynamic association table across available catalogs (exclude current catalog) ===
            # Determine which catalogs to include based on the candidate rows; always include
            # all catalogs that precede the current one in the merge order used for the PDF/merge
            # (gaia → 2MASS → wise → chandra → xmm). This ensures, e.g., that GAIA appears in
            # the WISE panel even if none of the sources in the panel have GAIA counterparts.
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
            # Ensure we include all catalogs up to (but excluding) the current one in merge order
            order = ['gaia', '2MASS', 'wise', 'chandra', 'xmm']
            try:
                pos = order.index(key)
            except ValueError:
                pos = 0
            for k_prev in order[:pos]:
                include.add(k_prev)
            display_cats = [k for k, _ in cat_col_map if (k in include and k != key)]

            # Build rows sorted by distance from the current catalog source.
            # Use Gaia positions (master frame) when available to avoid shifts
            # introduced by earlier ellipse-area updates of combined_df coords.
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
                        if (gid is not None) and (not _is_missing(gid)) and (gid in _gaia_df.index):
                            grow = _gaia_df.loc[gid]
                            dxm = (float(grow['ra_deg']) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                            dym = (float(grow['dec_deg']) - dec0) * 3600.0
                            dist_as = float(np.hypot(dxm, dym))
                        else:
                            dist_as = float(np.hypot(*_src_offset(m)))
                    else:
                        dist_as = float(np.hypot(*_src_offset(m)))
                except Exception:
                    dist_as = float(np.hypot(*_src_offset(m)))
                cand_with_dist.append((int(m), dist_as))
            cand_with_dist.sort(key=lambda t: t[1])

            def _fmt(val):
                return str(val) if val is not None else '—'

            # Build table columns and headers; draw G magnitude near the right panel edge (not in table)
            col_keys = list(display_cats)
            add_gmag = ('gaia' in display_cats) and ('Gmag' in combined_df.columns)
            headers = [k.upper() for k in col_keys]
            raw_rows = []  # list of (cells_list, is_matched, master_idx)
            g_values = []  # per-row G values (strings) to render at right edge
            rows_limited = []  # rows actually rendered in the table
            toggle_order = []
            row_positions = {}
            # De-duplicate table rows by a stable master key (prefer GAIA id when available)
            seen_keys_tbl: set = set()
            for m, _dist in cand_with_dist:
                # Skip the synthetic self-row for the current catalog (unmatched self only)
                skip_self = False
                try:
                    same_id = (str(combined_df.iloc[m][id_col]) == str(new_id))
                except Exception:
                    same_id = False
                try:
                    is_synth = ('gaia_id' in combined_df.columns) and (pd.isna(combined_df.iloc[m].get('gaia_id', np.nan)) or combined_df.iloc[m].get('gaia_id', -1) == -1)
                except Exception:
                    is_synth = False
                if same_id and is_synth:
                    skip_self = True
                if skip_self:
                    continue
                mrow = combined_df.iloc[m]
                # Build de-duplication key: prefer GAIA id; else fall back to tuple of displayed columns
                try:
                    if ('gaia_id' in combined_df.columns) and (not _is_missing(mrow.get('gaia_id', None))):
                        dedup_key = ('gaia', str(int(mrow['gaia_id'])))
                    else:
                        dedup_key = tuple(str(mrow.get(col_for[k], '')) for k in col_keys)
                except Exception:
                    dedup_key = tuple(str(mrow.get(col_for[k], '')) for k in col_keys)
                if dedup_key in seen_keys_tbl:
                    continue
                seen_keys_tbl.add(dedup_key)
                cells = []
                for kcat in col_keys:
                    colname = col_for[kcat]
                    if colname in combined_df.columns and not _is_missing(mrow.get(colname, None)):
                        cells.append(_fmt(_short_of(kcat, mrow[colname])))
                    else:
                        cells.append(_fmt(None))
                # Determine matched state robustly by direct equality on this row
                is_matched = False
                try:
                    is_matched = (str(combined_df.iloc[m][id_col]) == str(new_id))
                except Exception:
                    is_matched = (m in matched_set)
                raw_rows.append((cells, is_matched, int(m)))
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

            # Skip printing the table if it is effectively empty (all entries are '—')
            all_empty = True
            try:
                for cells, _is_matched, _m_idx in raw_rows:
                    if any(str(c).strip() != '—' for c in cells):
                        all_empty = False
                        break
            except Exception:
                all_empty = False

            if headers and not all_empty:
                # Compute column widths from headers and rows (in characters)
                table_data = [headers] + [cells for (cells, _, _) in raw_rows]
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
                has_match_row = any(is_m for (_, is_m, _) in raw_rows)
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
                for r_idx, (cells, is_matched, master_idx) in enumerate(rows_limited):
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
                    row_positions[int(master_idx)] = r_idx

                # If there are more rows than fit, add a small ellipsis line
                if len(rows_limited) < len(raw_rows):
                    y = y0 - (len(rows_limited) + 1) * line_h
                    trans = offset_copy(ax.transAxes, fig=fig, x=x_offsets[0], y=0, units='points')
                    ax.text(0.05, y, f"… (+{len(raw_rows) - len(rows_limited)} more)", transform=trans,
                            va='top', ha='left', family='monospace', fontsize=font_size, color='gray')

                if row_positions:
                    toggle_order = [idx for idx, _ in sorted(row_positions.items(), key=lambda kv: kv[1])]

            profiler.toc('panel_table')

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
            if (best_sep is not None) and (best_d2 is not None) and np.isfinite(best_sep) and np.isfinite(best_d2):
                # Use Unicode double-prime for arcsec and superscript 2 for D²
                title += f" — sep={best_sep:.2f}″, D²={best_d2:.2f}"
                if (best_dt is not None) and np.isfinite(best_dt):
                    title += f", Δt={best_dt:.1f}y"
            elif (nearest_sep_arcsec is not None) or (closest_sep_pm_arcsec is not None):
                # No matches: report separation and D² for the nearest master source
                # Recompute D² robustly in arcsec units using the same ellipses as matching
                # (do this unconditionally to avoid carrying forward any earlier coarse values).
                try:
                        # 1) Find nearest master index if not available
                        mnear = None
                        if ('_nearest_idx' in locals()) and (_nearest_idx is not None):
                            mnear = int(_nearest_idx)
                        else:
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
                                # Prefer GAIA original positions when available
                                gaia_df2 = all_catalogs.get('gaia', {}).get('data') if isinstance(all_catalogs, dict) else None
                                best = None
                                for ii in idxs:
                                    if (gaia_df2 is not None) and ('gaia_id' in combined_df.columns):
                                        gid = combined_df.iloc[ii].get('gaia_id', None)
                                        if (gid is not None) and (gid in gaia_df2.index):
                                            grow = gaia_df2.loc[gid]
                                            dxm = (float(grow['ra_deg']) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                                            dym = (float(grow['dec_deg']) - dec0) * 3600.0
                                            d = float(np.hypot(dxm, dym))
                                        else:
                                            d = float(np.hypot(*_src_offset(ii)))
                                    else:
                                        d = float(np.hypot(*_src_offset(ii)))
                                    if (best is None) or (d < best[1]):
                                        best = (ii, d)
                                if best is not None:
                                    mnear = int(best[0])
                        if mnear is None:
                            raise RuntimeError('nearest index not found')
                        # Compute using the exact same math as the matched case
                        try:
                            _res = _compute_d2_with_master_idx(mnear)
                        except Exception:
                            _res = None
                        if _res is not None:
                            nearest_d2 = float(_res['D2'])
                            try:
                                closest_sep_pm_arcsec = float(_res['sep'])
                            except Exception:
                                pass
                            try:
                                nearest_dt = float(_res['dt']) if (_res['dt'] is not None) else nearest_dt
                            except Exception:
                                pass
                            if debug_env and debug_info_lines is not None:
                                try:
                                    evals, _ = np.linalg.eigh(_res['Csum'])
                                    evals = [float(x) for x in evals]
                                except Exception:
                                    evals = []
                                try:
                                    _gid_dbg = combined_df.iloc[mnear].get('gaia_id', None)
                                except Exception:
                                    _gid_dbg = None
                                debug_info_lines.append(
                                    f"D2 nearest (recomp-matchmath) master={mnear} gaia_id={_gid_dbg} vec=({_res['vec'][0]:.4f},{_res['vec'][1]:.4f}) sep={_res['sep']:.4f}"
                                )
                                debug_info_lines.append(
                                    f"  Cnew_as: (maj={_res['maj_use']:.2f}, min={_res['min_use']:.2f}, pa={_res['pa_use']:.1f})"
                                )
                                debug_info_lines.append("    Cnew_as_mat=" + _fmt_mat(_res['Cnew']))
                                debug_info_lines.append(
                                    f"  Cmas_as: (maj={_res['maj_m']:.2f}, min={_res['min_m']:.2f}, pa={_res['pa_m']:.1f})"
                                )
                                debug_info_lines.append("    Cmas_as_mat=" + _fmt_mat(_res['Cmas']))
                                debug_info_lines.append(
                                    f"  Csum_as eig={evals} D2={nearest_d2:.6f}"
                                )
                        try:
                            nearest_dt = float(_res['dt']) if (_res['dt'] is not None) else nearest_dt
                        except Exception:
                            pass
                        # Ensure Gaia master overlay (cross + PM arrow) is visible even when not matched
                        try:
                            gaia_cat = all_catalogs.get('gaia', {}).get('data') if isinstance(all_catalogs, dict) else None
                        except Exception:
                            gaia_cat = None
                        gid_overlay = None
                        try:
                            if 'gaia_id' in combined_df.columns:
                                val_overlay = combined_df.iloc[mnear].get('gaia_id', None)
                                if not _is_missing(val_overlay):
                                    gid_overlay = int(val_overlay)
                        except Exception:
                            gid_overlay = None
                        if (gaia_cat is not None) and (gid_overlay is not None) and (gid_overlay in gaia_cat.index):
                            key_gaia = ('gaia', gid_overlay)
                            if key_gaia not in drawn_components:
                                crow = gaia_cat.loc[gid_overlay]
                                dxm = (float(crow['ra_deg']) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                                dym = (float(crow['dec_deg']) - dec0) * 3600.0
                                try:
                                    ax.plot(dxm, dym, marker='+', color=CAT_COLORS['gaia'],
                                            markersize=7, markeredgewidth=1.2, linestyle='None', zorder=11)
                                except Exception:
                                    pass
                                drawn_components.add(key_gaia)
                                try:
                                    short_val = _short_of('gaia', gid_overlay)
                                    if short_val is not None and ('gaia', short_val) not in seen_ids:
                                        ids_parts.append(f"gaia:{short_val}")
                                        seen_ids.add(('gaia', short_val))
                                except Exception:
                                    pass
                                pm_dx = pm_dy = 0.0
                                if epoch_row is not None and ('pmra' in crow.index) and ('pmdec' in crow.index) and ('ref_epoch' in crow.index):
                                    try:
                                        dt_arrow = float(epoch_row - float(crow['ref_epoch']))
                                        pm_dx = float(crow.get('pmra', 0.0)) * 0.001 * dt_arrow
                                        pm_dy = float(crow.get('pmdec', 0.0)) * 0.001 * dt_arrow
                                    except Exception:
                                        pm_dx = pm_dy = 0.0
                                arrow_len = float(np.hypot(pm_dx, pm_dy))
                                if arrow_len > 0.0:
                                    try:
                                        ax.arrow(dxm, dym, pm_dx, pm_dy,
                                                 width=0.0, head_width=0.03, head_length=0.05,
                                                 length_includes_head=True,
                                                 color=CAT_COLORS['gaia'], alpha=0.9, zorder=11)
                                    except Exception:
                                        pass
                        # Fallback legacy path removed
                except Exception:
                    # Preserve previous value but log the error for visibility
                    if debug_env and debug_info_lines is not None:
                        import traceback as _tb
                        debug_info_lines.append("D2 nearest error: " + _tb.format_exc(limit=1).strip())
                if (nearest_d2 is not None) and np.isfinite(nearest_d2):
                    sep_to_show = (closest_sep_pm_arcsec
                                   if closest_sep_pm_arcsec is not None else nearest_sep_arcsec)
                    if sep_to_show is not None and np.isfinite(sep_to_show):
                        title += f" — nearest: sep={sep_to_show:.2f}″, D²={nearest_d2:.2f}"
                    else:
                        title += f" — nearest: D²={nearest_d2:.2f}"
                    if (nearest_dt is not None) and np.isfinite(nearest_dt):
                        title += f", Δt={nearest_dt:.1f}y"
                else:
                    sep_to_show = (closest_sep_pm_arcsec
                                   if closest_sep_pm_arcsec is not None else nearest_sep_arcsec)
                    if sep_to_show is not None and np.isfinite(sep_to_show):
                        title += f" — nearest: sep={sep_to_show:.2f}″"
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
            close_unmatched = 0
            close_unmatched_any = 0
            try:
                if master_indices:
                    for _m in master_indices:
                        skip_row = False
                        try:
                            skip_row = ('gaia_id' in combined_df.columns) and (
                                pd.isna(combined_df.iloc[_m].get('gaia_id', np.nan)) or
                                combined_df.iloc[_m].get('gaia_id', -1) == -1
                            )
                        except Exception:
                            skip_row = False
                        if skip_row:
                            continue
                        d_arcsec = float(np.hypot(*_src_offset(_m)))
                        if d_arcsec <= near_thr:
                            close_count += 1
                            try:
                                id_val = combined_df.iloc[_m].get(id_col, None) if id_col in combined_df.columns else None
                                if _is_missing(id_val):
                                    close_unmatched += 1
                            except Exception:
                                pass
            except Exception:
                pass
            try:
                for _cand in cand_indices:
                    try:
                        dx_c, dy_c = _src_offset(_cand)
                    except Exception:
                        continue
                    if np.hypot(dx_c, dy_c) <= near_thr:
                        try:
                            if str(combined_df.iloc[_cand].get(id_col, '')) == str(new_id):
                                continue
                        except Exception:
                            pass
                        try:
                            id_val = combined_df.iloc[_cand].get(id_col, None) if id_col in combined_df.columns else None
                            if _is_missing(id_val):
                                close_unmatched_any += 1
                        except Exception:
                            continue
            except Exception:
                pass
            is_problem = False
            # Case A: no identifications, but a nearby master is present
            if (num_valid_master == 0 and (near_sep is not None) and (near_sep <= near_thr)):
                is_problem = True
            # Case B: an identification exists and there are multiple nearby masters
            elif (num_valid_master >= 1 and close_count >= 2):
                is_problem = True
            # Case D: an identification exists but a nearby master appears unmatched (missing id)
            elif (num_valid_master >= 1 and close_unmatched >= 1):
                is_problem = True
            # Case E: any nearby master candidate (regardless of match list) is currently unmatched
            elif close_unmatched_any >= 1:
                is_problem = True
            # Case C: a computed Mahalanobis distance exists for the chosen match and exceeds the 1-sigma threshold
            try:
                if (not is_problem) and (best_d2 is not None) and np.isfinite(best_d2) and (best_d2 > 1.0):
                    is_problem = True
            except Exception:
                pass
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
            # Annotate ellipse parameters for master and current catalog (used values) only when debugging
            if debug_env:
                try:
                    new_pa = float(other_df.iloc[j]['errPA'])
                except Exception:
                    new_pa = 0.0
                # Ensure master ellipse values exist for annotation; if missing, derive from nearest GAIA
                if (best_master_maj is None) or (best_master_min is None) or (best_master_pa is None):
                    try:
                        # Prefer the previously chosen nearest index if available
                        mnear_ann = None
                        if '_nearest_idx' in locals() and (_nearest_idx is not None):
                            mnear_ann = int(_nearest_idx)
                        else:
                            # Recompute nearest using Gaia positions
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
                                # Prefer original Gaia coords when available
                                gaia_dfA = None
                                try:
                                    gaia_dfA = all_catalogs.get('gaia', {}).get('data') if isinstance(all_catalogs, dict) else None
                                except Exception:
                                    gaia_dfA = None
                                best = None
                                for ii in idxs:
                                    if (gaia_dfA is not None) and ('gaia_id' in combined_df.columns):
                                        gid = combined_df.iloc[ii].get('gaia_id', None)
                                        if (gid is not None) and (gid in gaia_dfA.index):
                                            grow = gaia_dfA.loc[gid]
                                            dxm = (float(grow['ra_deg']) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
                                            dym = (float(grow['dec_deg']) - dec0) * 3600.0
                                            d = float(np.hypot(dxm, dym))
                                        else:
                                            d = float(np.hypot(*_src_offset(ii)))
                                    else:
                                        d = float(np.hypot(*_src_offset(ii)))
                                    if (best is None) or (d < best[1]):
                                        best = (ii, d)
                                if best is not None:
                                    mnear_ann = int(best[0])
                        if mnear_ann is not None:
                            mrowA = combined_df.iloc[mnear_ann]
                            try:
                                best_master_maj = float(mrowA['errMaj'])
                                best_master_min = float(mrowA['errMin'])
                                best_master_pa  = float(mrowA['errPA'])
                            except Exception:
                                best_master_maj = best_master_min = None; best_master_pa = None
                    except Exception:
                        pass
                try:
                    txt_master = (f"M: maj={best_master_maj:.2f}″ min={best_master_min:.2f}″ PA={best_master_pa:.1f}°"
                                  if (best_master_maj is not None and best_master_min is not None and best_master_pa is not None)
                                  else "M: maj=— min=— PA=—")
                except Exception:
                    txt_master = "M: maj=— min=— PA=—"
                # Show both native (catalog) uncertainties and the scaled/clipped values used for matching
                try:
                    # Prefer fetching raw (catalog-native) uncertainties by ID to avoid any index drift
                    id_val = other_df.iloc[j][id_col]
                    raw_maj = raw_min = None
                    try:
                        _cat = all_catalogs.get(key, {}) if isinstance(all_catalogs, dict) else {}
                        _frame = _cat.get('frame') if isinstance(_cat, dict) else None
                        if _frame is not None and id_col in _frame.columns:
                            _row = _frame.loc[_frame[id_col] == id_val]
                            if len(_row) > 0:
                                raw_maj = float(_row.iloc[0]['errMaj'])
                                raw_min = float(_row.iloc[0]['errMin'])
                    except Exception:
                        raw_maj = raw_min = None
                    if raw_maj is None or raw_min is None:
                        # Fallback to positional iloc if lookup by id failed
                        raw_maj = float(other_df.iloc[j]['errMaj'])
                        raw_min = float(other_df.iloc[j]['errMin'])
                    # Scale by per-catalog factor and clip to per-catalog min_radius (arcsec)
                    row_factor_dbg = float(row_factors[j]) if j < len(row_factors) else factor_default
                    row_min_dbg = float(row_min_arcsec[j]) if j < len(row_min_arcsec) else min_radius_arcsec_default
                    maj_sc = max(raw_maj * row_factor_dbg, row_min_dbg)
                    min_sc = max(raw_min * row_factor_dbg, row_min_dbg)
                    txt_new = (
                        f"{key}: maj={maj_sc:.2f}″ min={min_sc:.2f}″ PA={new_pa:.1f}° "
                        f"(raw: {raw_maj:.2f}/{raw_min:.2f})"
                    )
                except Exception:
                    txt_new = f"{key}: maj=— min=— PA=—"
                try:
                    # Also note which master catalogs were used for D² (if applicable)
                    if used_master_idx is not None:
                        try:
                            mrow = combined_df.iloc[int(used_master_idx)]
                            parts = []
                            for kcat, colname, label in [
                                ('gaia','gaia_id','gaia'),
                                ('2MASS','2MASS','2MASS'),
                                ('wise','wise_id','wise'),
                                ('chandra','chandra_id','chandra'),
                                ('xmm','xmm_id','xmm'),
                            ]:
                                if (colname in combined_df.columns) and (not _is_missing(mrow.get(colname, None))):
                                    short = _short_of(kcat, mrow[colname])
                                    if short is not None:
                                        parts.append(f"{label} #{short}")
                            used_info = (" [master: " + ", ".join(parts) + "]") if parts else ""
                        except Exception:
                            used_info = ""
                    else:
                        used_info = ""
                    ax.text(0.5, 0.965, txt_master + " | " + txt_new + used_info,
                            transform=ax.transAxes, ha='center', va='top', fontsize=8, color='#444444', zorder=15,
                            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.75))
                except Exception:
                    pass

            # Optional: dump D² internals directly on the panel and to stdout
            if PDF_DEBUG_D2 and debug_info_lines:
                try:
                    y_dbg = 0.935
                    for line in debug_info_lines:
                        print(f"[D2] {line}")
                        ax.text(0.02, y_dbg, line, transform=ax.transAxes, ha='left', va='top',
                                fontsize=7, family='monospace', color='#333333',
                                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.65))
                        y_dbg -= 0.028
                except Exception:
                    pass

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
            legend_handles.append(
                Line2D([], [], color='#c0c0c0', marker='D', linestyle='None',
                       markerfacecolor='none', label='master pos')
            )
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
                    try:
                        pa_other = float(other_df.iloc[j]['errPA'])
                    except Exception:
                        pa_other = 0.0
                    pa_draw = pa_other - 90.0
                    try:
                        pa_draw = ((pa_draw + 540.0) % 360.0) - 180.0
                    except Exception:
                        pass
                    reg.append(
                        'ellipse(%.7f,%.7f,%.3f\",%.3f\",%.1f) # color=%s width=2' % (
                            ra0, dec0, oth_maj_plot[j], oth_min_plot[j],
                            pa_draw,
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
                        % (ra0, dec0, oth_maj_plot[j], oth_min_plot[j], pa_draw)
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
                                pa_g_draw = pa_g - 90.0
                                try:
                                    pa_g_draw = ((pa_g_draw + 540.0) % 360.0) - 180.0
                                except Exception:
                                    pass
                                ajs.append('draw blue ellipse(%.7f %.7f %.3farcsec %.3farcsec %.1f)'
                                           % (r_g, d_g, maj_g, min_g, pa_g_draw))
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
                    close_unmatched = 0
                    close_unmatched_any = 0
                    try:
                        if master_indices:
                            for _m in master_indices:
                                # Skip synthetic rows if present
                                is_synthetic = False
                                try:
                                    is_synthetic = ('gaia_id' in combined_df.columns) and (
                                        pd.isna(combined_df.iloc[_m].get('gaia_id', np.nan)) or
                                        combined_df.iloc[_m].get('gaia_id', -1) == -1
                                    )
                                except Exception:
                                    is_synthetic = False
                                if is_synthetic:
                                    continue
                                d_arcsec = float(np.hypot(*_src_offset(_m)))
                                if d_arcsec <= near_thr:
                                    close_count += 1
                                    try:
                                        id_val = combined_df.iloc[_m].get(id_col, None) if id_col in combined_df.columns else None
                                        if _is_missing(id_val):
                                            close_unmatched += 1
                                    except Exception:
                                        pass
                    except Exception:
                        pass
                    try:
                        for _cand in cand_indices:
                            try:
                                dx_c, dy_c = _src_offset(_cand)
                            except Exception:
                                continue
                            if np.hypot(dx_c, dy_c) <= near_thr:
                                try:
                                    if str(combined_df.iloc[_cand].get(id_col, '')) == str(new_id):
                                        continue
                                except Exception:
                                    pass
                                try:
                                    id_val = combined_df.iloc[_cand].get(id_col, None) if id_col in combined_df.columns else None
                                    if _is_missing(id_val):
                                        close_unmatched_any += 1
                                except Exception:
                                    continue
                    except Exception:
                        pass
                    is_problem = False
                    # Problematic if: (no identifications and a nearby master) or (an identification exists and another close-by master is present)
                    if (num_valid_master == 0 and (near_sep is not None) and (near_sep <= near_thr)):
                        is_problem = True
                        logger.debug(f"[{base}] Problem detected: Case A (no identifications, nearby master) - "
                                   f"num_valid_master=0, near_sep={near_sep:.3f}, near_thr={near_thr:.3f}")
                    elif (num_valid_master >= 1 and close_count >= 2):
                        is_problem = True
                        logger.debug(f"[{base}] Problem detected: Case B (multiple nearby masters) - "
                                   f"num_valid_master={num_valid_master}, close_count={close_count}")
                    elif (num_valid_master >= 1 and close_unmatched >= 1):
                        is_problem = True
                        logger.debug(f"[{base}] Problem detected: Case D (close unmatched candidate) - "
                                   f"num_valid_master={num_valid_master}, close_unmatched={close_unmatched}")
                    elif close_unmatched_any >= 1:
                        is_problem = True
                        logger.debug(f"[{base}] Problem detected: Case E (nearby unmatched candidate) - "
                                   f"close_unmatched_any={close_unmatched_any}")
                    # Also mark as problematic when a selected match has D^2 > 1
                    try:
                        if (not is_problem) and (best_d2 is not None) and np.isfinite(best_d2) and (best_d2 > 1.0):
                            is_problem = True
                            logger.debug(f"[{base}] Problem detected: Case C (D² > 1) - best_d2={best_d2:.3f}")
                    except Exception:
                        pass

                    # Log all panels for debugging
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug(f"[{base}] Panel j={j}, new_id={new_id}, is_problem={is_problem}, "
                                   f"num_valid_master={num_valid_master}, near_sep={near_sep}, "
                                   f"near_thr={near_thr:.3f}, close_count={close_count}, best_d2={best_d2}")

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
                nav_js = NAV_JS_TEMPLATE
                # Always (over)write nav_js so updated bounds logic applies to all pages
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

            # Log page-level problem detection
            if logger.isEnabledFor(logging.DEBUG):
                problem_panels = [a.get('base', '?') for a in page_areas if a.get('problem')]
                logger.debug(f"Page {page_num}: {len(page_areas)} panels, {len(problem_panels)} problematic, "
                           f"page_has_problem={page_has_problem}")
                if problem_panels:
                    logger.debug(f"  Problematic panels on page {page_num}: {', '.join(problem_panels)}")

            if page_has_problem:
                if page_num not in problem_page_set:
                    pages_marked_problem.add(page_num)
                problem_page_set.add(page_num)
                logger.info(f"Page {page_num} marked problematic (total problem pages now {len(problem_page_set)})")
            else:
                if page_num in problem_page_set:
                    pages_cleared_problem.add(page_num)
                problem_page_set.discard(page_num)
            skip_html = (_html_mode == 'none') or (_html_mode in ('problems', 'problem', 'problem-only') and not page_has_problem)

            profiler.tic('page_output')
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
            for area in page_areas:
                # Compute overlays from axes window extents with fixed insets
                ax_index = area.get('ax_index')
                try:
                    ax = axes_flat[ax_index]
                except Exception:
                    continue
                # Use axes position in normalized figure coordinates for robust, DPI‑independent mapping
                try:
                    pos = ax.get_position()
                except Exception:
                    pos = None
                if pos is None:
                    # Fallback: treat as full-figure panel
                    x0f = 0.0; y0f = 0.0; x1f = 1.0; y1f = 1.0
                else:
                    x0f, y0f, x1f, y1f = float(pos.x0), float(pos.y0), float(pos.x1), float(pos.y1)
                # Convert to pixel coordinates of the displayed image area
                ax_l = int(round(x0f * disp_w))
                ax_r = int(round(x1f * disp_w))
                ax_top = int(round((1.0 - y1f) * disp_h))
                ax_bot = int(round((1.0 - y0f) * disp_h))
                ax_w = max(0, ax_r - ax_l)
                ax_h = max(0, ax_bot - ax_top)
                # Try to align HTML overlays with the actual text bboxes from the figure
                # for 'Aladin' (left) and 'Aladin Lite' (right). If not found, fall back to
                # transAxes anchors near the corners.
                try:
                    renderer = fig.canvas.get_renderer()
                except Exception:
                    renderer = None

                def _find_text_bbox(ax_obj, label: str):
                    if renderer is None:
                        return None
                    cand = None
                    for t in getattr(ax_obj, 'texts', []):
                        try:
                            if str(t.get_text()).strip().lower().startswith(label):
                                cand = t
                                break
                        except Exception:
                            continue
                    if cand is None:
                        return None
                    try:
                        vis = cand.get_visible()
                        if not vis:
                            cand.set_visible(True)
                        bb = cand.get_window_extent(renderer=renderer)
                        if not vis:
                            cand.set_visible(False)
                        return bb
                    except Exception:
                        return None

                # Button target sizes (fallbacks); enlarge slightly and allow env overrides
                try:
                    _minw_env = int(str(os.environ.get('HTML_BTN_MIN_W', '160')))
                except Exception:
                    _minw_env = 160
                try:
                    _minh_env = int(str(os.environ.get('HTML_BTN_MIN_H', '36')))
                except Exception:
                    _minh_env = 36
                try:
                    _pad_env = int(str(os.environ.get('HTML_BTN_PAD', '8')))
                except Exception:
                    _pad_env = 8
                BTN_W_FALLBACK = int(round(_minw_env * HTML_SCALE))
                BTN_H_FALLBACK = int(round(_minh_env * HTML_SCALE))
                pad_px = int(round(_pad_env * HTML_SCALE))

                # Helpers to find exact vs prefix matches
                def _find_text_bbox_exact(ax_obj, label_exact: str):
                    if renderer is None:
                        return None
                    label_exact = label_exact.strip().lower()
                    for t in getattr(ax_obj, 'texts', []):
                        try:
                            if str(t.get_text()).strip().lower() == label_exact:
                                vis = t.get_visible()
                                if not vis:
                                    t.set_visible(True)
                                bb = t.get_window_extent(renderer=renderer)
                                if not vis:
                                    t.set_visible(False)
                                return bb
                        except Exception:
                            continue
                    return None

                # Left: 'Aladin' (prefer exact match to avoid catching 'Aladin Lite')
                bb_left = _find_text_bbox_exact(ax, 'aladin')
                if bb_left is None:
                    # Try a stricter prefix with trailing space to avoid 'lite'
                    bb_left = _find_text_bbox(ax, 'aladin ')
                if bb_left is not None:
                    lx1 = int(round(bb_left.x0 * xscale)) - pad_px
                    ly1 = int(round((fig_h - bb_left.y1) * yscale)) - pad_px
                    lw = int(round(bb_left.width * xscale)) + 2 * pad_px
                    lh = int(round(bb_left.height * yscale)) + 2 * pad_px
                    # Clamp into the axes rectangle
                    left_x1 = max(ax_l, min(lx1, ax_r))
                    left_y1 = max(ax_top, min(ly1, ax_bot))
                    left_x2 = min(ax_r, left_x1 + lw)
                    left_y2 = min(ax_bot, left_y1 + lh)
                    # Ensure minimum size
                    if (left_x2 - left_x1) < BTN_W_FALLBACK:
                        left_x2 = min(ax_r, left_x1 + BTN_W_FALLBACK)
                    if (left_y2 - left_y1) < BTN_H_FALLBACK:
                        left_y2 = min(ax_bot, left_y1 + BTN_H_FALLBACK)
                else:
                    # Fallback to corner anchor
                    w_frac = max(0.0, x1f - x0f)
                    h_frac = max(0.0, y1f - y0f)
                    lx_frac = x0f + 0.02 * w_frac
                    ly_frac = y1f - 0.02 * h_frac
                    left_x1 = int(round(lx_frac * disp_w))
                    left_y1 = int(round((1.0 - ly_frac) * disp_h))
                    left_x1 = max(ax_l, min(left_x1, max(ax_l, ax_r - BTN_W_FALLBACK)))
                    left_y1 = max(ax_top, min(left_y1, max(ax_top, ax_bot - BTN_H_FALLBACK)))
                    left_x2 = min(ax_r, left_x1 + BTN_W_FALLBACK)
                    left_y2 = min(ax_bot, left_y1 + BTN_H_FALLBACK)

                # Right: 'Aladin Lite'
                bb_right = _find_text_bbox(ax, 'aladin lite')
                if bb_right is not None:
                    rx2 = int(round(bb_right.x1 * xscale)) + pad_px
                    ry1 = int(round((fig_h - bb_right.y1) * yscale)) - pad_px
                    rw = int(round(bb_right.width * xscale)) + 2 * pad_px
                    rh = int(round(bb_right.height * yscale)) + 2 * pad_px
                    right_x2 = min(ax_r, max(ax_l, rx2))
                    right_y1 = max(ax_top, min(ry1, ax_bot))
                    right_x1 = max(ax_l, right_x2 - rw)
                    right_y2 = min(ax_bot, right_y1 + rh)
                    # Ensure minimum size
                    if (right_x2 - right_x1) < BTN_W_FALLBACK:
                        right_x1 = max(ax_l, right_x2 - BTN_W_FALLBACK)
                    if (right_y2 - right_y1) < BTN_H_FALLBACK:
                        right_y2 = min(ax_bot, right_y1 + BTN_H_FALLBACK)
                else:
                    w_frac = max(0.0, x1f - x0f)
                    h_frac = max(0.0, y1f - y0f)
                    rx_frac = x0f + 0.98 * w_frac
                    ry_frac = y1f - 0.02 * h_frac
                    right_x2 = int(round(rx_frac * disp_w))
                    right_y1 = int(round((1.0 - ry_frac) * disp_h))
                    right_x2 = min(ax_r, max(ax_l + BTN_W_FALLBACK, right_x2))
                    right_y1 = max(ax_top, min(right_y1, max(ax_top, ax_bot - BTN_H_FALLBACK)))
                    right_x1 = max(ax_l, right_x2 - BTN_W_FALLBACK)
                    right_y2 = min(ax_bot, right_y1 + BTN_H_FALLBACK)
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
            debug_attr = ' checked' if debug_env else ''
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
                '.btn.tgl{background:#666;padding:1px 6px;font-weight:bold;line-height:1;border-radius:3px}'
                '.panelbox .delbtn{position:absolute;bottom:6px;left:50%;transform:translateX(-50%);padding:4px 10px;'
                'font-size:12px;line-height:1;border:none;border-radius:4px;background:#c0392b;color:#fff;'
                'cursor:pointer;z-index:13;opacity:.92}'
                '.panelbox .delbtn:hover{opacity:1}'
                '.nav{margin-bottom:8px;font:14px -apple-system,BlinkMacSystemFont,Segoe UI,Arial}'
                '.kbdbar{margin-top:12px;padding:10px 12px;background:#fff;border:1px solid rgba(0,0,0,.12);'
                'border-radius:6px;box-shadow:0 1px 3px rgba(0,0,0,.06);display:flex;flex-wrap:wrap;gap:12px;'
                'font:13px/1.5 -apple-system,BlinkMacSystemFont,Segoe UI,Arial;color:#333}'
                '.kbdbar span{display:flex;align-items:center;gap:6px}'
                '.kbdbar kbd{background:#f3f3f3;border:1px solid rgba(0,0,0,.18);border-radius:4px;padding:2px 6px;'
                'font-size:12px;font-family:"SFMono-Regular",Consolas,"Liberation Mono",monospace;'
                'box-shadow:inset 0 -1px 0 rgba(0,0,0,.08)}'
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
            # Filter toggle controls
            lines.append(' | <label><input type="checkbox" id="only_prob"> Show only to-check</label>')
            lines.append(f' | <label><input type="checkbox" id="show_debug"{debug_attr}> Show D² debug</label>')
            lines.append(' | <label><input type="checkbox" id="auto_aladin"> Auto Aladin</label>')
            # Top-right link to go back one level (master index). Default points to static index
            # two levels up; when served via dyn_server, a small script rewrites it to '/'.
            lines.append('<span style="float:right"><a id="back_link" href="../../aladin_index.html" title="Back to index">Back</a></span>')
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
                    # Use relative path so dynamic server can handle the request; JS will rewrite
                    # to the configured SAMP_URL when available.
                    default_http = f'/run_samp?file={base_ajs}'
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
                    other_attr = escape(str(new_id))
                    cat_attr = escape(str(key))
                    del_attr = '1' if is_deleted else '0'
                    lines.append(
                        f'<div class="panelbox" data-base="{base_id}" data-problem="{prob}" data-cat="{cat_attr}" data-other="{other_attr}" data-deleted="{del_attr}" '
                        f'style="position:absolute;left:{lp:.3f}%;top:{tp:.3f}%;width:{wp:.3f}%;height:{hp:.3f}%;z-index:8;">'
                        f'<button type="button" class="delbtn" onclick="return deleteSource(this);">Delete</button>'
                        f'<input type="checkbox" class="skipbox" data-base="{base_id}" title="Skip" '
                        f'style="position:absolute;left:50%;top:4px;transform:translateX(-50%);z-index:12;" onclick="return toggleSkip(this);" />'
                        f'<div class="mask" style="display:none;position:absolute;left:0;top:0;width:100%;height:100%;background:rgba(255,255,255,0.85);"></div>'
                        f'</div>'
                    )
                    # Add small clickable toggle links for each master row in this panel (top-left of panel)
                    try:
                        # Identify master rows for this current source
                        try:
                            mask_mid = (combined_df[id_col].astype(str) == str(new_id))
                        except Exception:
                            mask_mid = np.zeros(len(combined_df), dtype=bool)
                        m_idxs = [int(ix) for ix in list(combined_df.index[mask_mid])]
                        # Compute anchors regardless of current matches so buttons persist after detach
                        if True:
                            # Compute pixel anchors aligned with the table printed inside the panel
                            panel_left, panel_top, panel_right, panel_bottom = a['coords']
                            panel_w = max(1, panel_right - panel_left)
                            panel_h = max(1, panel_bottom - panel_top)
                            # Reproduce the table layout used in the figure to align the toggle button
                            base_font = 9.0
                            # Scale font with grid like in the figure
                            area_scale = (4.0 / float(max(1, ncols * nrows))) ** 0.5
                            font_size = float(int(round(np.clip(base_font * area_scale, 6, 12))))
                            line_h = 0.040 * (font_size / base_font)  # in axes (0–1) units
                            # Table header top anchor in axes coords, then first row is one line lower
                            y0 = 0.95 - 0.02 * max(0, max(nrows, ncols) - 2)
                            y0 = max(0.85, min(0.95, y0))
                            first_row_y = y0 - line_h
                            # Convert table anchor to pixel coordinates inside this panel
                            table_x_px = panel_left + 0.05 * panel_w
                            table_y_px = panel_top + (1.0 - first_row_y) * panel_h
                            # Place a compact toggle just to the right of the left axis (inside the plot),
                            # to avoid covering tick labels. Use ~2% of panel width from the left, plus a tiny pad.
                            # Move slightly further left to avoid overlapping tick marks
                            x_btn_px = panel_left + (0.01 * panel_w) + 2.0
                            # Vertical step per row in pixels
                            y_step_px = line_h * panel_h
                            # helper to choose a stable master key and label
                            def _mk_of_row(mrow):
                                try:
                                    if ('gaia_id' in combined_df.columns) and (not pd.isna(mrow.get('gaia_id', np.nan))) and (int(mrow.get('gaia_id', -1)) != -1):
                                        return f"gaia:{int(mrow.get('gaia_id'))}", f"gaia {int(mrow.get('gaia_id'))}"
                                except Exception:
                                    pass
                                try:
                                    if ('2MASS' in combined_df.columns) and (not pd.isna(mrow.get('2MASS', np.nan))):
                                        return f"2MASS:{str(mrow.get('2MASS'))}", f"2MASS {str(mrow.get('2MASS'))}"
                                except Exception:
                                    pass
                                try:
                                    if ('wise_id' in combined_df.columns) and (not pd.isna(mrow.get('wise_id', np.nan))):
                                        return f"wise:{str(mrow.get('wise_id'))}", f"wise {str(mrow.get('wise_id'))}"
                                except Exception:
                                    pass
                                try:
                                    if ('chandra_id' in combined_df.columns) and (not pd.isna(mrow.get('chandra_id', np.nan))):
                                        return f"chandra:{str(mrow.get('chandra_id'))}", f"chandra {str(mrow.get('chandra_id'))}"
                                except Exception:
                                    pass
                                try:
                                    if ('xmm_id' in combined_df.columns) and (not pd.isna(mrow.get('xmm_id', np.nan))):
                                        return f"xmm:{str(mrow.get('xmm_id'))}", f"xmm {str(mrow.get('xmm_id'))}"
                                except Exception:
                                    pass
                                return None, None
                            # Render up to N toggle links (matched masters first, then nearest others)
                            N = 8
                            shown = 0
                            seen_keys = set()
                            def _emit_for_index(mi):
                                nonlocal shown
                                if shown >= N:
                                    return
                                mrow = combined_df.iloc[mi]
                                # Skip synthetic self for current catalog
                                try:
                                    _same = (str(mrow.get(id_col, '')) == str(new_id))
                                except Exception:
                                    _same = False
                                try:
                                    _synth = ('gaia_id' in combined_df.columns) and (pd.isna(mrow.get('gaia_id', np.nan)) or mrow.get('gaia_id', -1) == -1)
                                except Exception:
                                    _synth = False
                                if _same and _synth:
                                    return
                                mkey, mlabel = _mk_of_row(mrow)
                                if not mkey or mkey in seen_keys:
                                    return
                                seen_keys.add(mkey)
                                from urllib.parse import quote
                                _id_enc = quote(str(new_id), safe='')
                                _mk_enc = quote(str(mkey), safe='')
                                href = f"/api/edit_link?cat={key}&id={_id_enc}&master={_mk_enc}&action=toggle&page={page_num}&row={mi}"
                                lx = x_btn_px
                                # Slight vertical offset so the toggle aligns with the text baseline
                                ly = table_y_px + shown * y_step_px + 2.0
                                btn_label = "±"
                                lines.append(
                                    f'<a class="btn tgl" style="left:{(lx/disp_w)*100:.3f}%;top:{(ly/disp_h)*100:.3f}%;width:auto;height:auto;transform:scale(0.9);transform-origin:left top;" '
                                    f'href="{href}" title="Toggle identification for {mlabel}">{btn_label}</a>'
                                )
                                shown += 1
                            # 1) Current master rows first
                            ordered_indices = toggle_order if toggle_order else []
                            used_match_idxs = set()
                            for mi in ordered_indices:
                                _emit_for_index(mi)
                                used_match_idxs.add(mi)
                                if shown >= N:
                                    break
                            if shown < N:
                                for mi in m_idxs:
                                    if mi in used_match_idxs:
                                        continue
                                    _emit_for_index(mi)
                                    if shown >= N:
                                        break
                            # 2) Then other master rows in panel proximity
                            if shown < N:
                                for mi, _dist in cand_with_dist:
                                    _emit_for_index(mi)
                                    if shown >= N:
                                        break
                    except Exception:
                        pass
            lines.append('</div>')
            lines.append(
                '<div class="kbdbar" role="contentinfo">'
                '<span><kbd>←</kbd><kbd>→</kbd> page</span>'
                '<span><kbd>a</kbd> Aladin</span>'
                '<span><kbd>g</kbd> goto</span>'
                '<span><kbd>1</kbd><kbd>2</kbd><kbd>3</kbd> toggle</span>'
                '</div>'
            )
            # JS helper to call local SAMP link server when available
            lines += [
                '<script>\n'
                f'const ALADIN_ITEMS = {aladin_items_js!s};\n'
                'var __orig = (window.location && window.location.origin) || "";\n'
                'var __def = (__orig && /^https?:/i.test(__orig)) ? __orig : "http://127.0.0.1:8765";\n'
                'const SAMP_URL = (window.localStorage && localStorage.getItem("ALADIN_LINK_URL")) || __def;\n'
                
                f'const PAGE_KEY = {key!r}; const PAGE_NUM = {page_num}; const TOTAL_PAGES = {total_pages};\n'
                'try{ window.PAGE_KEY = PAGE_KEY; window.PAGE_NUM = PAGE_NUM; window.TOTAL_PAGES = TOTAL_PAGES; }catch(e){}\n'
                '// Extract actual PAGE_KEY from current URL (handles radius suffix like _r0.1)\n'
                '(function(){ try{ var m = String(location.pathname || "").match(/([^\\/]+)_page\\d+\\.html$/); if(m && m[1]){ window.PAGE_KEY = m[1]; } }catch(e){} })();\n'
                '// Fix Back link when served via dyn_server (HTTP), leave static default otherwise\n'
                '(function(){ try{ var a=document.getElementById("back_link"); if(!a) return; var isHttp=/^https?:/i.test(location.protocol); var p=location.pathname||""; if(isHttp && p.indexOf("/aladin_scripts/html/")!==-1){ a.setAttribute("href","/"); } }catch(e){} })();\n'
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
                'function deleteSource(btn){ try{ const panel=btn.closest(".panelbox"); if(!panel) return false; const cat=panel.getAttribute("data-cat"); const other=panel.getAttribute("data-other"); if(!cat||!other) return false; if(!(window.confirm("Remove this source from the master catalog?"))){ return false; } btn.disabled=true; fetch("/api/delete_source",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({cat:cat,id:other,page:PAGE_NUM})}).then(function(resp){ if(!resp.ok){ throw new Error(resp.status); } return resp.json().catch(function(){ return {ok:true}; }); }).then(function(payload){ if(payload && payload.ok){ location.reload(); } else { btn.disabled=false; alert((payload && payload.error) || "Delete failed."); } }).catch(function(err){ btn.disabled=false; alert("Delete failed: "+err); }); }catch(e){ alert("Delete failed."); } return false; }\n'
                'function _readSkip(){ try{ const k = "SKIP_"+PAGE_KEY; const s = localStorage.getItem(k); return s? JSON.parse(s) : {}; }catch(e){ return {}; } }\n'
                'function _writeSkip(obj){ try{ const k = "SKIP_"+PAGE_KEY; localStorage.setItem(k, JSON.stringify(obj)); }catch(e){} }\n'
                'function _readOnly(){ try{ const k = "ONLY_PROB_"+PAGE_KEY; return localStorage.getItem(k)==='+'"1"'+'; }catch(e){ return false; } }\n'
                'function _writeOnly(v){ try{ const k = "ONLY_PROB_"+PAGE_KEY; localStorage.setItem(k, v?"1":"0"); }catch(e){} }\n'
                'function toggleSkip(cb){ try{ const base = cb.getAttribute("data-base"); if(!base) return false; const sk=_readSkip(); if(cb.checked){ sk[base]=1; } else { delete sk[base]; } _writeSkip(sk); const pb = cb.closest(".panelbox"); if(pb){ pb.setAttribute("data-skip", cb.checked?"1":"0"); } updatePanels(); }catch(e){} return true; }\n'
                'function updatePanels(){ try{ const only = document.getElementById("only_prob"); const onlyOn = !!(only && only.checked); const pb = document.querySelectorAll(".panelbox"); pb.forEach(function(p){ const prob = p.getAttribute("data-problem")==="1"; const skip = p.getAttribute("data-skip")==="1"; const mask = p.querySelector(".mask"); let show = true; if(onlyOn){ show = prob && !skip; } else { show = !skip; } if(mask){ mask.style.display = show? "none" : "block"; } }); }catch(e){} return false; }\n'
                'function _hasUnskippedProblems(){ try{ const pb = document.querySelectorAll(".panelbox"); for(let i=0;i<pb.length;i++){ const p=pb[i]; if(p.getAttribute("data-problem")==="1" && p.getAttribute("data-skip")!=="1"){ return true; } } }catch(e){} return false; }\n'
                'function _getSeek(){ try{ const h=String(location.hash||""); const m=h.match(/seek=(next|prev)/); return m? m[1] : "next"; }catch(e){ return "next"; } }\n'
                'function _navToPage(n, dir, opts){ if(typeof TOTAL_PAGES!=="number") return; if(n<1||n>TOTAL_PAGES) return; const key = (typeof window.PAGE_KEY==="string" && window.PAGE_KEY) ? window.PAGE_KEY : PAGE_KEY; try{ const auto = document.getElementById("auto_aladin"); const onlyBox = document.getElementById("only_prob"); if(auto && onlyBox && auto.checked && onlyBox.checked){ sessionStorage.setItem("AUTO_ALADIN_SUPPRESS","1"); } }catch(e){} let params; try{ params = new URLSearchParams(window.location.search||""); }catch(e){ params = null; } if(!params){ params = new URLSearchParams(); } let onlyOn = false; try{ const box = document.getElementById("only_prob"); onlyOn = !!(box && box.checked); }catch(e){} if(onlyOn){ params.set("only","to_check"); params.set("img","skip"); } else { params.delete("only"); if(params.get("img")==="skip"){ params.delete("img"); } params.delete("bg"); } if(opts && Object.prototype.hasOwnProperty.call(opts,"img")){ if(opts.img === null){ params.delete("img"); } else { params.set("img", String(opts.img)); } } if(opts && Object.prototype.hasOwnProperty.call(opts,"bg")){ if(opts.bg === null){ params.delete("bg"); } else { params.set("bg", String(opts.bg)); } } if(opts && Object.prototype.hasOwnProperty.call(opts,"debug")){ params.set("debug", opts.debug ? "1" : "0"); } const query = params.toString(); let href = key + "_page" + n + ".html"; if(query){ href += "?" + query; } href += "#seek=" + (dir||"next"); if(opts && opts.replace){ window.location.replace(href); } else { window.location.href = href; } }\n'
                'function _maybeAutoAdvance(){ try{ if(typeof maybeAutoAdvance === "function"){ return maybeAutoAdvance(window.PROBLEM_PAGES || []); } }catch(e){} return false; }\n'
                'document.addEventListener("DOMContentLoaded", function(){ try{ const sk=_readSkip(); document.querySelectorAll(".panelbox").forEach(function(p){ const b=p.getAttribute("data-base"); if(sk[b]){ p.setAttribute("data-skip","1"); const cb=p.querySelector(".skipbox"); if(cb) cb.checked=true; } }); const only=document.getElementById("only_prob"); if(only){ only.checked = _readOnly(); only.addEventListener("change", function(){ _writeOnly(only.checked); updatePanels(); if(only.checked){ var hopped = _maybeAutoAdvance(); if(!hopped){ maybeTriggerAuto(); } } }); } updatePanels(); // tag Prev/Next anchors with seek\n'
                '  document.querySelectorAll(".nav a").forEach(function(a){ const label=(a.textContent||"").toLowerCase(); if(label.indexOf("next page")!==-1){ if(a.href.indexOf("#seek=")===-1) a.href += "#seek=next"; } if(label.indexOf("prev page")!==-1 || label.indexOf("« prev page")!==-1){ if(a.href.indexOf("#seek=")===-1) a.href += "#seek=prev"; } a.addEventListener("click", function(ev){ try{ if(typeof _navToPage !== "function") return; const raw = a.getAttribute("href") || ""; const match = raw.match(/_page(\\d+)\\.html/); if(!match) return; ev.preventDefault(); const target = parseInt(match[1],10); const seek = (raw.indexOf("#seek=prev")!==-1 || label.indexOf("prev page")!==-1) ? "prev" : "next"; _navToPage(target, seek); }catch(err){} }); }); const dbg=document.getElementById("show_debug"); if(dbg){ dbg.addEventListener("change", function(){ try{ const on = !!dbg.checked; let target = "/page/" + PAGE_KEY + "/" + PAGE_NUM + "?debug=" + (on ? "1" : "0"); try{ if(window.location.search && window.location.search.indexOf("img=1") !== -1){ target += "&img=1"; } }catch(e){} const hash = window.location.hash || ""; if(hash && hash.indexOf("seek=") >= 0){ window.location.href = target + hash; } else { window.location.href = target; } }catch(e){ window.location.href = "/page/" + PAGE_KEY + "/" + PAGE_NUM + "?debug=" + (dbg.checked ? "1" : "0"); } }); } }catch(e){} });\n'
                '</script>'
            ]
            # Remove bottom per-source nav (was unreliable without local server)
            try:
                _prob_json = _json.dumps(sorted(problem_page_set))
            except Exception:
                _prob_json = '[]'
            lines.append(f'<script>window.PAGE_HAS_PROBLEM = {"true" if page_has_problem else "false"}; window.PROBLEM_PAGES = {_prob_json};</script>')
            # Load complete problem pages list from JSON index
            lines.append(f'''<script>
(function(){{
  try{{
    var idxKey = (typeof window.PAGE_KEY === 'string' && window.PAGE_KEY) ? window.PAGE_KEY : {key!r};
    if(!idxKey){{
      var m = String(window.location && window.location.pathname || '').match(/([^\\/]+)_page\\d+\\.html$/);
      if(m && m[1]){{ idxKey = m[1]; }}
    }}
    if(!idxKey){{ idxKey = {key!r}; }}
    var idxUrl = idxKey + '_index.json?_t=' + Date.now();
    fetch(idxUrl)
      .then(function(r){{ return r.json(); }})
      .then(function(data){{
        if(data && Array.isArray(data.problem_pages)){{
          window.PROBLEM_PAGES = data.problem_pages;
          console.log('[nav_keys] Loaded PROBLEM_PAGES from index:', window.PROBLEM_PAGES.length, 'pages');
        }}
      }})
      .catch(function(e){{ console.warn('[nav_keys] Could not load index:', e); }});
  }}catch(e){{ console.warn('[nav_keys] Index load error:', e); }}
}})();
</script>''')
            lines.append(f'<script src="nav_keys.js?v={NAV_KEYS_VERSION}"></script>')
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
        profiler.toc('page_output')
        profiler.toc('page_total')
    # Close SAMP connection
    if samp_client is not None:
        try:
            samp_client.disconnect()
        except Exception:
            pass
    if pdf:
        pdf.close()

    # Write a minimal per-catalog JSON index to accelerate navigation
    try:
        updated_list = sorted(problem_page_set)
        with open(idx_path, 'w') as f:
            _json.dump({'total_pages': int(total), 'problem_pages': updated_list}, f)
        print(f"[index] wrote {idx_path} with {len(updated_list)} problem pages")
        logger.info(f"Wrote index {idx_path}: total_pages={total}, problem_pages={len(updated_list)}")
        if pages_marked_problem:
            logger.info(f"  Newly marked pages this run: {sorted(pages_marked_problem)}")
        if pages_cleared_problem:
            logger.info(f"  Cleared problem flag for pages: {sorted(pages_cleared_problem)}")
    except Exception as _e:
        print(f"[index] could not write per-catalog JSON index: {_e}")

    if profiler.enabled:
        profiler.emit(f"pages={total_pages} panels={total}")


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
    parser.add_argument('--radius-deg', type=float, default=0.02,
                        help='Search radius around NGC 2264 center in degrees (default: 0.02)')
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
        'gaia':    {'factor': 2.5, 'min_radius': 0.10, 'epoch': 2016.0},
        '2MASS':   {'factor': 2.5, 'min_radius': 0.20, 'epoch': 2000.0},
        'wise':    {'factor': 2.0, 'min_radius': 1.00, 'epoch': 2010.5},
        'chandra': {'factor': 2.0, 'min_radius': 5.00, 'epoch': 2002.0},
        'xmm':     {'factor': 2.0, 'min_radius': 5.00, 'epoch': 2003.0}
    }
    # Define the central coordinate of NGC 2264 and search radius
    center = SkyCoord(ra=100.25 * u.deg, dec=9.883333 * u.deg, frame='icrs')
    # Use user-provided radius (degrees)
    try:
        radius_val = float(args.radius_deg)
    except Exception:
        radius_val = 1.00
    global CURRENT_EDITS_RADIUS
    CURRENT_EDITS_RADIUS = radius_val
    radius = radius_val * u.deg

    print("Querying Gaia DR3 ...")
    gaia = query_gaia(center, radius)
    print(f"Retrieved {len(gaia)} Gaia sources (raw)")
    gaia = filter_gaia_quality(gaia)

    print("Querying 2MASS ...")
    tmass = query_2mass(center, radius)
    print(f"Retrieved {len(tmass)} 2MASS sources (raw)")
    tmass = filter_2mass_quality(tmass)
    if tmass is not None and not tmass.empty:
        tmass = tmass.copy()
        default_factor = catalog_params['2MASS']['factor']
        default_min_radius = catalog_params['2MASS']['min_radius']
        tmass['match_factor'] = default_factor
        tmass['match_min_radius'] = default_min_radius
        q_series = pd.Series('', index=tmass.index, dtype=object)
        for col in ('Qflg', 'qflg', 'ph_qual'):
            if col in tmass.columns:
                q_series = tmass[col].astype(str).str.upper()
                break
        q_short = q_series.str.slice(0, 3)
        has_A = q_short.str.contains('A', na=False)
        bad_quality = ~has_A
        special_mask = bad_quality
        if special_mask.any():
            tmass.loc[special_mask, 'match_factor'] = 3.5
            tmass.loc[special_mask, 'match_min_radius'] = 0.5

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
        'frame': gaia.copy(),
        **catalog_params['gaia']
    }

    edits_all = _load_edits(radius_val)
    for df_other, id_col in [
        (tmass,     '2MASS'),
        (wise,      'wise_id'),
        (chan,      'chandra_id'),
        (xmm,       'xmm_id'),
    ]:
        # Strip '_id' suffix for catalog key
        key = id_col.replace('_id', '')
        params = catalog_params[key]
        edits_cat = edits_all.get(key, {}) if isinstance(edits_all, dict) else {}
        if not isinstance(edits_cat, dict):
            edits_cat = {}
        delete_ids: set[str] = set()
        for entry in (edits_cat.get('delete') or []):
            if entry is None:
                continue
            if isinstance(entry, dict):
                val = entry.get('id', entry.get('other'))
            else:
                val = entry
            if val is not None:
                delete_ids.add(str(val))

        def _row_match_params(row: pd.Series) -> tuple[float, float]:
            try:
                rf = float(row.get('match_factor', params['factor']))
            except Exception:
                rf = float(params['factor'])
            try:
                rm = float(row.get('match_min_radius', params['min_radius']))
            except Exception:
                rm = float(params['min_radius'])
            return rf, rm

        def _row_match_params(row: pd.Series) -> tuple[float, float]:
            try:
                rf = float(row.get('match_factor', params['factor']))
            except Exception:
                rf = float(params['factor'])
            try:
                rm = float(row.get('match_min_radius', params['min_radius']))
            except Exception:
                rm = float(params['min_radius'])
            return rf, rm

        def _row_match_params(row: pd.Series) -> tuple[float, float]:
            try:
                rf = float(row.get('match_factor', params['factor']))
            except Exception:
                rf = float(params['factor'])
            try:
                rm = float(row.get('match_min_radius', params['min_radius']))
            except Exception:
                rm = float(params['min_radius'])
            return rf, rm
        # Always include RA and Dec in each catalog
        cols: list[str] = ['ra_deg', 'dec_deg']
        for ellipse in ['errMaj', 'errMin', 'errPA']:
            if ellipse in df_other.columns:
                cols.append(ellipse)
        for extra_col in ('match_factor', 'match_min_radius'):
            if extra_col in df_other.columns and extra_col not in cols:
                cols.append(extra_col)
        for extra_col in ('match_factor', 'match_min_radius'):
            if extra_col in df_other.columns and extra_col not in cols:
                cols.append(extra_col)
        for extra_col in ('match_factor', 'match_min_radius'):
            if extra_col in df_other.columns and extra_col not in cols:
                cols.append(extra_col)
        for extra_col in ('match_factor', 'match_min_radius'):
            if extra_col in df_other.columns and extra_col not in cols:
                cols.append(extra_col)
        for extra_col in ('match_factor', 'match_min_radius'):
            if extra_col in df_other.columns and extra_col not in cols:
                cols.append(extra_col)
        catalogs[key] = {
            'data': df_other.set_index(id_col)[cols],
            'frame': df_other.copy(),
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
        # Apply user edits (remove/force/delete) for this catalog before merging
        matches = _apply_edits_to_matches(combined_df, matches, id_col, edits_cat)
        # Record which other ids have been matched in the main Gaia-anchored pass
        matched_other_ids = set()
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
        def _safe_lookup_other_row(key):
            trials = [key]
            try:
                trials.append(int(key))
            except Exception:
                pass
            s = str(key)
            trials.append(s)
            if id_col == '2MASS' and ' ' in s:
                trials.append(s.replace(' ', '+'))
                trials.append(s.replace(' ', '-'))
            for k in trials:
                try:
                    return df_other_indexed.loc[k]
                except KeyError:
                    continue
            return None
        # Helper: robust row lookup for cases where IDs in edits/matches
        # may contain minor typos (e.g., 2MASS id with a space instead of '+').
        def _safe_lookup_other_row(key):
            # Try native key, int(key), and str(key)
            trials = [key]
            try:
                trials.append(int(key))
            except Exception:
                pass
            s = str(key)
            trials.append(s)
            # Special-case normalization for 2MASS IDs
            if id_col == '2MASS':
                if ' ' in s:
                    trials.append(s.replace(' ', '+'))
                    trials.append(s.replace(' ', '-'))
            for k in trials:
                try:
                    return df_other_indexed.loc[k]
                except KeyError:
                    continue
            return None
        # First, update existing rows for matches
        for i in combined_df.index:
            matched_ids = matches[i]
            if not matched_ids:
                continue
            # Handle first match in place
            first_oid = matched_ids[0]
            if str(first_oid) in deleted_ids:
                continue
            matched_other_ids.add(first_oid)
            # Robust lookup of the matched other-row (handle mixed id types)
            other_row = _safe_lookup_other_row(first_oid)
            if other_row is None:
                # Skip this match if we cannot retrieve the source row
                continue
            # Scale other ellipse and enforce minimum radius
            row_factor, row_min_radius = _row_match_params(other_row)
            scaled_maj = max(other_row['errMaj'] * row_factor,
                             row_min_radius)
            scaled_min = max(other_row['errMin'] * row_factor,
                             row_min_radius)
            scaled_pa  = other_row['errPA']
            # Choose smaller AREA ellipse between current master and matched catalog
            try:
                area_master = float(combined_df.at[i, 'errMaj']) * float(combined_df.at[i, 'errMin'])
            except Exception:
                area_master = float('inf')
            area_other = float(scaled_maj) * float(scaled_min)
            if not np.isfinite(area_master) or (area_other < area_master):
                combined_df.at[i, 'ra_deg'] = float(other_row['ra_deg'])
                combined_df.at[i, 'dec_deg'] = float(other_row['dec_deg'])
                combined_df.at[i, 'errMaj']  = float(scaled_maj)
                combined_df.at[i, 'errMin']  = float(scaled_min)
                combined_df.at[i, 'errPA']   = float(scaled_pa)
            combined_df.at[i, id_col] = first_oid
            # Duplicate for any additional matches
            for oid in matched_ids[1:]:
                if str(oid) in deleted_ids:
                    continue
                matched_other_ids.add(oid)
                # Use robust lookup for additional matches as well
                other_row = _safe_lookup_other_row(oid)
                if other_row is None:
                    # Skip if not found in the catalog index
                    continue
                # Scale other ellipse and enforce minimum radius
                row_factor_extra, row_min_radius_extra = _row_match_params(other_row)
                scaled_maj = max(other_row['errMaj'] * row_factor_extra,
                                 row_min_radius_extra)
                scaled_min = max(other_row['errMin'] * row_factor_extra,
                                 row_min_radius_extra)
                scaled_pa  = other_row['errPA']
                # Determine which ellipse to use by area
                try:
                    area_master = float(combined_df.at[i, 'errMaj']) * float(combined_df.at[i, 'errMin'])
                except Exception:
                    area_master = float('inf')
                area_other = float(scaled_maj) * float(scaled_min)
                if (not np.isfinite(area_master)) or (area_other < area_master):
                    ra_sel, dec_sel = float(other_row['ra_deg']), float(other_row['dec_deg'])
                    maj_sel, min_sel, pa_sel = float(scaled_maj), float(scaled_min), float(scaled_pa)
                else:
                    ra_sel  = float(combined_df.at[i, 'ra_deg'])
                    dec_sel = float(combined_df.at[i, 'dec_deg'])
                    maj_sel = float(combined_df.at[i, 'errMaj'])
                    min_sel = float(combined_df.at[i, 'errMin'])
                    pa_sel  = float(combined_df.at[i, 'errPA'])
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
        # Chain 'other' ids to unmatched master rows (gaia_id == -1) using direct ellipse overlap
        combined_df = chain_other_to_unmatched_master(
            combined_df,
            df_other,
            id_col=id_col,
            factor=float(params['factor']),
            min_radius_arcsec=float(params['min_radius'])
        )

        # Then add sources in other_df with no matches (after chaining),
        # based on reverse matching against the updated combined_df
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
                if str(other_row[id_col]) in deleted_ids:
                    continue
                # Scale ellipse and enforce minimum radius
                row_factor_rev, row_min_radius_rev = _row_match_params(other_row)
                maj   = max(other_row['errMaj'] * row_factor_rev,
                            row_min_radius_rev)
                min_  = max(other_row['errMin'] * row_factor_rev,
                            row_min_radius_rev)
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
    combined_df.reset_index(drop=True, inplace=True)
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
                       include_catalogs=None,
                       radius_deg: float | None = None,
                       edits: Dict[str, Any] | None = None):
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
        'gaia':    {'factor': 2.5, 'min_radius': 0.10, 'epoch': 2016.0},
        '2MASS':   {'factor': 2.5, 'min_radius': 0.20, 'epoch': 2000.0},
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
    # Choose radius: explicit param > env DYN_RADIUS_DEG > default 1.00 deg
    try:
        _r_env = os.environ.get('DYN_RADIUS_DEG')
        _r_val = float(radius_deg if radius_deg is not None else (_r_env if _r_env is not None else 1.00))
    except Exception:
        _r_val = 1.00
    radius = _r_val * u.deg
    global CURRENT_EDITS_RADIUS
    CURRENT_EDITS_RADIUS = _r_val

    gaia = query_gaia(center, radius)
    gaia = filter_gaia_quality(gaia)
    tmass = wise = chan = xmm = None
    if (include_set is None) or ('2mass' in include_set or '2mass_j' in include_set or '2mass-j' in include_set or '2MASS' in (include_catalogs if include_catalogs else [])):
        tmass = query_2mass(center, radius)
        tmass = filter_2mass_quality(tmass)
        if tmass is not None and not tmass.empty:
            tmass = tmass.copy()
            default_factor = catalog_params['2MASS']['factor']
            default_min_radius = catalog_params['2MASS']['min_radius']
            tmass['match_factor'] = default_factor
            tmass['match_min_radius'] = default_min_radius
            q_series = pd.Series('', index=tmass.index, dtype=object)
            for col in ('Qflg', 'qflg', 'ph_qual'):
                if col in tmass.columns:
                    q_series = tmass[col].astype(str).str.upper()
                    break
            q_short = q_series.str.slice(0, 3)
            has_A = q_short.str.contains('A', na=False)
            bad_quality = ~has_A
            special_mask = bad_quality
            if special_mask.any():
                tmass.loc[special_mask, 'match_factor'] = 3.5
                tmass.loc[special_mask, 'match_min_radius'] = 0.5
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
        edits_cat = (edits.get(key) if isinstance(edits, dict) else {}) if edits is not None else {}
        if not isinstance(edits_cat, dict):
            edits_cat = {}
        delete_ids: set[str] = set()
        for entry in (edits_cat.get('delete') or []):
            if entry is None:
                continue
            if isinstance(entry, dict):
                val = entry.get('id', entry.get('other'))
            else:
                val = entry
            if val is not None:
                delete_ids.add(str(val))
        def _row_match_params(row: pd.Series) -> tuple[float, float]:
            try:
                rf = float(row.get('match_factor', params['factor']))
            except Exception:
                rf = float(params['factor'])
            try:
                rm = float(row.get('match_min_radius', params['min_radius']))
            except Exception:
                rm = float(params['min_radius'])
            return rf, rm
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
        # Apply user edits (remove/force/delete) for this catalog before merging
        matches = _apply_edits_to_matches(combined_df, matches, id_col, edits_cat)
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
        def _safe_lookup_other_row(key):
            trials = [key]
            try:
                trials.append(int(key))
            except Exception:
                pass
            s = str(key)
            trials.append(s)
            if id_col == '2MASS' and ' ' in s:
                trials.append(s.replace(' ', '+'))
                trials.append(s.replace(' ', '-'))
            for k in trials:
                try:
                    return df_other_indexed.loc[k]
                except KeyError:
                    continue
            return None
        for i in combined_df.index:
            matched_ids = matches[i]
            if not matched_ids:
                continue
            first_oid = matched_ids[0]
            if str(first_oid) in delete_ids:
                continue
            other_row = _safe_lookup_other_row(first_oid)
            if other_row is None:
                print(f"[web] warn: couldn't resolve {id_col} id {first_oid!r} in catalog index; skipping")
                continue
            row_factor, row_min_radius = _row_match_params(other_row)
            scaled_maj = max(other_row['errMaj'] * row_factor, row_min_radius)
            scaled_min = max(other_row['errMin'] * row_factor, row_min_radius)
            scaled_pa  = other_row['errPA']
            if scaled_maj < combined_df.at[i, 'errMaj']:
                combined_df.at[i, 'ra_deg'] = other_row['ra_deg']
                combined_df.at[i, 'dec_deg'] = other_row['dec_deg']
                combined_df.at[i, 'errMaj']  = scaled_maj
                combined_df.at[i, 'errMin']  = scaled_min
                combined_df.at[i, 'errPA']   = scaled_pa
            combined_df.at[i, id_col] = first_oid
            for oid in matched_ids[1:]:
                if str(oid) in delete_ids:
                    continue
                try:
                    lookup_key = int(oid)
                except ValueError:
                    lookup_key = oid
                other_row = _safe_lookup_other_row(lookup_key)
                if other_row is None:
                    print(f"[web] warn: couldn't resolve {id_col} id {lookup_key!r} in catalog index; skipping extra match")
                    continue
                row_factor_extra, row_min_radius_extra = _row_match_params(other_row)
                scaled_maj = max(other_row['errMaj'] * row_factor_extra, row_min_radius_extra)
                scaled_min = max(other_row['errMin'] * row_factor_extra, row_min_radius_extra)
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
        # Chain 'other' ids to unmatched master rows (gaia_id == -1) first, then
        # append synthetic rows for any remaining unmatched others
        combined_df = chain_other_to_unmatched_master(
            combined_df,
            df_other,
            id_col=id_col,
            factor=float(params['factor']),
            min_radius_arcsec=float(params['min_radius'])
        )
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
                if str(other_row[id_col]) in delete_ids:
                    continue
                row_factor_rev, row_min_radius_rev = _row_match_params(other_row)
                maj = max(other_row['errMaj'] * row_factor_rev, row_min_radius_rev)
                min_ = max(other_row['errMin'] * row_factor_rev, row_min_radius_rev)
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

    try:
        baseline_df = combined_df.copy(deep=True)
        combined_df = apply_link_edits_to_combined(baseline_df, catalogs, edits)
    except Exception:
        # Fall back to the pre-edit catalogue if anything goes wrong
        pass

    combined_df.reset_index(drop=True, inplace=True)
    return combined_df, catalogs


if __name__ == '__main__':
    combined_df, catalogs = main()
