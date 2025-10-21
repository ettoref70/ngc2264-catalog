from __future__ import annotations

import base64
import logging
import os
from pathlib import Path
from typing import Dict

from flask import Flask, abort, redirect, render_template_string, request, send_from_directory, url_for
import sys
import math
import threading
from concurrent.futures import ThreadPoolExecutor

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


BASE_DIR = Path(__file__).resolve().parent.parent
# Ensure project root is on sys.path so we can import local modules when launched as a submodule
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

# Local modules (after adjusting sys.path)
import cross_match_ngc2264 as cm  # type: ignore
from aladin_link_server import run_samp_script  # type: ignore
SCRIPTS_DIR = (BASE_DIR / "aladin_scripts").resolve()
HTML_DIR = (SCRIPTS_DIR / "html").resolve()
EDITS_DIR = (BASE_DIR / "edits").resolve()


def _is_safe_path(root: Path, target: Path) -> bool:
    """Check if target path is safely within root directory.

    Guards against path traversal attacks and symlink escapes.
    """
    try:
        root = root.resolve(strict=True)
        try:
            target = target.resolve(strict=True)
        except FileNotFoundError:
            # For non-existent files, validate parent exists and is safe
            try:
                resolved_parent = target.parent.resolve(strict=True)
                target = resolved_parent / target.name
            except Exception:
                return False
        # Python 3.9+: use is_relative_to if available, fallback to string check
        try:
            return target.is_relative_to(root)
        except AttributeError:
            return str(target).startswith(str(root) + os.sep)
    except Exception:
        return False


def _load_edits_for_radius(radius_deg: float | None):
    import json as _json
    path = Path(cm._edits_path(radius_deg))
    fallback = None
    if radius_deg is not None:
        fallback_candidate = Path(cm._edits_path(None))
        if fallback_candidate != path:
            fallback = fallback_candidate
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.warning(f"Could not create edits directory: {e}")
    data = {}
    try:
        if path.exists():
            with open(path, 'r') as f:
                try:
                    data = _json.load(f)
                except _json.JSONDecodeError as e:
                    logger.warning(f"Invalid JSON in {path}: {e}")
                    data = {}
        elif fallback and fallback.exists():
            with open(fallback, 'r') as f:
                try:
                    data = _json.load(f)
                except _json.JSONDecodeError as e:
                    logger.warning(f"Invalid JSON in {fallback}: {e}")
                    data = {}
    except Exception as e:
        logger.error(f"Error loading edits: {e}")
        data = {}
    return path, data


def _set_debug_mode_from_request(request) -> None:
    """Helper to set debug mode from request parameters."""
    debug_arg = request.args.get('debug')
    if debug_arg is not None:
        val = '1' if debug_arg.lower() in ('1', 'true', 'yes', 'on') else '0'
        os.environ['PDF_DEBUG_D2'] = val
        os.environ['ALADIN_SHOW_DEBUG'] = val
    else:
        os.environ['PDF_DEBUG_D2'] = '0'
        os.environ['ALADIN_SHOW_DEBUG'] = '0'


def create_app() -> Flask:
    app = Flask(__name__)

    try:
        _prefetch_workers = max(1, int(os.environ.get('DYN_PREFETCH_WORKERS', '2')))
    except Exception:
        _prefetch_workers = 2
    PREFETCH_EXEC = ThreadPoolExecutor(max_workers=_prefetch_workers)
    PREFETCH_LOCK = threading.Lock()
    PREFETCH_TASKS: set[tuple[str, int]] = set()

    # Cache datasets per requested catalog key (e.g., '2MASS', 'wise')
    # Each entry keeps the pre-edit baseline so we can reapply edits quickly.
    DATA: Dict[str, Dict[str, object]] = {}
    DATA_LOCK = threading.RLock()  # Reentrant lock for nested access

    def _try_build_from_combined(key: str):
        import pandas as _pd
        csv_path = (BASE_DIR / 'ngc2264_combined.csv')
        if not csv_path.exists():
            return None
        try:
            df = _pd.read_csv(csv_path)
        except Exception:
            return None
        # Minimal params
        params = {
            'gaia':    {'factor': 2.0, 'min_radius': 0.05, 'epoch': 2016.0},
            '2MASS':   {'factor': 2.0, 'min_radius': 0.10, 'epoch': 2000.0},
            'wise':    {'factor': 2.0, 'min_radius': 1.00, 'epoch': 2010.5},
            'chandra': {'factor': 2.0, 'min_radius': 5.00, 'epoch': 2002.0},
            'xmm':     {'factor': 2.0, 'min_radius': 5.00, 'epoch': 2003.0},
        }
        # Build catalogs subset: always include gaia + requested key (if available)
        catalogs: Dict[str, dict] = {}
        # Gaia data
        if 'gaia_id' in df.columns:
            gcols = [c for c in ['ra_deg','dec_deg','errMaj','errMin','errPA','pmra','pmdec','pmra_error','pmdec_error','ref_epoch'] if c in df.columns]
            gdf = df[~df['gaia_id'].isna() & (df['gaia_id'] != -1)].drop_duplicates('gaia_id')
            if not gdf.empty:
                catalogs['gaia'] = {'data': gdf.set_index('gaia_id')[gcols], 'frame': gdf[gcols + ['gaia_id']], **params['gaia']}
        # Other catalog
        key_norm = key.lower()
        id_map = {'2mass': '2MASS', 'wise': 'wise_id', 'chandra': 'chandra_id', 'xmm': 'xmm_id'}
        id_col = id_map.get(key_norm, key if key in df.columns else None)
        if id_col is None or id_col not in df.columns:
            return None
        ocols = [c for c in ['ra_deg','dec_deg','errMaj','errMin','errPA', id_col] if c in df.columns]
        odf = df[~df[id_col].isna() & (df[id_col] != -1)].drop_duplicates(id_col)
        if odf.empty:
            return None
        catalogs[key if key != 'wise_id' else 'wise'] = {
            'data': odf.set_index(id_col)[[c for c in ['ra_deg','dec_deg','errMaj','errMin','errPA'] if c in odf.columns]],
            'frame': odf[ocols],
            **params.get(key if key != 'wise_id' else 'wise', params['2MASS'])
        }
        return df, catalogs

    def _load_data_for_key(key: str):
        with DATA_LOCK:
            refresh_flag = (os.environ.get('DYN_REFRESH', '0') in ('1', 'true', 'yes'))
            entry = DATA.get(key)
            _r_env = os.environ.get('DYN_RADIUS_DEG')
            try:
                _r_val = float(_r_env) if _r_env is not None else None
            except Exception:
                _r_val = None
            same_radius = False
            if isinstance(entry, dict):
                stored_radius = entry.get('radius_deg')
                try:
                    if stored_radius is None and _r_val is None:
                        same_radius = True
                    elif (stored_radius is not None) and (_r_val is not None):
                        same_radius = abs(float(stored_radius) - float(_r_val)) <= 1e-9
                except Exception:
                    same_radius = False
                if (not refresh_flag) and same_radius:
                    if 'radius_suffix' not in entry:
                        entry['radius_suffix'] = cm._radius_suffix(entry.get('radius_deg'))
                    return entry['combined'], entry['catalogs']

        # Build from queries/caches to ensure per-catalog ellipses and positions are correct.
        refresh = refresh_flag
        # By default, dynamic pages should mirror the PDF for that catalog:
        # include only Gaia + previous catalogs up to and including `key` in the
        # merge order used by main(): [2MASS, wise, chandra, xmm].
        # You can override with DYN_INCLUDE_CATALOGS=gaia,2MASS,... if desired.
        inc_env = os.environ.get('DYN_INCLUDE_CATALOGS')
        if inc_env:
            include = [s.strip() for s in inc_env.split(',') if s.strip()]
        else:
            order = ['gaia', '2MASS', 'wise', 'chandra', 'xmm']
            k_norm = key.strip()
            # Normalize some aliases
            if k_norm.lower() in ('2mass', '2mass_j', '2massj'):
                k_norm = '2MASS'
            # Find position; if not found, default to Gaia only
            try:
                pos = order.index(k_norm)
            except ValueError:
                pos = 0
            include = order[:pos+1]
        # Load edit links and pass to builder
        edits_path, edits = _load_edits_for_radius(_r_val)

        # Get Chandra CSV path from environment (required)
        chandra_csv = os.environ.get('CHANDRA_CSV_PATH')
        if not chandra_csv:
            logger.warning("CHANDRA_CSV_PATH not set, Chandra catalog will be unavailable")

        combined_df, catalogs = cm.build_data_for_web(refresh=refresh,
                                                      chandra_csv_path=chandra_csv,
                                                      include_catalogs=include,
                                                      radius_deg=_r_val,
                                                      edits=edits)
        baseline_df = combined_df.attrs.get('_baseline_unforced')
        if baseline_df is None:
            try:
                baseline_df = combined_df.copy(deep=True)
            except Exception:
                baseline_df = combined_df
        radius_suffix = cm._radius_suffix(_r_val)

        with DATA_LOCK:
            DATA[key] = {
                'combined': combined_df,
                'baseline': baseline_df,
                'catalogs': catalogs,
                'radius_deg': _r_val,
                'edits_path': edits_path,
                'radius_suffix': radius_suffix,
            }
        return combined_df, catalogs

    def _schedule_prefetch(key: str, current_page: int, total_pages: int, draw_images: bool) -> None:
        # Prefetching is disabled to keep page generation strictly on-demand.
        return

    def _ensure_pages_for(key: str, page_num: int, *, draw_images: bool | None = None, prefetch: bool = True) -> None:
        # Respect existing PDF_DEBUG_D2 setting; default to '0' (overlay off)
        if 'PDF_DEBUG_D2' not in os.environ:
            os.environ['PDF_DEBUG_D2'] = '0'
        combined_df, catalogs = _load_data_for_key(key)
        entry = DATA.get(key)
        radius_suffix = ''
        if isinstance(entry, dict):
            radius_suffix = entry.get('radius_suffix') or cm._radius_suffix(entry.get('radius_deg'))
            entry['radius_suffix'] = radius_suffix
        page_prefix = f"{key}{radius_suffix or ''}"
        if key not in catalogs:
            abort(404)
        # Normalize key for lookup when built from CSV
        key_lookup = key
        if key not in catalogs and key == 'wise' and 'wise_id' in catalogs:
            key_lookup = 'wise'
        if key_lookup not in catalogs:
            abort(404)
        meta = catalogs[key_lookup]
        df_other = meta.get('frame')  # original with id column
        if df_other is None:
            abort(404)
        id_col = '2MASS' if key == '2MASS' else (f"{key}_id" if not key.endswith('_id') and key not in ('gaia',) else key)
        if id_col == 'gaia':
            # We don't render pages for the master catalog
            abort(404)
        if id_col not in df_other.columns:
            # try common names
            name_map = {
                'wise': 'wise_id',
                'chandra': 'chandra_id',
                'xmm': 'xmm_id',
                '2MASS': '2MASS',
            }
            id_col = name_map.get(key, id_col)
        # Make sure output directory exists
        HTML_DIR.mkdir(parents=True, exist_ok=True)
        # Respect typical grid and invert settings (env tunable)
        grid = os.environ.get('DYN_GRID', '1x1')
        try:
            gc, gr = grid.lower().split('x')
            ncols, nrows = max(1, int(gc)), max(1, int(gr))
        except Exception:
            ncols, nrows = 2, 2
        invert = (os.environ.get('DYN_INVERT_CMAP', '1') in ('1', 'true', 'yes'))
        # Default to draw images unless explicitly disabled or overridden via arg
        if draw_images is None:
            draw_images = not (os.environ.get('DYN_NO_IMAGES', '0') in ('1', 'true', 'yes'))
        per_page = max(1, ncols * nrows)
        total_pages = max(1, math.ceil(len(df_other) / per_page))
        # Generate only the requested page for this catalog (writes to aladin_scripts/html)
        # Use env var understood by the plotter to limit generation to a single page
        os.environ['ALADIN_PAGE_ONLY'] = str(int(page_num))
        prev_suffix_env = os.environ.get('ALADIN_PAGE_SUFFIX')
        try:
            if radius_suffix:
                os.environ['ALADIN_PAGE_SUFFIX'] = radius_suffix
            else:
                os.environ.pop('ALADIN_PAGE_SUFFIX', None)
        except Exception:
            pass
        try:
            if isinstance(entry, dict):
                entry['last_draw_images'] = bool(draw_images)
            cm.plot_after_merge(
                combined_df,
                df_other,
                id_col,
                catalogs,
                pdf_path=None,
                plot_mode='match',
                ncols=ncols,
                nrows=nrows,
                invert_cmap=invert,
                draw_images=draw_images,
                aladin_dir=str(SCRIPTS_DIR),
                samp_enabled=False,
                samp_addr=os.environ.get('SAMP_ADDR', '127.0.0.1')
            )
            if prefetch:
                _schedule_prefetch(key, page_num, total_pages, draw_images)
            if radius_suffix:
                try:
                    import shutil
                    base_prefix = key
                    target_prefix = f"{key}{radius_suffix}"
                    file_pairs = []
                    base_html = HTML_DIR / f"{base_prefix}_page{page_num}.html"
                    target_html = HTML_DIR / f"{target_prefix}_page{page_num}.html"
                    if base_html.exists():
                        file_pairs.append((base_html, target_html))
                    for ext in ('jpg', 'png', 'webp'):
                        base_img = HTML_DIR / f"{base_prefix}_page{page_num}.{ext}"
                        target_img = HTML_DIR / f"{target_prefix}_page{page_num}.{ext}"
                        if base_img.exists():
                            file_pairs.append((base_img, target_img))
                    base_index = HTML_DIR / f"{base_prefix}_index.json"
                    target_index = HTML_DIR / f"{target_prefix}_index.json"
                    if base_index.exists():
                        file_pairs.append((base_index, target_index))
                    for src_path, dest_path in file_pairs:
                        try:
                            if src_path.resolve() == dest_path.resolve():
                                continue
                        except Exception:
                            pass
                        dest_path.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(src_path, dest_path)
                except Exception:
                    pass
        finally:
            try:
                if prev_suffix_env is None:
                    os.environ.pop('ALADIN_PAGE_SUFFIX', None)
                else:
                    os.environ['ALADIN_PAGE_SUFFIX'] = prev_suffix_env
            except Exception:
                pass
        # Clean up the limiter so subsequent calls can generate other pages if needed
        try:
            del os.environ['ALADIN_PAGE_ONLY']
        except Exception:
            pass

    @app.get('/')
    def index():
        # Lightweight index that does not trigger data loading.
        # Pages will be generated on-demand when you click a catalog.
        items = ['2MASS', 'wise', 'chandra', 'xmm']
        html = [
            '<!doctype html>',
            '<meta charset="utf-8">',
            '<title>Dynamic Pages</title>',
            '<style>body{font:14px -apple-system,BlinkMacSystemFont,Segoe UI,Arial} li{margin:6px 0}</style>',
            '<h2>Dynamic Aladin Pages</h2>',
            '<ul>'
        ]
        for key in items:
            html.append(f'<li><a href="{url_for("page", key=key, page=1)}">{key}</a></li>')
        html.append('</ul>')
        return render_template_string('\n'.join(html))

    @app.get('/page/<key>/<int:page>')
    def page(key: str, page: int):
        # Ensure that static HTML pages exist for this catalog; then serve the requested page
        want_img = request.args.get('img', '').lower() in ('1','true','yes','on')
        _set_debug_mode_from_request(request)
        _ensure_pages_for(key, page, draw_images=want_img or None)
        entry = DATA.get(key)
        radius_suffix = ''
        if isinstance(entry, dict):
            radius_suffix = entry.get('radius_suffix') or cm._radius_suffix(entry.get('radius_deg'))
        page_prefix = f"{key}{radius_suffix or ''}"
        fname = f"{page_prefix}_page{page}.html"
        fpath = HTML_DIR / fname
        if not fpath.exists():
            abort(404)
        # Redirect to the canonical static path so relative links resolve
        return redirect(url_for('serve_generated', subpath=fname), code=302)

    # ---- Support for "show only to-check" page navigation ----
    def _load_page_index(key: str) -> dict:
        """Load the per-catalog index JSON (handles optional radius suffix)."""
        import json as _json
        entry = DATA.get(key)
        radius_suffix = ''
        if isinstance(entry, dict):
            radius_suffix = entry.get('radius_suffix') or cm._radius_suffix(entry.get('radius_deg'))
        base = f"{key}{radius_suffix or ''}"
        idx_path = (HTML_DIR / f"{base}_index.json")
        try:
            with open(idx_path, 'r') as f:
                return _json.load(f) or {}
        except Exception:
            return {}

    def _index_total_pages(index: dict) -> int:
        for k in ('total_pages', 'pages_total', 'num_pages', 'count_pages'):
            try:
                v = int(index.get(k))
                if v > 0:
                    return v
            except Exception:
                pass
        # Fallback: infer from explicit list
        try:
            pages = index.get('pages')
            if isinstance(pages, list):
                # if list of dicts with 'page'
                nums = []
                for it in pages:
                    try:
                        nums.append(int(it.get('page')))
                    except Exception:
                        continue
                if nums:
                    return max(nums)
        except Exception:
            pass
        return 0

    def _to_check_set(index: dict) -> set[int]:
        """Return the set of pages flagged as to-check.
        Recognizes multiple schemas: keys 'to_check_pages', 'tocheck_pages', 'problematic_pages',
        or a 'pages' array with per-item boolean flags under 'to_check'/'tocheck'/'problematic'.
        """
        pages = set()
        # direct lists
        for k in ('to_check_pages', 'tocheck_pages', 'problematic_pages'):
            try:
                lst = index.get(k)
                if isinstance(lst, list):
                    for p in lst:
                        try:
                            pages.add(int(p))
                        except Exception:
                            pass
            except Exception:
                pass
        # structured list
        try:
            arr = index.get('pages')
            if isinstance(arr, list):
                for it in arr:
                    try:
                        p = int(it.get('page'))
                    except Exception:
                        continue
                    flag = bool(it.get('to_check') or it.get('tocheck') or it.get('problematic'))
                    if flag:
                        pages.add(p)
        except Exception:
            pass
        return pages

    @app.get('/page/<key>/<int:page>/jump')
    def page_jump(key: str, page: int):
        """Redirect to next/prev page; supports limiting to pages marked "to-check".

        Query params:
          dir=next|prev           (default: next)
          only=to_check|all       (default: all)
          img=1|0                 passthrough to /page
          debug=1|0               passthrough to /page
        """
        direction = (request.args.get('dir') or 'next').lower()
        only = (request.args.get('only') or 'all').lower()
        want_img = request.args.get('img', '')
        debug = request.args.get('debug', '')

        # Touch current page so that *_index.json exists/updates for this catalog & radius
        try:
            _ensure_pages_for(key, max(1, page), draw_images=(want_img.lower() in ('1','true','yes','on')) or None, prefetch=False)
        except Exception:
            pass

        index = _load_page_index(key)
        total = _index_total_pages(index)
        if total <= 0:
            # Fallback: move by ±1 within [1, ...]
            target = max(1, page + (-1 if direction == 'prev' else 1))
            return redirect(url_for('page', key=key, page=target, img=want_img, debug=debug), code=302)

        if only == 'to_check':
            tset = sorted(_to_check_set(index))
            if not tset:
                target = page  # nothing to check; stay
            else:
                if direction == 'prev':
                    prevs = [p for p in tset if p < page]
                    target = (prevs[-1] if prevs else tset[-1])
                else:
                    nexts = [p for p in tset if p > page]
                    target = (nexts[0] if nexts else tset[0])
        else:
            # all pages with wrap-around
            if direction == 'prev':
                target = (page - 1) if page > 1 else total
            else:
                target = (page + 1) if page < total else 1

        return redirect(url_for('page', key=key, page=target, img=want_img, debug=debug), code=302)

    # Serve generated assets under /aladin_scripts/html/...
    @app.get('/aladin_scripts/html/<path:subpath>')
    def serve_generated(subpath: str):
        # If a requested page doesn't exist yet or we want to always refresh pages,
        # generate it on the fly by parsing the pattern <key>_page<N>.html
        import re as _re
        m = _re.match(r"^([^/]+)_page(\d+)\.(html|jpg|png|webp)$", subpath)
        target_subpath = subpath
        fpath = (HTML_DIR / target_subpath)
        if m:
            prefix_token = m.group(1)
            page_num = int(m.group(2))
            extension = m.group(3)
            want_img = request.args.get('img', '').lower() in ('1','true','yes','on')
            _set_debug_mode_from_request(request)
            refresh_env = os.environ.get('DYN_REFRESH', '0').lower() in ('1','true','yes')
            refresh_qs = request.args.get('refresh', '').lower() in ('1','true','yes','on')

            base_key = prefix_token
            if base_key not in DATA and '_r' in prefix_token:
                base_candidate = prefix_token.split('_r', 1)[0]
                if base_candidate:
                    base_key = base_candidate
            if base_key not in DATA:
                try:
                    _load_data_for_key(base_key)
                except Exception:
                    pass
            entry = DATA.get(base_key)
            radius_suffix = ''
            if isinstance(entry, dict):
                radius_suffix = entry.get('radius_suffix') or cm._radius_suffix(entry.get('radius_deg'))
                entry['radius_suffix'] = radius_suffix
            expected_prefix = f"{base_key}{radius_suffix or ''}"
            target_name = f"{expected_prefix}_page{page_num}.{extension}"
            target_subpath = target_name
            fpath = HTML_DIR / target_subpath

            # Desired draw mode mirrors the page handler: explicit query overrides env.
            # Note: want_img already set above, debug mode already set by earlier call
            draw_expectation = want_img or None
            if draw_expectation is None:
                draw_expectation = not (os.environ.get('DYN_NO_IMAGES', '0').lower() in ('1','true','yes','on'))

            needs_refresh = (not fpath.exists()) or refresh_env or refresh_qs
            entry = DATA.get(base_key)
            if isinstance(entry, dict):
                last_mode = entry.get('last_draw_images')
                if last_mode is not None and bool(last_mode) != bool(draw_expectation):
                    needs_refresh = True
            if needs_refresh:
                try:
                    _ensure_pages_for(base_key, page_num, draw_images=want_img or None)
                except Exception:
                    pass
                fpath = HTML_DIR / target_subpath

            if target_subpath != subpath:
                if _is_safe_path(HTML_DIR, fpath) and fpath.exists():
                    return redirect(url_for('serve_generated', subpath=target_subpath), code=302)
                # fall through to final safety check

        if not _is_safe_path(HTML_DIR, fpath) or not fpath.exists():
            abort(404)
        return send_from_directory(HTML_DIR, target_subpath)

    # Also expose the entire aladin_scripts tree (e.g., .ajs and regions)
    @app.get('/aladin_scripts/<path:subpath>')
    def serve_scripts(subpath: str):
        fpath = (SCRIPTS_DIR / subpath)
        if not _is_safe_path(SCRIPTS_DIR, fpath) or not fpath.exists():
            abort(404)
        return send_from_directory(SCRIPTS_DIR, subpath)

    # Edit-link API: force/remove/toggle links between current catalog source and a master row
    @app.route('/api/edit_link', methods=['GET', 'POST'])
    def api_edit_link():
        import json as _json
        cat = (request.values.get('cat') or '').strip()
        other = (request.values.get('id') or '').strip()
        master = (request.values.get('master') or '').strip()
        action = (request.values.get('action') or 'toggle').strip().lower()
        row = (request.values.get('row') or '').strip()
        page = int(request.values.get('page') or '1')
        if not cat or not other or not master:
            abort(400)
        EDITS_DIR.mkdir(parents=True, exist_ok=True)
        canonical_remove_row = None
        canonical_force_row = None
        radius_hint = None
        cached = DATA.get(cat)
        if isinstance(cached, dict):
            radius_hint = cached.get('radius_deg')
        edits_path, data = _load_edits_for_radius(radius_hint)
        try:
            sect = data.get(cat) or {'force': [], 'remove': []}
            data[cat] = sect
            # Canonicalize known id formats from URLs that may collapse '+' into space
            if cat == '2MASS' and ' ' in other and '+' not in other and '-' not in other:
                other = other.replace(' ', '+')
            # Helper to check presence
            def _contains(lst, other_id, mkey):
                for it in lst:
                    if str(it.get('other')) == str(other_id) and str(it.get('master')) == str(mkey):
                        return True
                return False
            # Parse master spec to locate the intended combined_df row(s)
            def _parse_master(mkey: str):
                try:
                    pref, val = mkey.split(':', 1)
                except ValueError:
                    return None, None
                p = pref.strip().lower()
                cmap = {
                    'gaia': 'gaia_id',
                    '2mass': '2MASS', '2mass_j': '2MASS', '2massj': '2MASS', '2mass-j': '2MASS', '2massj-': '2MASS',
                    '2mass_id': '2MASS', '2massj_id': '2MASS',
                    'wise': 'wise_id',
                    'chandra': 'chandra_id',
                    'xmm': 'xmm_id',
                }
                col = cmap.get(p, None)
                return col, val

            # Determine real desired action on current dataset using the clicked row as primary hint:
            # if (other) is attached on the hinted row -> remove; otherwise -> force on that row,
            # regardless of attachments elsewhere.
            if action == 'toggle':
                try:
                    # Load current in-memory dataset for this catalog (without applying the new edit yet)
                    combined_df, catalogs = _load_data_for_key(cat)
                    # Validate that we have data to work with
                    if combined_df is None or combined_df.empty:
                        logger.warning(f"Empty dataframe for catalog {cat}, falling back to force mode")
                        action = 'force'
                        raise ValueError("Empty dataframe")
                    # Identify id column for the catalog
                    id_col_map = {'2MASS': '2MASS', 'wise': 'wise_id', 'chandra': 'chandra_id', 'xmm': 'xmm_id'}
                    id_col = id_col_map.get(cat, cat if cat in combined_df.columns else None)
                    if id_col is None or id_col not in combined_df.columns:
                        id_col = '2MASS' if cat == '2MASS' else (f"{cat}_id" if f"{cat}_id" in combined_df.columns else None)
                    mcol, mval = _parse_master(master)
                    # Resolve the hinted row (if valid and consistent with master key)
                    # Use iloc for positional access as row index may not be contiguous
                    hinted_row = None
                    if row and row.isdigit():
                        irow = int(row)
                        if (irow >= 0) and (irow < len(combined_df)):
                            actual_index = combined_df.index[irow]
                            hinted_row = actual_index
                            # If a master column is specified, ensure the hint matches the same master
                            if mcol and (mcol in combined_df.columns):
                                try:
                                    if str(combined_df.at[actual_index, mcol]) != str(mval):
                                        hinted_row = None
                                except Exception:
                                    hinted_row = None
                    # Evaluate attachment for this master/other pair
                    attached_row = None
                    if hinted_row is not None and id_col in combined_df.columns:
                        try:
                            same_master = True
                            if mcol and (mcol in combined_df.columns):
                                same_master = (str(combined_df.at[hinted_row, mcol]) == str(mval))
                            if same_master and str(combined_df.at[hinted_row, id_col]) == str(other):
                                attached_row = hinted_row
                        except Exception:
                            attached_row = None
                    if attached_row is None and id_col in combined_df.columns:
                        try:
                            mask = combined_df[id_col].astype(str) == str(other)
                            if mcol and (mcol in combined_df.columns):
                                mask = mask & (combined_df[mcol].astype(str) == str(mval))
                            idxs = [int(ix) for ix in combined_df.index[mask]]
                            if idxs:
                                attached_row = idxs[0]
                        except Exception:
                            attached_row = None
                    if attached_row is not None:
                        action = 'remove'
                        canonical_remove_row = str(attached_row)
                    else:
                        action = 'force'
                        force_row = None
                        if hinted_row is not None:
                            force_row = hinted_row
                            if mcol and (mcol in combined_df.columns):
                                try:
                                    if str(combined_df.at[hinted_row, mcol]) != str(mval):
                                        force_row = None
                                except Exception:
                                    force_row = None
                        if force_row is None and mcol and (mcol in combined_df.columns):
                            try:
                                m_mask = combined_df[mcol].astype(str) == str(mval)
                                idxs = [int(ix) for ix in combined_df.index[m_mask]]
                                if idxs:
                                    force_row = idxs[0]
                            except Exception:
                                force_row = None
                        if force_row is not None:
                            canonical_force_row = str(force_row)
                except Exception:
                    # Fallback to legacy behaviour
                    action = 'force'
            # Apply action
            if action == 'toggle':
                # Legacy toggle: fall back to force/remove list flipping
                if _contains(sect.get('force', []), other, master):
                    sect['force'] = [it for it in sect.get('force', []) if not (str(it.get('other')) == str(other) and str(it.get('master')) == str(master))]
                elif _contains(sect.get('remove', []), other, master):
                    sect['remove'] = [it for it in sect.get('remove', []) if not (str(it.get('other')) == str(other) and str(it.get('master')) == str(master))]
                else:
                    rec = {'other': other, 'master': master}
                    if row:
                        rec['row'] = row
                    sect.setdefault('force', []).append(rec)
            elif action == 'force':
                if not _contains(sect.get('force', []), other, master):
                    rec = {'other': other, 'master': master}
                    row_hint = row or canonical_force_row or None
                    if row_hint:
                        rec['row'] = str(row_hint)
                    sect.setdefault('force', []).append(rec)
                # Remove any remove entry
                sect['remove'] = [it for it in sect.get('remove', []) if not (str(it.get('other')) == str(other) and str(it.get('master')) == str(master))]
            elif action == 'remove':
                if not _contains(sect.get('remove', []), other, master):
                    rec = {'other': other, 'master': master}
                    row_hint = row or canonical_remove_row or None
                    if row_hint:
                        rec['row'] = str(row_hint)
                    sect.setdefault('remove', []).append(rec)
                # Remove any force entry
                sect['force'] = [it for it in sect.get('force', []) if not (str(it.get('other')) == str(other) and str(it.get('master')) == str(master))]
            # Persist
            with open(edits_path, 'w') as f:
                _json.dump(data, f, indent=2)
        except Exception as e:
            logger.exception("Error persisting edit links")
            return {'ok': False, 'error': 'Failed to save edit'}, 500

        # Reapply edits and respond
        try:
            for cache_key, entry in list(DATA.items()):
                base_df = entry.get('baseline')
                catalogs_entry = entry.get('catalogs')
                if base_df is None or catalogs_entry is None:
                    continue
                try:
                    updated = cm.apply_link_edits_to_combined(base_df, catalogs_entry, data)
                except Exception:
                    updated = base_df.copy(deep=True)
                    updated.attrs['_baseline_unforced'] = base_df
                entry['combined'] = updated
            # After edits, force a refresh of the page content
            _ensure_pages_for(cat, page, draw_images=True, prefetch=False)
            wants_json = False
            try:
                best = request.accept_mimetypes.best
            except Exception:
                best = None
            try:
                wants_json = (
                    (best == 'application/json') or
                    (request.headers.get('X-Requested-With') == 'XMLHttpRequest') or
                    (request.accept_mimetypes.get('application/json', 0) >= request.accept_mimetypes.get('text/html', 0))
                )
            except Exception:
                wants_json = False
            if wants_json:
                return {'ok': True}
            ref_url = request.referrer
            if ref_url:
                return redirect(ref_url, code=303)
            return {'ok': True}
        except Exception as e:
            logger.exception("Error applying edit links")
            return {'ok': False, 'error': 'Failed to apply edits'}, 500

    # Provide a compatibility route for the static index reference used by generated pages
    @app.get('/aladin_index.html')
    def serve_static_index():
        # If a static aladin_index.html exists at project root, serve it; otherwise redirect to '/'
        idx = (BASE_DIR / 'aladin_index.html')
        if idx.exists():
            return send_from_directory(str(BASE_DIR), 'aladin_index.html')
        return redirect(url_for('index'), code=302)

    @app.teardown_appcontext
    def shutdown_executor(exception=None):
        """Cleanup resources when app context ends."""
        try:
            PREFETCH_EXEC.shutdown(wait=False)
        except Exception as e:
            logger.error(f"Error shutting down executor: {e}")


    @app.route('/api/key_debug', methods=['GET', 'POST'])
    def api_key_debug():
        key = (request.values.get('key') or '').strip()
        code = (request.values.get('code') or '').strip()
        ctrl = request.values.get('ctrl', '0')
        alt = request.values.get('alt', '0')
        shift = request.values.get('shift', '0')
        meta = request.values.get('meta', '0')
        ts = request.values.get('ts', '')
        try:
            print(f"[key_debug] key={key!r} code={code!r} ctrl={ctrl} alt={alt} shift={shift} meta={meta} ts={ts}")
        except Exception:
            pass
        return {'ok': True}

    # Optional convenience: forward dynamic SAMP sends through this server
    @app.get('/run_samp')
    def run_samp():
        base = (request.args.get('file') or '').strip()
        if not base:
            abort(400)
        if os.path.sep in base or (os.path.altsep and os.path.altsep in base):
            abort(400)
        ajs_path = (SCRIPTS_DIR / base)
        if not _is_safe_path(SCRIPTS_DIR, ajs_path) or not ajs_path.exists():
            abort(404)
        ok, msg = run_samp_script(str(ajs_path), addr=os.environ.get('SAMP_ADDR', '127.0.0.1'))
        if ok:
            return render_template_string('<html><body><h3>Sent to Aladin</h3><pre>{{name}}</pre></body></html>', name=base)
        else:
            return render_template_string('<html><body><h3>Failed</h3><pre>{{msg}}</pre></body></html>', msg=msg), 500

    return app


if __name__ == '__main__':
    # Simple CLI for host/port and radius (falls back to environment variables)
    import argparse
    parser = argparse.ArgumentParser(description='Dynamic NGC 2264 pages server')
    parser.add_argument('--host', default=os.environ.get('HOST', '127.0.0.1'),
                        help='Host/IP to bind (default from HOST env or 127.0.0.1)')
    parser.add_argument('--port', type=int, default=int(os.environ.get('PORT', '5050')),
                        help='Port to listen on (default from PORT env or 5050)')
    parser.add_argument('--radius-deg', type=float, default=None,
                        help='Cone search radius in degrees for all catalogs (overrides DYN_RADIUS_DEG env)')
    parser.add_argument('--refresh', action='store_true', default=False,
                        help='Ignore cached CSVs on first build (sets DYN_REFRESH=1)')
    # Image/HTML rendering quality controls
    parser.add_argument('--html-dpi', type=int, default=None,
                        help='Figure DPI used for rasterization (env HTML_DPI)')
    parser.add_argument('--html-scale', type=float, default=None,
                        help='Display scale factor in HTML (env HTML_SCALE)')
    parser.add_argument('--html-raster-scale', type=float, default=None,
                        help='Additional raster downscale factor (env HTML_RASTER_SCALE)')
    parser.add_argument('--html-image-format', choices=['png','jpg','jpeg','webp'], default=None,
                        help='Raster format for page images (env HTML_IMAGE_FORMAT)')
    parser.add_argument('--html-jpeg-quality', type=int, default=None,
                        help='JPEG quality 20–95 (env HTML_JPEG_QUALITY)')
    parser.add_argument('--html-webp-quality', type=int, default=None,
                        help='WEBP quality 30–95 (env HTML_WEBP_QUALITY)')
    # Background 2MASS J image resolution controls
    parser.add_argument('--tmass-fetch', choices=['auto','hips','skyview','off'], default='hips',
                        help='Select image provider for 2MASS J backgrounds (env TMASS_FETCH, default hips)')
    parser.add_argument('--tmass-arcsec-per-pix', type=float, default=None,
                        help='2MASS J HiPS pixel scale in arcsec/pixel (env TMASS_ARCSEC_PER_PIX)')
    parser.add_argument('--tmass-max-size', type=int, default=None,
                        help='2MASS J HiPS maximum output size (px) (env TMASS_MAX_SIZE)')
    parser.add_argument('--no-images', action='store_true', default=False,
                        help='Disable background images (sets DYN_NO_IMAGES=1)')
    args = parser.parse_args()

    # Export CLI options to environment so downstream builder sees them
    if args.radius_deg is not None:
        os.environ['DYN_RADIUS_DEG'] = str(args.radius_deg)
    if args.refresh:
        os.environ['DYN_REFRESH'] = '1'
    if args.no_images:
        os.environ['DYN_NO_IMAGES'] = '1'
    if args.html_dpi is not None:
        os.environ['HTML_DPI'] = str(int(args.html_dpi))
    if args.html_scale is not None:
        os.environ['HTML_SCALE'] = str(float(args.html_scale))
    if args.html_raster_scale is not None:
        os.environ['HTML_RASTER_SCALE'] = str(float(args.html_raster_scale))
    if args.html_image_format is not None:
        os.environ['HTML_IMAGE_FORMAT'] = args.html_image_format
    if args.html_jpeg_quality is not None:
        os.environ['HTML_JPEG_QUALITY'] = str(int(args.html_jpeg_quality))
    if args.html_webp_quality is not None:
        os.environ['HTML_WEBP_QUALITY'] = str(int(args.html_webp_quality))
    if args.tmass_fetch is not None:
        os.environ['TMASS_FETCH'] = args.tmass_fetch
    if args.tmass_arcsec_per_pix is not None:
        os.environ['TMASS_ARCSEC_PER_PIX'] = str(float(args.tmass_arcsec_per_pix))
    if args.tmass_max_size is not None:
        os.environ['TMASS_MAX_SIZE'] = str(int(args.tmass_max_size))

    # Set higher-quality defaults unless overridden by CLI or existing env
    # These defaults match your preferred settings: 300 DPI, PNG, no extra downscale.
    os.environ.setdefault('HTML_DPI', '300')
    os.environ.setdefault('HTML_IMAGE_FORMAT', 'png')
    os.environ.setdefault('HTML_RASTER_SCALE', '1.0')

    port = int(args.port)
    host = str(args.host)
    # Make sure the expected env for the generator points to our scripts dir
    os.environ.setdefault('ALADIN_SCRIPTS_DIR', str(SCRIPTS_DIR))
    app = create_app()

    # Only enable debug mode if explicitly requested
    debug_mode = os.environ.get('FLASK_DEBUG', '0') in ('1', 'true', 'yes')

    print(f"[dyn] Serving dynamic pages at http://{host}:{port}/ (scripts={SCRIPTS_DIR}, debug={debug_mode})")
    app.run(host=host, port=port, debug=debug_mode)
