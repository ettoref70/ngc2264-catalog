from __future__ import annotations

import base64
import os
from pathlib import Path
from typing import Dict

from flask import Flask, abort, redirect, render_template_string, request, send_from_directory, url_for
import sys


BASE_DIR = Path(__file__).resolve().parent.parent
# Ensure project root is on sys.path so we can import local modules when launched as a submodule
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

# Local modules (after adjusting sys.path)
import cross_match_ngc2264 as cm  # type: ignore
from aladin_link_server import run_samp_script  # type: ignore
SCRIPTS_DIR = (BASE_DIR / "aladin_scripts").resolve()
HTML_DIR = (SCRIPTS_DIR / "html").resolve()


def _is_safe_path(root: Path, target: Path) -> bool:
    try:
        root = root.resolve()
        target = target.resolve()
    except FileNotFoundError:
        target = target.parent.resolve()
    return str(target).startswith(str(root))


def create_app() -> Flask:
    app = Flask(__name__)

    # Cache datasets per requested catalog key (e.g., '2MASS', 'wise')
    DATA: Dict[str, tuple] = {}

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
        # Return cached subset if present
        if key in DATA:
            return DATA[key]
        # Try fast path from prebuilt combined CSV
        fast = _try_build_from_combined(key)
        if fast is not None:
            DATA[key] = fast
            return fast
        # Fallback to building from queries, limited to gaia + requested key
        refresh = (os.environ.get('DYN_REFRESH', '0') in ('1', 'true', 'yes'))
        include = ['gaia', '2MASS'] if key.lower() in ('2mass', '2mass_j', '2massj') else ['gaia', key]
        combined_df, catalogs = cm.build_data_for_web(refresh=refresh,
                                                      chandra_csv_path=os.environ.get('CHANDRA_CSV_PATH', '/Users/ettoref/ASTRONOMY/DATA/N2264_XMM_alt/N2264_acis12.csv'),
                                                      include_catalogs=include)
        DATA[key] = (combined_df, catalogs)
        return DATA[key]

    def _ensure_pages_for(key: str, page_num: int) -> None:
        combined_df, catalogs = _load_data_for_key(key)
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
        # Default to no images for responsiveness unless explicitly enabled.
        draw_images = not (os.environ.get('DYN_NO_IMAGES', '1') in ('1', 'true', 'yes'))
        # Generate only the requested page for this catalog (writes to aladin_scripts/html)
        # Use env var understood by the plotter to limit generation to a single page
        os.environ['ALADIN_PAGE_ONLY'] = str(int(page_num))
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
        _ensure_pages_for(key, page)
        fname = f"{key}_page{page}.html"
        fpath = HTML_DIR / fname
        if not fpath.exists():
            abort(404)
        # Redirect to the canonical static path so relative links resolve
        return redirect(url_for('serve_generated', subpath=fname), code=302)

    # Serve generated assets under /aladin_scripts/html/...
    @app.get('/aladin_scripts/html/<path:subpath>')
    def serve_generated(subpath: str):
        fpath = (HTML_DIR / subpath)
        if not _is_safe_path(HTML_DIR, fpath) or not fpath.exists():
            abort(404)
        return send_from_directory(HTML_DIR, subpath)

    # Also expose the entire aladin_scripts tree (e.g., .ajs and regions)
    @app.get('/aladin_scripts/<path:subpath>')
    def serve_scripts(subpath: str):
        fpath = (SCRIPTS_DIR / subpath)
        if not _is_safe_path(SCRIPTS_DIR, fpath) or not fpath.exists():
            abort(404)
        return send_from_directory(SCRIPTS_DIR, subpath)

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
    port = int(os.environ.get('PORT', '5050'))
    host = os.environ.get('HOST', '127.0.0.1')
    # Make sure the expected env for the generator points to our scripts dir
    os.environ.setdefault('ALADIN_SCRIPTS_DIR', str(SCRIPTS_DIR))
    app = create_app()
    print(f"[dyn] Serving dynamic pages at http://{host}:{port}/ (scripts={SCRIPTS_DIR})")
    app.run(host=host, port=port, debug=True)
