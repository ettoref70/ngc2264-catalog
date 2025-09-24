from __future__ import annotations

import os
from pathlib import Path
from typing import List

from flask import (
    Flask,
    abort,
    redirect,
    render_template,
    request,
    send_from_directory,
    url_for,
)
import mimetypes


BASE_DIR = Path(__file__).resolve().parent.parent
FITS_ROOT = (BASE_DIR / "catalogs").resolve()


def _is_safe_path(root: Path, target: Path) -> bool:
    try:
        root = root.resolve()
        target = target.resolve()
    except FileNotFoundError:
        # If target doesn't exist yet, check parent
        target = target.parent.resolve()
    return str(target).startswith(str(root))


def _list_dirs_and_fits(rel_dir: str) -> tuple[List[str], List[str]]:
    """Return (subdirectories, fits files) for a directory under FITS_ROOT.

    Paths are returned relative to FITS_ROOT using POSIX separators.
    """
    base = FITS_ROOT / rel_dir
    if not _is_safe_path(FITS_ROOT, base):
        abort(403)
    if not base.exists() or not base.is_dir():
        abort(404)
    dirs: List[str] = []
    fits: List[str] = []
    for entry in sorted(base.iterdir()):
        if entry.name.startswith('.'):
            continue
        if entry.is_dir():
            # ensure dir has some fits within (optional optimization)
            dirs.append((Path(rel_dir) / entry.name).as_posix())
        else:
            n = entry.name.lower()
            if n.endswith('.fits') or n.endswith('.fits.gz'):
                fits.append((Path(rel_dir) / entry.name).as_posix())
    return dirs, fits


def create_app() -> Flask:
    app = Flask(__name__)
    # Ensure correct MIME for WebAssembly so JS9 can load astroemw.wasm
    try:
        mimetypes.add_type('application/wasm', '.wasm')
    except Exception:
        pass

    @app.get("/")
    def index():
        # List top-level subdirectories under FITS_ROOT to start browsing
        if not FITS_ROOT.exists():
            return render_template("error.html", title="Missing catalogs/", message=f"Not found: {FITS_ROOT}")
        top_dirs = [p.name for p in sorted(FITS_ROOT.iterdir()) if p.is_dir() and not p.name.startswith('.')]  # noqa: E501
        return render_template("index.html", top_dirs=top_dirs)

    @app.get("/browse")
    def browse():
        rel_dir = request.args.get("dir", "").strip()
        page = max(int(request.args.get("page", 1) or 1), 1)
        per_page = min(max(int(request.args.get("per", 200) or 200), 10), 1000)
        dirs, fits = _list_dirs_and_fits(rel_dir)
        total = len(fits)
        start = (page - 1) * per_page
        end = start + per_page
        page_fits = fits[start:end]
        has_next = end < total
        has_prev = start > 0
        return render_template(
            "browse.html",
            rel_dir=rel_dir,
            subdirs=dirs,
            files=page_fits,
            page=page,
            per=per_page,
            total=total,
            has_next=has_next,
            has_prev=has_prev,
        )

    @app.get("/view")
    def view():
        rel_path = request.args.get("file", "").strip()
        if not rel_path:
            return redirect(url_for("index"))
        fpath = (FITS_ROOT / rel_path)
        if not _is_safe_path(FITS_ROOT, fpath):
            abort(403)
        if not fpath.exists():
            abort(404)
        # Build an absolute URL so clients don't misresolve relative paths
        fits_url = url_for("serve_fits", subpath=rel_path, _external=True)
        return render_template("view.html", rel_path=rel_path, fits_url=fits_url)

    @app.get("/fits/<path:subpath>")
    def serve_fits(subpath: str):
        # Serve FITS from catalogs/ with safety checks
        fpath = (FITS_ROOT / subpath)
        if not _is_safe_path(FITS_ROOT, fpath):
            abort(403)
        if not fpath.exists() or not fpath.is_file():
            abort(404)
        # Send using directory+filename to avoid exposing full path
        directory = fpath.parent
        filename = fpath.name
        # Let JS9 use range requests; avoid caching issues during development
        resp = send_from_directory(directory, filename, as_attachment=False, conditional=True)
        # Help JS9/worker with range requests and content type
        resp.headers.setdefault("Accept-Ranges", "bytes")
        # FITS is typically 'application/fits' though 'application/octet-stream' also works
        if filename.lower().endswith((".fits", ".fits.gz")):
            resp.headers.setdefault("Content-Type", "application/fits")
        resp.headers.setdefault("Cache-Control", "public, max-age=3600")
        return resp

    # Convenience: if someone requests a FITS path without the /fits prefix,
    # and it exists under the catalogs/ tree, redirect to the canonical URL.
    @app.get("/<path:subpath>")
    def redirect_plain_fits(subpath: str):
        # Only handle .fits or .fits.gz
        low = subpath.lower()
        if not (low.endswith('.fits') or low.endswith('.fits.gz')):
            abort(404)
        # Map to catalogs/ and validate
        fpath = (FITS_ROOT / subpath)
        if not _is_safe_path(FITS_ROOT, fpath) or not fpath.exists():
            abort(404)
        return redirect(url_for('serve_fits', subpath=subpath, _external=True), code=302)

    return app


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5000"))
    host = os.environ.get("HOST", "127.0.0.1")
    app = create_app()
    print(f"[webapp] Serving FITS from: {FITS_ROOT}")
    app.run(host=host, port=port, debug=True)
