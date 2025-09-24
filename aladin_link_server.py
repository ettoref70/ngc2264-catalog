#!/usr/bin/env python3
"""
aladin_link_server.py
---------------------

Tiny local HTTP server to bridge HTML clicks to Aladin via SAMP.

- Serves /run_samp?file=<ajs-basename>, loads the .ajs from --scripts-dir
  and sends it to Aladin Desktop through SAMP (Aladin must have SAMP connected).

Usage:
  python aladin_link_server.py --port 8765 --scripts-dir aladin_scripts [--addr 127.0.0.1]

Set ALADIN_LINK_URL in browser localStorage to the server URL if you want
to change it from the default (http://127.0.0.1:8765) used by the HTML pages:

  localStorage.setItem('ALADIN_LINK_URL','http://127.0.0.1:8765')
"""
from __future__ import annotations

import argparse
import html
import json
import os
import sys
import urllib.parse
import mimetypes
import io
import hashlib
from http.server import BaseHTTPRequestHandler, HTTPServer
from socketserver import ThreadingMixIn
import tempfile


def run_samp_script(ajs_path: str, addr: str = "127.0.0.1") -> tuple[bool, str]:
    try:
        import subprocess
        py = os.environ.get("PYTHON", sys.executable or "python3")
        cmd = [py, os.path.abspath("aladin_samp_test.py"), "--mode", "script", "--send-ajs", ajs_path, "--addr", addr]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=15)
        if proc.returncode == 0:
            return True, "OK"
        else:
            return False, f"samp_helper_failed: rc={proc.returncode}, out={proc.stdout}\nerr={proc.stderr}"
    except Exception as e:
        return False, f"exception: {e}"


class ThreadingHTTPServer(ThreadingMixIn, HTTPServer):
    daemon_threads = True


class Handler(BaseHTTPRequestHandler):
    server_version = "AladinLink/1.1"

    def do_GET(self):  # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path in ("/health", "/ping"):
            return self._send_json({"ok": True, "msg": "alive"})
        if parsed.path == "/run_samp":
            return self._handle_run_samp(parsed)
        if parsed.path == "/aladin":
            return self._handle_aladin_dynamic(parsed)
        if parsed.path == "/img":
            return self._handle_render_image(parsed)
        if parsed.path in ("/", "/index", "/index.html"):
            return self._handle_root()
        if parsed.path == "/legacy-index":
            return self._handle_legacy_index()
        if parsed.path.startswith("/aladin_scripts/"):
            return self._serve_aladin_static(parsed.path)
        self.send_response(404)
        self.end_headers()

    def log_message(self, fmt, *args):  # quiet
        sys.stderr.write("%s - - %s\n" % (self.address_string(), fmt % args))

    def _send_html(self, body: str, code: int = 200):
        data = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _send_json(self, obj: dict, code: int = 200):
        data = json.dumps(obj).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _send_bytes(self, data: bytes, content_type: str, code: int = 200):
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _handle_root(self):
        scripts_dir = os.path.abspath(getattr(self.server, "scripts_dir", "aladin_scripts"))  # type: ignore[attr-defined]
        legacy = os.path.abspath(os.path.join(os.getcwd(), "aladin_index.html"))
        has_legacy = os.path.exists(legacy)
        body = f"""
        <html><head><meta charset=\"utf-8\"><title>Aladin Link</title>
        <style>body{{font:14px -apple-system,BlinkMacSystemFont,Segoe UI,Arial}} code{{background:#f4f4f4;padding:2px 4px;border-radius:3px}}</style>
        </head><body>
        <h2>Aladin Link Server</h2>
        <p>Endpoints:</p>
        <ul>
          <li><code>/health</code> — server health</li>
          <li><code>/aladin</code> — dynamic AJS sender (params: ra, dec, fov, fov_unit, survey, fits, dryrun)</li>
          <li><code>/run_samp?file=NAME.ajs</code> — send an existing script from <code>{html.escape(scripts_dir)}</code></li>
          <li><code>/aladin_scripts/html/...</code> — serve existing HTML pages under <code>{html.escape(scripts_dir)}/html</code></li>
        </ul>
        <p>Quick test (dryrun preview):<br>
           <a href=\"/aladin?ra=100.25&dec=9.8833&fov=6&survey=CDS/P/2MASS/J&dryrun=1\">Center on RA=100.25 Dec=9.8833 (2MASS/J)</a>
        </p>
        {('<p>Legacy pages: <a href=\"/legacy-index\">aladin_index.html</a></p>' if has_legacy else '')}
        </body></html>
        """
        return self._send_html(body)

    def _handle_legacy_index(self):
        p = os.path.abspath(os.path.join(os.getcwd(), "aladin_index.html"))
        if not os.path.exists(p):
            return self._send_html("No aladin_index.html present", 404)
        try:
            with open(p, "rb") as f:
                data = f.read()
            return self._send_bytes(data, "text/html; charset=utf-8")
        except Exception as e:
            return self._send_html(f"Failed to read aladin_index.html: {html.escape(str(e))}", 500)

    def _serve_aladin_static(self, req_path: str):
        # Map /aladin_scripts/... to the configured scripts_dir safely
        base = os.path.abspath(getattr(self.server, "scripts_dir", "aladin_scripts"))  # type: ignore[attr-defined]
        # strip prefix
        rel = req_path[len("/aladin_scripts/"):]
        local = os.path.abspath(os.path.join(base, rel))
        if not local.startswith(base + os.path.sep):
            return self._send_html("Forbidden", 403)
        if not os.path.exists(local) or not os.path.isfile(local):
            return self._send_html("Not found", 404)
        ctype = mimetypes.guess_type(local)[0] or "application/octet-stream"
        if ctype.startswith("text/"):
            ctype += "; charset=utf-8"
        try:
            with open(local, "rb") as f:
                data = f.read()
            return self._send_bytes(data, ctype)
        except Exception as e:
            return self._send_html(f"Failed to read file: {html.escape(str(e))}", 500)

    def _handle_run_samp(self, parsed):
        qs = urllib.parse.parse_qs(parsed.query)
        base = (qs.get("file", [""])[0] or "").strip()
        scripts_dir = self.server.scripts_dir  # type: ignore[attr-defined]
        addr = getattr(self.server, "samp_addr", "127.0.0.1")  # type: ignore[attr-defined]
        if not base:
            return self._send_html("<b>Missing</b> ?file= parameter", 400)
        if os.path.sep in base or (os.path.altsep and os.path.altsep in base):
            return self._send_html("Invalid file parameter", 400)
        ajs_path = os.path.abspath(os.path.join(scripts_dir, base))
        if not ajs_path.startswith(os.path.abspath(scripts_dir) + os.path.sep):
            return self._send_html("Forbidden", 403)
        if not os.path.exists(ajs_path):
            return self._send_html("Not found: " + html.escape(base), 404)
        ok, msg = run_samp_script(ajs_path, addr=addr)
        if ok:
            return self._send_html(f"<html><body><h3>Sent to Aladin</h3><pre>{html.escape(base)}</pre></body></html>")
        else:
            return self._send_html(f"<html><body><h3>Failed</h3><pre>{html.escape(msg)}</pre></body></html>", 500)

    # ----------------------- Dynamic AJS generator -----------------------
    def _handle_aladin_dynamic(self, parsed):
        qs = urllib.parse.parse_qs(parsed.query)
        def getf(name: str, default: float | None = None) -> float | None:
            v = qs.get(name, [None])[0]
            if v is None or v == "":
                return default
            try:
                return float(v)
            except Exception:
                return default

        ra = getf("ra")
        dec = getf("dec")
        # fov unit: arcmin by default; allow fov_unit=deg|arcmin|arcsec
        fov = getf("fov", None)
        fov_unit = (qs.get("fov_unit", ["arcmin"])[0] or "arcmin").lower()
        survey = (qs.get("survey", ["CDS/P/DSS2/color"])[0] or "CDS/P/DSS2/color").strip()
        fits_param = (qs.get("fits", [""])[0] or "").strip()
        dryrun = (qs.get("dryrun", ["0"])[0] or "0") in ("1", "true", "yes")

        # Build AJS script lines
        ajs: list[str] = []
        ajs.append("reset")
        if survey:
            ajs.append(f"get hips(\"{survey}\")")
        # Center, if given
        if (ra is not None) and (dec is not None):
            ajs.append(f"{ra:.7f} {dec:.7f}")
        # Zoom, if given
        if fov is not None:
            unit = "arcmin" if fov_unit not in ("deg", "arcsec", "arcmin") else fov_unit
            ajs.append(f"zoom {fov:.3f} {unit}")

        # Optional FITS load (local file under fits_root or remote URL)
        if fits_param:
            load_arg = None
            # Allow http/https URLs
            if fits_param.lower().startswith(("http://", "https://")):
                load_arg = fits_param
            else:
                # Treat as path under configured fits_root
                fits_root = getattr(self.server, "fits_root", None)  # type: ignore[attr-defined]
                if fits_root:
                    requested = os.path.abspath(os.path.join(fits_root, fits_param))
                    fits_root_abs = os.path.abspath(fits_root)
                    if requested.startswith(fits_root_abs + os.path.sep) and os.path.exists(requested):
                        load_arg = requested
                    else:
                        return self._send_html("Invalid or not found FITS path", 400)
                else:
                    return self._send_html("FITS loading disabled (no --fits-root configured)", 400)
            ajs.append(f"load \"{load_arg}\"")

        ajs_text = "\n".join(ajs) + "\n"

        if dryrun:
            body = "<html><body><h3>Dynamic AJS (dryrun)</h3><pre>" + html.escape(ajs_text) + "</pre>" \
                   + "<p>Append &dryrun=0 to send to Aladin.</p></body></html>"
            return self._send_html(body)

        # Write to temp file and send via SAMP helper
        try:
            fd, tmpname = tempfile.mkstemp(prefix="aladin_dynamic_", suffix=".ajs")
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                f.write(ajs_text)
        except Exception as e:
            return self._send_html(f"Failed to create temp AJS: {html.escape(str(e))}", 500)

        addr = getattr(self.server, "samp_addr", "127.0.0.1")  # type: ignore[attr-defined]
        ok, msg = run_samp_script(tmpname, addr=addr)
        try:
            os.unlink(tmpname)
        except Exception:
            pass
        if ok:
            info = {
                "ok": True,
                "sent": True,
                "survey": survey,
                "center": {"ra": ra, "dec": dec} if (ra is not None and dec is not None) else None,
                "fov": {"value": fov, "unit": fov_unit} if fov is not None else None,
                "fits": fits_param or None,
            }
            return self._send_html("<html><body><h3>Sent to Aladin</h3><pre>" + html.escape(ajs_text) + "</pre></body></html>")
        else:
            return self._send_html(f"<html><body><h3>Failed</h3><pre>{html.escape(msg)}</pre><pre>{html.escape(ajs_text)}</pre></body></html>", 500)

    # ----------------------- FITS -> JPEG renderer -----------------------
    def _handle_render_image(self, parsed):
        """Render a FITS file to a JPEG/PNG with optional DS9 region overlays.

        Query parameters:
          - fits: relative path under fits_root (required)
          - reg:  region file name (repeatable), resolved under scripts_dir by default
          - size: output image width in pixels (default 1024)
          - fmt:  'jpg' (default) or 'png'
          - stretch: 'asinh' (default) or 'linear'
          - percent: percentile for scaling (default 99.5)
          - cmap: matplotlib colormap name (default 'gray')
        """
        qs = urllib.parse.parse_qs(parsed.query)
        fits_rel = (qs.get("fits", [""])[0] or "").strip()
        if not fits_rel:
            return self._send_html("Missing ?fits= parameter", 400)
        fits_root = getattr(self.server, "fits_root", None)  # type: ignore[attr-defined]
        if not fits_root:
            return self._send_html("FITS root not configured", 500)
        fits_path = os.path.abspath(os.path.join(fits_root, fits_rel))
        if not fits_path.startswith(os.path.abspath(fits_root) + os.path.sep):
            return self._send_html("Forbidden", 403)
        if not os.path.exists(fits_path):
            return self._send_html("Not found: " + html.escape(fits_rel), 404)

        # Optional regions (resolved under scripts_dir unless absolute)
        reg_list = qs.get("reg", [])
        reg_paths: list[str] = []
        scripts_dir = getattr(self.server, "scripts_dir", os.getcwd())  # type: ignore[attr-defined]
        for r in reg_list:
            r = (r or "").strip()
            if not r:
                continue
            if os.path.isabs(r) or r.lower().startswith(("http://", "https://")):
                p = r
            else:
                p = os.path.abspath(os.path.join(scripts_dir, r))
            if os.path.exists(p):
                reg_paths.append(p)

        # Glob pattern support: regpat=aladin_2MASS_7_*.reg (repeatable)
        import glob as _glob
        for pat in qs.get("regpat", []):
            pat = (pat or "").strip()
            if not pat:
                continue
            g = os.path.abspath(os.path.join(scripts_dir, pat))
            for p in sorted(_glob.glob(g)):
                if os.path.isfile(p):
                    reg_paths.append(p)

        try:
            size = int((qs.get("size", ["1024"])[0] or "1024"))
            size = max(128, min(size, 4096))
        except Exception:
            size = 1024
        fmt = (qs.get("fmt", ["jpg"])[0] or "jpg").lower()
        if fmt not in ("jpg", "jpeg", "png"):
            fmt = "jpg"
        stretch = (qs.get("stretch", ["asinh"])[0] or "asinh").lower()
        try:
            percent = float((qs.get("percent", ["99.5"])[0] or "99.5"))
        except Exception:
            percent = 99.5
        cmap = (qs.get("cmap", ["gray"])[0] or "gray")
        reg_color = (qs.get("regcolor", ["cyan"])[0] or "cyan")
        try:
            reg_width = float((qs.get("regwidth", ["1.5"])[0] or "1.5"))
        except Exception:
            reg_width = 1.5

        # Disk cache key
        cache_root = os.path.abspath(getattr(self.server, "render_cache", os.path.join(os.getcwd(), "docs_cache", "render")))  # type: ignore[attr-defined]
        os.makedirs(cache_root, exist_ok=True)
        key_src = json.dumps({
            "fits": fits_path,
            "regs": reg_paths,
            "size": size,
            "fmt": fmt,
            "stretch": stretch,
            "percent": percent,
            "cmap": cmap,
            "reg_color": reg_color,
            "reg_width": reg_width,
        }, sort_keys=True).encode("utf-8")
        key = hashlib.sha1(key_src).hexdigest()[:20]
        cache_file = os.path.join(cache_root, f"{key}.{fmt}")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "rb") as f:
                    data = f.read()
                ctype = "image/jpeg" if fmt in ("jpg", "jpeg") else "image/png"
                self.send_response(200)
                self.send_header("Content-Type", ctype)
                self.send_header("Cache-Control", "public, max-age=86400")
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                self.wfile.write(data)
                return
            except Exception:
                pass

        # Render on the fly
        try:
            img_bytes = self._render_fits_core(
                fits_path=fits_path,
                reg_paths=reg_paths,
                size=size,
                fmt=fmt,
                stretch=stretch,
                percent=percent,
                cmap=cmap,
                reg_color=reg_color,
                reg_width=reg_width,
            )
        except Exception as e:
            return self._send_html("Render failed: " + html.escape(str(e)), 500)

        try:
            with open(cache_file, "wb") as f:
                f.write(img_bytes)
        except Exception:
            pass
        ctype = "image/jpeg" if fmt in ("jpg", "jpeg") else "image/png"
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Cache-Control", "public, max-age=86400")
        self.send_header("Content-Length", str(len(img_bytes)))
        self.end_headers()
        self.wfile.write(img_bytes)

    def _render_fits_core(self, fits_path: str, reg_paths: list[str], size: int, fmt: str, stretch: str, percent: float, cmap: str, reg_color: str = "cyan", reg_width: float = 1.5) -> bytes:
        import numpy as np
        from astropy.io import fits as _fits
        from astropy.wcs import WCS as _WCS
        from astropy.visualization import ImageNormalize, AsinhStretch, LinearStretch, PercentileInterval, ZScaleInterval
        import matplotlib
        matplotlib.use("Agg")  # headless
        import matplotlib.pyplot as plt
        # Load data
        with _fits.open(fits_path, memmap=True) as hdul:
            hdu = None
            for h in hdul:
                if getattr(h, 'data', None) is not None and isinstance(h.data, np.ndarray) and h.data.ndim == 2:
                    hdu = h
                    break
            if hdu is None:
                raise RuntimeError("no 2D image HDU found")
            data = np.asarray(hdu.data, dtype=float)
            hdr = hdu.header
        # Clean data
        data = np.nan_to_num(data, copy=False)
        # Normalize
        if stretch == 'linear':
            inter = PercentileInterval(percent)
            st = LinearStretch()
        else:
            # default asinh + percentile, fallback to zscale
            try:
                inter = PercentileInterval(percent)
                st = AsinhStretch()
            except Exception:
                inter = ZScaleInterval()
                st = AsinhStretch()
        vmin, vmax = inter.get_limits(data)
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=st)
        # Figure sizing
        dpi = 100.0
        w = h = float(size) / dpi
        wcs = None
        try:
            wcs = _WCS(hdr)
        except Exception:
            wcs = None
        fig = plt.figure(figsize=(w, h), dpi=dpi)
        if wcs is not None and wcs.is_celestial:
            ax = fig.add_subplot(111, projection=wcs)
        else:
            ax = fig.add_subplot(111)
        im = ax.imshow(data, origin='lower', cmap=cmap, norm=norm, interpolation='nearest')
        ax.set_axis_off()
        # Overlay regions when available
        if reg_paths:
            try:
                from regions import Regions
            except Exception as e:
                raise RuntimeError("regions package not installed; install with 'pip install regions'") from e
            regs = []
            for rp in reg_paths:
                try:
                    rset = Regions.read(rp, format='ds9')
                    regs.extend(list(rset))
                except Exception:
                    continue
            # Build a list of pixel regions to draw
            pixel_regions = []
            for r in regs:
                # Try to convert sky regions to pixel using WCS; keep pixel regions as-is
                converted = None
                if wcs is not None:
                    try:
                        converted = r.to_pixel(wcs)  # returns PixelRegion for SkyRegion; PixelRegion may raise
                    except Exception:
                        converted = None
                if converted is not None:
                    pixel_regions.append(converted)
                else:
                    # If it's already a pixel region, use it directly
                    try:
                        # crude duck-typing: PixelRegion usually has to_sky and as_artist(origin=...)
                        if hasattr(r, 'to_sky') and hasattr(r, 'as_artist'):
                            pixel_regions.append(r)
                    except Exception:
                        pass
            # Draw pixel regions
            for pr in pixel_regions:
                try:
                    artist = pr.as_artist(origin='lower')
                    if hasattr(artist, 'set_zorder'):
                        artist.set_zorder(1000)
                    # style
                    try:
                        if hasattr(artist, 'set_edgecolor'):
                            artist.set_edgecolor(reg_color)
                        if hasattr(artist, 'set_linewidth'):
                            artist.set_linewidth(reg_width)
                        if hasattr(artist, 'set_facecolor'):
                            artist.set_facecolor('none')
                        if hasattr(artist, 'set_color'):
                            artist.set_color(reg_color)
                        if hasattr(artist, 'set_alpha'):
                            artist.set_alpha(0.9)
                    except Exception:
                        pass
                    ax.add_artist(artist)
                except Exception:
                    pass
        buf = io.BytesIO()
        if fmt in ("jpg", "jpeg"):
            # Some Matplotlib versions do not accept a 'quality' kw; use defaults
            fig.savefig(buf, format='jpeg', dpi=dpi, bbox_inches='tight', pad_inches=0)
        else:
            fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
        return buf.getvalue()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Local link server for Aladin HTML pages")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--scripts-dir", default="aladin_scripts")
    ap.add_argument("--addr", default="127.0.0.1", help="SAMP host address for the client")
    ap.add_argument("--fits-root", default="catalogs", help="Root directory for local FITS files (used by /aladin?fits=...) ")
    args = ap.parse_args(argv)

    scripts_dir = os.path.abspath(args.scripts_dir)
    if not os.path.isdir(scripts_dir):
        print(f"[error] scripts-dir not found: {scripts_dir}")
        return 2

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    setattr(server, "scripts_dir", scripts_dir)
    setattr(server, "samp_addr", args.addr)
    setattr(server, "fits_root", os.path.abspath(args.fits_root))
    print(f"[aladin-link] http://{args.host}:{args.port} (scripts={scripts_dir}, fits_root={getattr(server, 'fits_root', '')})")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
