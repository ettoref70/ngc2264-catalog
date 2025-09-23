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
from http.server import BaseHTTPRequestHandler, HTTPServer
from socketserver import ThreadingMixIn


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


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Local link server for Aladin HTML pages")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--scripts-dir", default="aladin_scripts")
    ap.add_argument("--addr", default="127.0.0.1", help="SAMP host address for the client")
    args = ap.parse_args(argv)

    scripts_dir = os.path.abspath(args.scripts_dir)
    if not os.path.isdir(scripts_dir):
        print(f"[error] scripts-dir not found: {scripts_dir}")
        return 2

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    setattr(server, "scripts_dir", scripts_dir)
    setattr(server, "samp_addr", args.addr)
    print(f"[aladin-link] http://{args.host}:{args.port} (scripts={scripts_dir})")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

