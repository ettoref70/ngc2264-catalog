FITS Web Viewer (on-the-fly)
===========================

This lightweight Flask app serves the already-downloaded FITS files under
`catalogs/` and provides a browser UI with the JS9 FITS viewer. Nothing is
pre-generated: the pages and previews are built on the fly.

Quick start
-----------

1) Install dependencies (in your environment):

   pip install flask

   (JS9 is loaded from a CDN in the HTML; no Python dependency is needed.)

2) Run the server from the repo root:

   python webapp/app.py

   By default it serves at http://127.0.0.1:5000 and exposes `catalogs/` at
   `http://127.0.0.1:5000/fits/...`.

3) Open the UI:

   http://127.0.0.1:5000/

Notes
-----

- JS9 fetches FITS via HTTP range requests, so files load efficiently without
  creating PNG/JPEGs server-side.
- The Browse page paginates large directories (default 200 files per page).
- Paths are restricted to the `catalogs/` subtree for safety.
- To change host/port: `HOST=0.0.0.0 PORT=8000 python webapp/app.py`.

