#!/usr/bin/env python3
"""
aladin_samp_test.py

Minimal SAMP test to drive Aladin Desktop from Python.

Usage examples:

  python aladin_samp_test.py --ra 100.25 --dec 9.883333 --fov 6 \
      --survey "CDS/P/DSS2/color"

Notes:
- Start Aladin Desktop first and ensure SAMP is connected (Interop/SAMP panel).
- Requires astropy (for astropy.samp): pip install astropy
"""

from __future__ import annotations

import sys
import argparse

try:
    from astropy.samp import SAMPIntegratedClient  # type: ignore
except Exception:
    print("[error] astropy.samp not available. Install with: pip install astropy")
    sys.exit(1)


def find_aladin_clients(client: SAMPIntegratedClient) -> list[str]:
    cids: list[str] = []
    for cid in client.get_registered_clients():
        try:
            meta = client.get_metadata(cid)
            name = " ".join(str(meta.get(k, "")) for k in ("samp.name", "application.name", "author")).lower()
            if ("aladin" in name) or ("cds" in name):
                cids.append(cid)
        except Exception:
            continue
    return cids


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Send a simple script to Aladin via SAMP")
    ap.add_argument("--ra", type=float, default=100.25, help="ICRS RA in degrees")
    ap.add_argument("--dec", type=float, default=9.883333, help="ICRS Dec in degrees")
    ap.add_argument("--fov", type=float, default=6.0, help="FoV in arcminutes")
    ap.add_argument(
        "--survey",
        default="CDS/P/DSS2/color",
        help="HiPS survey id (e.g. 'CDS/P/2MASS/J' or 'CDS/P/AllWISE/Color')",
    )
    ap.add_argument("--verbose", action="store_true", help="Print extra diagnostics")
    ap.add_argument(
        "--mode",
        choices=["center", "script"],
        default="center",
        help="Use safe center-only SAMP message (default) or send a script (may be blocked)",
    )
    ap.add_argument("--overlay", action="store_true",
                    help="Send a VOTable overlay (points) via SAMP after centering")
    ap.add_argument("--overlay-file", default="", help="Path to an existing VOTable to send via SAMP")
    ap.add_argument("--gaia-ra", type=float, default=None, help="Optional Gaia RA for overlay (deg)")
    ap.add_argument("--gaia-dec", type=float, default=None, help="Optional Gaia Dec for overlay (deg)")
    ap.add_argument("--pmra", type=float, default=None, help="pmRA* mas/yr for PM endpoint (optional)")
    ap.add_argument("--pmdec", type=float, default=None, help="pmDec mas/yr for PM endpoint (optional)")
    ap.add_argument("--dt", type=float, default=None, help="Time baseline in years for PM endpoint (optional)")
    ap.add_argument(
        "--addr",
        default="127.0.0.1",
        help="Local address to bind the SAMP client to (try '127.0.0.1' or '::1')",
    )
    ap.add_argument(
        "--send-ajs",
        default="",
        help="Path to a local .ajs file to send to Aladin via SAMP as a script",
    )
    args = ap.parse_args(argv)

    # Try to bind explicitly to a local interface to avoid 'Can't assign requested address'
    cli = SAMPIntegratedClient(name="Aladin-SAMP-Test", addr=args.addr)
    try:
        cli.connect()
    except Exception as e:
        if args.addr != "::1":
            print(f"[warn] Connect failed on addr='{args.addr}': {e}. Retrying with IPv6 '::1' …")
            try:
                cli = SAMPIntegratedClient(name="Aladin-SAMP-Test", addr="::1")
                cli.connect()
            except Exception as e2:
                print(f"[error] Could not connect to SAMP hub: {e2}")
                print("- Ensure Aladin Desktop is open and SAMP is Connected (Interop panel)")
                print("- Some systems require binding to 127.0.0.1 or ::1; try --addr 127.0.0.1 or --addr ::1")
                return 2
        else:
            print(f"[error] Could not connect to SAMP hub: {e}")
            print("- Ensure Aladin Desktop is open and SAMP is Connected (Interop panel)")
            print("- Some systems require binding to 127.0.0.1 or ::1; try --addr 127.0.0.1 or --addr ::1")
            return 2

    if args.verbose:
        print("[SAMP] Registered clients:")
        for cid in cli.get_registered_clients():
            try:
                print(" -", cid, cli.get_metadata(cid).get("samp.name"))
            except Exception:
                print(" -", cid)

    aladin_cids = find_aladin_clients(cli)
    if not aladin_cids:
        print("[warn] No Aladin client detected via SAMP. Is Aladin open and SAMP connected?")
        cli.disconnect()
        return 1

    sent = False
    if args.mode == "script":
        # Prepare script content: either from --send-ajs or a small generated script
        if args.send_ajs:
            try:
                with open(args.send_ajs, "r") as f:
                    script_text = f.read()
            except Exception as e:
                print(f"[error] Could not read --send-ajs '{args.send_ajs}': {e}")
                script_text = ""
        else:
            ajs = [
                "reset",
                f"info \"SAMP test script\"",
                f"get hips(\"{args.survey}\")",
                f"{args.ra:.7f} {args.dec:.7f}",
                f"zoom {args.fov:.3f} arcmin",
            ]
            script_text = "\n".join(ajs)

        if script_text:
            for mtype in ("script.aladin.send", "aladin.script", "clientcmd"):
                try:
                    msg = {"samp.mtype": mtype, "samp.params": {"script": script_text}}
                    for cid in aladin_cids:
                        cli.notify(cid, msg)
                    if args.verbose:
                        print(f"[SAMP] Sent script via mtype='{mtype}' to {aladin_cids}")
                    sent = True
                    break
                except Exception as e:
                    if args.verbose:
                        print(f"[SAMP] mtype '{mtype}' failed: {e}")

    if not sent:
        # Safe default: just center on the target (widely supported)
        try:
            msg = {"samp.mtype": "coord.pointAt.sky", "samp.params": {"ra": f"{args.ra:.7f}", "dec": f"{args.dec:.7f}"}}
            for cid in aladin_cids:
                cli.notify(cid, msg)
            print("[SAMP] Sent coord.pointAt.sky to", aladin_cids)
            sent = True
        except Exception as e:
            print("[error] Could not send any SAMP message to Aladin:", e)
            cli.disconnect()
            return 3

    # Optionally send an overlay as a VOTable (points only)
    if args.overlay:
        import os, tempfile
        vot_path = args.overlay_file
        if not vot_path:
            # Build a minimal VOTable with 1–2 points: center and optional PM endpoint or Gaia
            ra_c = float(args.gaia_ra) if (args.gaia_ra is not None) else float(args.ra)
            de_c = float(args.gaia_dec) if (args.gaia_dec is not None) else float(args.dec)
            rows = [("center", f"{args.ra:.7f}", f"{args.dec:.7f}")]
            # Add GAIA point if different from center
            if (abs(ra_c - args.ra) > 1e-8) or (abs(de_c - args.dec) > 1e-8):
                rows.append(("gaia", f"{ra_c:.7f}", f"{de_c:.7f}"))
            # PM endpoint
            if (args.pmra is not None) and (args.pmdec is not None) and (args.dt is not None):
                # Convert pm to deg over dt; pmra is mu_alpha* (RA*cosDec)
                dra_deg = (float(args.pmra) / 3600e3) * (1.0 / max(1e-8, __import__('math').cos(__import__('math').radians(de_c)))) * float(args.dt)
                ddec_deg = (float(args.pmdec) / 3600e3) * float(args.dt)
                rows.append(("pm_end", f"{ra_c + dra_deg:.7f}", f"{de_c + ddec_deg:.7f}"))

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
            for name, ra_s, de_s in rows:
                vot.append(f"     <TR><TD>{name}</TD><TD>{ra_s}</TD><TD>{de_s}</TD></TR>")
            vot += [
                "    </TABLEDATA>",
                "   </DATA>",
                "  </TABLE>",
                " </RESOURCE>",
                "</VOTABLE>",
            ]
            fd, tmpname = tempfile.mkstemp(prefix="aladin_overlay_", suffix=".vot")
            with os.fdopen(fd, "w") as f:
                f.write("\n".join(vot))
            vot_path = tmpname

        url = "file://" + os.path.abspath(vot_path)
        try:
            msg = {"samp.mtype": "table.load.votable", "samp.params": {"url": url}}
            for cid in aladin_cids:
                cli.notify(cid, msg)
            print("[SAMP] Sent table.load.votable:", url)
        except Exception as e:
            print("[warn] Could not send VOTable overlay:", e)

    cli.disconnect()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
