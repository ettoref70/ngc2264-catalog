import argparse
import os
from typing import Iterable

from webapp.dyn_server import create_app


def iter_pages(start: int, end: int | None, count: int | None) -> Iterable[int]:
    if end is not None:
        if end < start:
            raise ValueError("end must be >= start")
        return range(start, end + 1)
    if count is None or count < 1:
        raise ValueError("count must be a positive integer when end is not provided")
    return range(start, start + count)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Pre-generate Aladin HTML pages using the dynamic renderer."
    )
    parser.add_argument(
        "--catalogs",
        nargs="+",
        default=["2MASS"],
        help="Catalog keys to pre-render (default: 2MASS).",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=1,
        help="First page number to render (default: 1).",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--end",
        type=int,
        help="Last page number to render (inclusive).",
    )
    group.add_argument(
        "--count",
        type=int,
        default=400,
        help="Number of pages to render starting from --start when --end is not provided (default: 400).",
    )
    parser.add_argument(
        "--no-refresh-first",
        action="store_true",
        help="Do not force a refresh on the first page (assumes caches already populated).",
    )
    parser.add_argument(
        "--radius-deg",
        type=float,
        help="Override catalog radius in degrees (sets DYN_RADIUS_DEG).",
    )
    args = parser.parse_args()

    # Match dynamic server defaults so raster quality is consistent
    os.environ.setdefault("HTML_IMAGE_FORMAT", "png")
    os.environ.setdefault("HTML_RASTER_SCALE", "1.0")
    os.environ.setdefault("HTML_JPEG_QUALITY", "90")
    if args.radius_deg is not None:
        os.environ["DYN_RADIUS_DEG"] = str(args.radius_deg)

    app = create_app()
    ensure_pages_for = getattr(app, "ensure_pages_for", None)
    load_data_for_key = getattr(app, "load_data_for_key", None)
    if ensure_pages_for is None or load_data_for_key is None:
        raise RuntimeError("create_app() did not expose page generation helpers")

    with app.app_context():
        pages = list(iter_pages(args.start, args.end, args.count))
        if not pages:
            print("No pages to build.")
            return

        for catalog in args.catalogs:
            print(f"== {catalog}: preparing {len(pages)} page(s) ==")
            try:
                if args.no_refresh_first:
                    os.environ["DYN_REFRESH"] = "0"
                else:
                    os.environ["DYN_REFRESH"] = "1"
                    load_data_for_key(catalog)
                    os.environ["DYN_REFRESH"] = "0"
            except Exception as exc:
                print(f"{catalog}: failed to load data ({exc})")
                continue

            for idx, page in enumerate(pages):
                try:
                    ensure_pages_for(catalog, page, draw_images=True, prefetch=False)
                    print(f"{catalog} page {page} built")
                except Exception as exc:
                    print(f"{catalog} page {page} failed: {exc}")
                    break


if __name__ == "__main__":
    main()
