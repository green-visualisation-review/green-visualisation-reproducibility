from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple, Dict, Any

import fitz  # PyMuPDF
import cv2
import numpy as np
import pandas as pd
from tqdm import tqdm


@dataclass(frozen=True)
class GreenHSVConfig:
    # Hue range in DEGREES (0..360). OpenCV uses hue 0..179 internally.
    hue_low_deg: float = 70.0      # ~yellow-green boundary
    hue_high_deg: float = 170.0    # includes green->teal/cyan-ish

    # Thresholds to avoid counting white/very dark pixels as green
    green_s_min: int = 40          # 0..255
    green_v_min: int = 40          # 0..255

    # "White background" definition (to exclude backgrounds)
    white_s_max: int = 25          # low saturation
    white_v_min: int = 230         # high value

    # "Colored pixels" definition (to compute green share among colors)
    colored_s_min: int = 40
    colored_v_min: int = 40


def _deg_to_cv_hue(h_deg: float) -> int:
    """Convert degrees (0..360) -> OpenCV hue (0..179)."""
    h = int(round(h_deg / 2.0))
    return max(0, min(179, h))


def _hue_mask(h: np.ndarray, low: int, high: int) -> np.ndarray:
    """
    Hue mask for OpenCV hue channel in [0..179].
    Handles wrap-around ranges (e.g., 170..10).
    """
    if low <= high:
        return (h >= low) & (h <= high)
    # wrap-around case
    return (h >= low) | (h <= high)


def render_page_rgb(page: fitz.Page, dpi: int = 150, max_dim: int | None = 1600) -> np.ndarray:
    """
    Render a PDF page to an RGB numpy array.
    - dpi controls rendering detail
    - max_dim optionally downsizes to keep processing fast
    """
    zoom = dpi / 72.0
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, alpha=False)  # RGB
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)

    if pix.n == 4:
        img = img[:, :, :3]

    # PyMuPDF gives RGB already
    rgb = img

    if max_dim is not None:
        h, w = rgb.shape[:2]
        scale = max(h / max_dim, w / max_dim)
        if scale > 1.0:
            new_w = int(round(w / scale))
            new_h = int(round(h / scale))
            rgb = cv2.resize(rgb, (new_w, new_h), interpolation=cv2.INTER_AREA)

    return rgb


def compute_green_metrics(rgb: np.ndarray, cfg: GreenHSVConfig) -> Dict[str, float]:
    """
    Compute green metrics from an RGB image using HSV thresholds.

    Returns:
      - green_coverage_all: green pixels / all pixels
      - green_share_foreground: green pixels / non-white pixels
      - green_share_colored: green pixels / "colored" pixels
    """
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    h, s, v = cv2.split(hsv)

    low_h = _deg_to_cv_hue(cfg.hue_low_deg)
    high_h = _deg_to_cv_hue(cfg.hue_high_deg)

    is_green_hue = _hue_mask(h, low_h, high_h)
    is_green = is_green_hue & (s >= cfg.green_s_min) & (v >= cfg.green_v_min)

    is_white_bg = (s <= cfg.white_s_max) & (v >= cfg.white_v_min)
    is_foreground = ~is_white_bg

    is_colored = (s >= cfg.colored_s_min) & (v >= cfg.colored_v_min)

    total = rgb.shape[0] * rgb.shape[1]
    green_count = int(is_green.sum())
    fg_count = int(is_foreground.sum())
    colored_count = int(is_colored.sum())

    green_coverage_all = green_count / total if total else 0.0
    green_share_foreground = green_count / fg_count if fg_count else 0.0
    green_share_colored = green_count / colored_count if colored_count else 0.0

    return {
        "green_coverage_all": float(green_coverage_all),
        "green_share_foreground": float(green_share_foreground),
        "green_share_colored": float(green_share_colored),
    }


def iter_pdfs(folder: Path) -> Iterable[Path]:
    yield from folder.rglob("*.pdf")


def process_folder(
    folder: Path,
    cfg: GreenHSVConfig,
    dpi: int = 150,
    max_dim: int | None = 1600,
    max_pages_per_pdf: int | None = None,
    page_stride: int = 1,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      df_pages: one row per processed page
      df_reports: one row per PDF (summary + baseline)
    """
    rows: list[Dict[str, Any]] = []

    pdfs = sorted(iter_pdfs(folder))
    for pdf_path in tqdm(pdfs, desc="PDFs"):
        try:
            doc = fitz.open(pdf_path)
        except Exception as e:
            rows.append({
                "pdf": str(pdf_path),
                "page_index": None,
                "error": f"open_failed: {e}",
            })
            continue

        n_pages = doc.page_count
        page_indices = list(range(0, n_pages, page_stride))
        if max_pages_per_pdf is not None:
            page_indices = page_indices[:max_pages_per_pdf]

        for i in page_indices:
            try:
                page = doc.load_page(i)
                rgb = render_page_rgb(page, dpi=dpi, max_dim=max_dim)
                metrics = compute_green_metrics(rgb, cfg)
                rows.append({
                    "pdf": str(pdf_path),
                    "page_index": int(i),
                    "pages_in_pdf": int(n_pages),
                    "error": "",
                    **metrics,
                })
            except Exception as e:
                rows.append({
                    "pdf": str(pdf_path),
                    "page_index": int(i),
                    "pages_in_pdf": int(n_pages),
                    "error": f"page_failed: {e}",
                })

        doc.close()

    df_pages = pd.DataFrame(rows)

    # Keep only successful pages for summary stats
    ok = df_pages["error"].fillna("") == ""
    df_ok = df_pages[ok].copy()

    if df_ok.empty:
        df_reports = pd.DataFrame(columns=[
            "pdf", "pages_processed",
            "baseline_green_share_colored_median",
            "baseline_green_share_foreground_median",
            "green_share_colored_mean",
            "green_share_foreground_mean",
            "green_coverage_all_mean",
        ])
        return df_pages, df_reports

    # Per-PDF baseline: median across pages (robust to a few very-green pages)
    g = df_ok.groupby("pdf", as_index=False)
    df_reports = g.agg(
        pages_processed=("page_index", "count"),
        baseline_green_share_colored_median=("green_share_colored", "median"),
        baseline_green_share_foreground_median=("green_share_foreground", "median"),
        green_share_colored_mean=("green_share_colored", "mean"),
        green_share_foreground_mean=("green_share_foreground", "mean"),
        green_coverage_all_mean=("green_coverage_all", "mean"),
    ).sort_values("baseline_green_share_colored_median", ascending=False)

    return df_pages, df_reports


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate a baseline level of green (HSV) from PDFs in a folder."
    )

    # ✅ folder 改成可选参数（nargs="?"），并给默认路径
    parser.add_argument(
        "folder",
        nargs="?",
        default="/Users/claire/Library/CloudStorage/OneDrive-Personal/Colour/chatpgt/bar/RR1/CSR report",
        type=str,
        help="Folder containing PDFs (will search recursively). Quote it if it contains spaces.",
    )

    parser.add_argument("--dpi", type=int, default=150, help="Render DPI for PDF pages.")
    parser.add_argument(
        "--max-dim",
        type=int,
        default=1600,
        help="Downscale longest side to this many pixels (0 disables).",
    )
    parser.add_argument(
        "--max-pages-per-pdf",
        type=int,
        default=0,
        help="0 means process all pages; otherwise limit.",
    )
    parser.add_argument("--page-stride", type=int, default=1, help="Process every Nth page (1 = all).")

    # Green definition
    parser.add_argument("--hue-low", type=float, default=70.0, help="Green hue low bound in degrees (0..360).")
    parser.add_argument("--hue-high", type=float, default=170.0, help="Green hue high bound in degrees (0..360).")
    parser.add_argument("--green-s-min", type=int, default=40, help="Min saturation for green pixels (0..255).")
    parser.add_argument("--green-v-min", type=int, default=40, help="Min value for green pixels (0..255).")

    # White background definition
    parser.add_argument("--white-s-max", type=int, default=25, help="Max saturation to treat as white background.")
    parser.add_argument("--white-v-min", type=int, default=230, help="Min value to treat as white background.")

    # Colored pixel definition
    parser.add_argument("--colored-s-min", type=int, default=40, help="Min saturation to treat as colored pixel.")
    parser.add_argument("--colored-v-min", type=int, default=40, help="Min value to treat as colored pixel.")

    parser.add_argument("--out-prefix", type=str, default="green_baseline", help="Prefix for output CSV files.")

    args = parser.parse_args()

    folder = Path(args.folder).expanduser()
    if not folder.exists():
        raise SystemExit(f"Folder not found: {folder}")

    cfg = GreenHSVConfig(
        hue_low_deg=args.hue_low,
        hue_high_deg=args.hue_high,
        green_s_min=args.green_s_min,
        green_v_min=args.green_v_min,
        white_s_max=args.white_s_max,
        white_v_min=args.white_v_min,
        colored_s_min=args.colored_s_min,
        colored_v_min=args.colored_v_min,
    )

    max_dim = None if args.max_dim <= 0 else args.max_dim
    max_pages = None if args.max_pages_per_pdf <= 0 else args.max_pages_per_pdf

    df_pages, df_reports = process_folder(
        folder=folder,
        cfg=cfg,
        dpi=args.dpi,
        max_dim=max_dim,
        max_pages_per_pdf=max_pages,
        page_stride=args.page_stride,
    )

    pages_csv = f"{args.out_prefix}_pages.csv"
    reports_csv = f"{args.out_prefix}_reports.csv"
    df_pages.to_csv(pages_csv, index=False)
    df_reports.to_csv(reports_csv, index=False)

    print(f"\nWrote:\n  {pages_csv}\n  {reports_csv}\n")

    if not df_reports.empty:
        ok = df_pages["error"].fillna("") == ""
        overall_baseline = float(df_pages.loc[ok, "green_share_colored"].median())
        print(f"Overall baseline (median green_share_colored across all pages): {overall_baseline:.6f}")


if __name__ == "__main__":
    main()

