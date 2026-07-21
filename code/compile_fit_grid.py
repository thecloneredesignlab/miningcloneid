"""
compile_fit_grid.py
────────────────────
Combines every mouse's final-fit plot into a single tiled 8x2 (8 rows,
2 columns = 16 mice) image, so all mice's fits can be compared side by side
at a glance.

By default it tiles observed_fit_tumor_burden.png — the fitted (MAP) model
simulated on the *same drugs the mouse actually received*, overlaid on the
observed tumor-burden data points. That is the plot to look at to judge
whether the model fits the data, because the predicted curve and the data
points are directly comparable.

(This is intentionally NOT baseline_tumor_burden.png. That plot runs the
fitted model forward on the beam-search *optimal* drug course — a
counterfactual "what treatment should we have given" — so its curve is not
expected to pass through the observed points. To build the grid from that
plot instead, pass --image-name baseline_tumor_burden.png.)

Unlike compile_gif_grid.py, this only cares about the *final* static fit,
not the gibbs_fitting_animation.gif that shows the MCMC fitting process
converging frame by frame. Output is therefore a single static image, not a
GIF.

If observed_fit_tumor_burden.png isn't found for a mouse, this falls back to
baseline_tumor_burden.png so the grid still fills for older result dirs that
predate that plot. Regenerate the observed-drug plots for already-completed
runs with regenerate_observed_fits.py.

Additional shared-y-bounds grid
───────────────────────────────
By default this ALSO writes a second grid (final_fit_grid_uniform.png) built
from observed_fit_tumor_burden_uniform.png — the per-mouse copies that
regenerate_observed_fits.py renders with identical y-axis bounds across every
mouse. In that grid the panels are directly comparable in height because they
share one y-axis range, whereas the main grid's panels are each autoscaled.
The uniform grid does NOT fall back to the baseline plot (falling back would
break the shared bounds); mice missing the uniform image are left blank, so
run regenerate_observed_fits.py first. Pass --no-uniform to skip this second
grid.

Where it looks for each mouse's plot (first match wins), for each candidate
image name in turn:
    1. results_joint/<harvest_name>/<image_name>
       (joint-fit mode — run_joint.py / submit_joint.sbatch)
    2. harvest_zips/<harvest_name>.zip containing <image_name>
       (parallel per-harvest mode — run_all_harvests.py, after zipping)
    3. results_<harvest_name>/<image_name>
       (per-harvest mode, before zipping / cleanup)

Usage:
    python compile_fit_grid.py
    python compile_fit_grid.py --rows 8 --cols 2 --tile-width 480
    python compile_fit_grid.py --output harvest_zips/final_fit_grid.png
    python compile_fit_grid.py --mice mouseA mouseB ...   # explicit order
    python compile_fit_grid.py --image-name baseline_tumor_burden.png  # old plot
    python compile_fit_grid.py --no-uniform               # skip shared-bounds grid
    python compile_fit_grid.py --uniform-output out.png   # rename shared-bounds grid
"""

import argparse
import csv
import io
import os
import zipfile
from typing import Optional

from PIL import Image, ImageDraw, ImageFont

# ── Paths (match run_all_harvests.py / beam_search_flip_rate_wgd.py) ────────
CSV_PATH        = "../data/InVivoData_Gemcitabine/harvest_ploidy_mapping.csv"
RESULTS_JOINT   = "results_joint"
HARVEST_ZIPS    = "harvest_zips"
RESULTS_PREFIX  = "results"   # results_<harvest>
IMAGE_NAME      = "observed_fit_tumor_burden.png"   # fitted model on the drugs actually given
FALLBACK_IMAGE  = "baseline_tumor_burden.png"       # used only if IMAGE_NAME is missing
UNIFORM_IMAGE   = "observed_fit_tumor_burden_uniform.png"  # shared-y-bounds copies (regenerate_observed_fits.py)


def get_matching_harvests(csv_path: str) -> list[str]:
    harvests = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row["has_match"].strip() == "True":
                harvests.append(row["harvest"].strip())
    return harvests


def locate_image_bytes(harvest: str, image_names: list[str]) -> tuple[Optional[bytes], Optional[str]]:
    """
    Find a mouse's fit plot and return (raw bytes, which image name matched).

    image_names is tried in order (e.g. the observed-drug fit first, then the
    baseline plot as a fallback); for each name every location is checked.
    """
    for image_name in image_names:
        # 1. Joint-fit mode
        p = os.path.join(RESULTS_JOINT, harvest, image_name)
        if os.path.isfile(p):
            with open(p, "rb") as fh:
                return fh.read(), image_name

        # 2. Zipped per-harvest results
        zpath = os.path.join(HARVEST_ZIPS, f"{harvest}.zip")
        if os.path.isfile(zpath):
            with zipfile.ZipFile(zpath) as zf:
                for name in zf.namelist():
                    if os.path.basename(name) == image_name:
                        return zf.read(name), image_name

        # 3. Unzipped per-harvest results dir
        p = os.path.join(f"{RESULTS_PREFIX}_{harvest}", image_name)
        if os.path.isfile(p):
            with open(p, "rb") as fh:
                return fh.read(), image_name

    return None, None


def make_label_font(size: int):
    try:
        return ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", size
        )
    except OSError:
        return ImageFont.load_default()


def build_grid(mice: list[str], image_names: list[str], args,
               output_path: str, header_text: str) -> bool:
    """
    Load each mouse's plot (first name in *image_names* that exists wins),
    tile them into an args.rows x args.cols grid, and save to *output_path*.

    Returns True if the grid was written, False if no source images were found
    (in which case nothing is saved and the caller can decide how to report it).
    """
    n_slots = args.rows * args.cols

    tiles: list[Optional[Image.Image]] = []
    tile_size: Optional[tuple[int, int]] = None

    font = make_label_font(max(12, args.tile_width // 24)) if args.label else None
    for harvest in mice:
        img_bytes, matched_name = locate_image_bytes(harvest, image_names)
        if img_bytes is None:
            print(f"  MISSING: {harvest} — none of {image_names} found, leaving slot blank")
            tiles.append(None)
            continue

        img = Image.open(io.BytesIO(img_bytes)).convert("RGB")
        note = "" if matched_name == image_names[0] else f"  (fallback: {matched_name})"
        print(f"  Loaded {harvest}: {img.size[0]}x{img.size[1]}{note}")

        if tile_size is None:
            w0, h0 = img.size
            tile_h = int(round(args.tile_width * h0 / w0))
            tile_size = (args.tile_width, tile_h)

        img = img.resize(tile_size, Image.LANCZOS)
        if args.label:
            img = img.copy()
            draw = ImageDraw.Draw(img)
            text = harvest
            pad = 4
            bbox = draw.textbbox((0, 0), text, font=font)
            tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
            draw.rectangle([0, 0, tw + 2 * pad, th + 2 * pad], fill=(0, 0, 0))
            draw.text((pad, pad), text, fill=(255, 255, 255), font=font)
        tiles.append(img)

    if tile_size is None:
        return False

    blank_tile = Image.new("RGB", tile_size, (30, 30, 30))
    for i, t in enumerate(tiles):
        if t is None:
            tiles[i] = blank_tile

    header_font = None
    if args.header:
        header_font = make_label_font(max(16, int(args.header_height * 0.42)))

    # ── Compose the grid ─────────────────────────────────────────────────
    canvas_w = args.cols * tile_size[0]
    grid_h   = args.rows * tile_size[1]
    header_h = args.header_height if args.header else 0
    canvas_h = grid_h + header_h

    canvas = Image.new("RGB", (canvas_w, canvas_h), (0, 0, 0))
    for idx in range(n_slots):
        row, col = divmod(idx, args.cols)
        tile = tiles[idx] if idx < len(tiles) else blank_tile
        canvas.paste(tile, (col * tile_size[0], header_h + row * tile_size[1]))

    if args.header:
        draw = ImageDraw.Draw(canvas)
        bbox = draw.textbbox((0, 0), header_text, font=header_font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        draw.text(((canvas_w - tw) / 2, (header_h - th) / 2 - bbox[1]),
                   header_text, fill=(255, 255, 255), font=header_font)

    # ── Save ──────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    canvas.save(output_path)
    print(f"Saved compilation image: {output_path}  "
          f"({canvas_w}x{canvas_h}px, {args.rows}x{args.cols} grid"
          f"{', with header' if args.header else ''})")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=CSV_PATH,
                     help="harvest_ploidy_mapping.csv path (used for default mouse order)")
    ap.add_argument("--mice", nargs="+", default=None,
                     help="Explicit ordered list of harvest names to include "
                          "(overrides --csv autodiscovery)")
    ap.add_argument("--rows", type=int, default=8)
    ap.add_argument("--cols", type=int, default=2)
    ap.add_argument("--tile-width", type=int, default=480,
                     help="Width in px each mouse's panel is resized to before tiling")
    ap.add_argument("--image-name", default=IMAGE_NAME,
                     help=f"Which per-mouse plot to tile (default: {IMAGE_NAME}, "
                          f"the fitted model on the drugs actually given)")
    ap.add_argument("--no-fallback", dest="fallback", action="store_false", default=True,
                     help=f"Don't fall back to {FALLBACK_IMAGE} when --image-name "
                          f"is missing for a mouse")
    ap.add_argument("--no-label", dest="label", action="store_false", default=True,
                     help="Disable the harvest-name label overlay on each panel "
                          "(the plot title already names the harvest, so this is "
                          "an extra top-left tag for quick scanning)")
    ap.add_argument("--no-header", dest="header", action="store_false", default=True,
                     help="Disable the top title banner")
    ap.add_argument("--header-height", type=int, default=70,
                     help="Height in px of the top title banner")
    ap.add_argument("--header-text",
                     default="Fitted Model on Observed Drugs vs Data — All Harvests",
                     help="Text shown in the top title banner")
    ap.add_argument("--output", default=os.path.join(HARVEST_ZIPS, "final_fit_grid.png"))
    ap.add_argument("--no-uniform", dest="uniform", action="store_false", default=True,
                     help="Skip the additional shared-y-bounds grid built from "
                          f"{UNIFORM_IMAGE}")
    ap.add_argument("--uniform-image-name", default=UNIFORM_IMAGE,
                     help=f"Per-mouse image to tile for the shared-bounds grid "
                          f"(default: {UNIFORM_IMAGE})")
    ap.add_argument("--uniform-output",
                     default=os.path.join(HARVEST_ZIPS, "final_fit_grid_uniform.png"),
                     help="Output path for the shared-y-bounds grid")
    ap.add_argument("--uniform-header-text",
                     default="Fitted Model on Observed Drugs vs Data — Shared Y-Axis — All Harvests",
                     help="Top-banner text for the shared-y-bounds grid")
    args = ap.parse_args()

    n_slots = args.rows * args.cols

    # ── Determine mouse order ────────────────────────────────────────────
    if args.mice:
        mice = args.mice
    else:
        mice = get_matching_harvests(args.csv)

    if len(mice) != n_slots:
        print(f"WARNING: found {len(mice)} mice but grid needs {n_slots} "
              f"({args.rows}x{args.cols}). ", end="")
        if len(mice) > n_slots:
            print(f"Truncating to first {n_slots}.")
            mice = mice[:n_slots]
        else:
            print("Remaining grid slots will be left blank.")

    # ── Main grid: candidate image names (primary first, optional fallback) ──
    image_names = [args.image_name]
    if args.fallback and FALLBACK_IMAGE != args.image_name:
        image_names.append(FALLBACK_IMAGE)

    print("Building main grid...")
    ok = build_grid(mice, image_names, args, args.output, args.header_text)
    if not ok:
        raise RuntimeError(
            f"No mouse fit images could be found (looked for {image_names}) — "
            f"nothing to compile. If your results predate observed_fit_tumor_burden.png, "
            f"run regenerate_observed_fits.py first, or pass "
            f"--image-name baseline_tumor_burden.png."
        )

    # ── Additional grid: same panels but every mouse on shared y-bounds ──────
    # No fallback here — falling back to a differently-scaled plot would defeat
    # the point of a shared y-axis, so missing panels are simply left blank.
    if args.uniform:
        print("\nBuilding shared-y-bounds grid...")
        uni_ok = build_grid(mice, [args.uniform_image_name], args,
                            args.uniform_output, args.uniform_header_text)
        if not uni_ok:
            print(f"  No {args.uniform_image_name} images found — shared-bounds grid "
                  f"not written. Run regenerate_observed_fits.py first to create them.")


if __name__ == "__main__":
    main()