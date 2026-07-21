"""
compile_gif_grid.py
────────────────────
Combines every mouse's gibbs_fitting_animation.gif into a single tiled
8x2 (8 rows, 2 columns = 16 mice) compilation GIF, so all mice's fits can be
watched converging side by side.

Because every mouse in this project is fit for the same number of MCMC
iterations, each individual gibbs_fitting_animation.gif has the same number
of animation frames, and frame i in every mouse's GIF corresponds to the
same iteration index. This script simply reads frame i from every mouse's
GIF, tiles them into a grid, and writes that as frame i of the output GIF —
so all 16 sub-animations advance through their iterations in lockstep.

Where it looks for the per-mouse GIFs (first match wins):
    1. results_joint/<harvest_name>/gibbs_fitting_animation.gif
       (joint-fit mode — run_joint.py / submit_joint.sbatch)
    2. harvest_zips/<harvest_name>.zip containing gibbs_fitting_animation.gif
       (parallel per-harvest mode — run_all_harvests.py, after zipping)
    3. results_<harvest_name>/gibbs_fitting_animation.gif
       (per-harvest mode, before zipping / cleanup)

Usage:
    python compile_gif_grid.py
    python compile_gif_grid.py --rows 8 --cols 2 --tile-width 480
    python compile_gif_grid.py --output harvest_zips/gibbs_fitting_grid.gif
    python compile_gif_grid.py --mice mouseA mouseB ...   # explicit order
"""

import argparse
import csv
import io
import os
import zipfile
from typing import Optional

import numpy as np
from PIL import Image, ImageDraw, ImageFont

# ── Paths (match run_all_harvests.py / beam_search_flip_rate_wgd.py) ────────
CSV_PATH        = "../data/InVivoData_Gemcitabine/harvest_ploidy_mapping.csv"
RESULTS_JOINT   = "results_joint"
HARVEST_ZIPS    = "harvest_zips"
RESULTS_PREFIX  = "results"   # results_<harvest>
GIF_NAME        = "gibbs_fitting_animation.gif"

# ── MCMC iteration bookkeeping (must match mcmc_fit.py) ─────────────────────
# Each per-mouse gibbs_fitting_animation.gif is built by sampling
# min(200, N_BURNIN + N_SAMPLES) evenly-spaced points out of the full MCMC
# trace via np.linspace — see beam_search_flip_rate_wgd.py's
# "Fitting animation" section and mcmc_fit.py's all_trace.append() loop.
# To show the true iteration number on each compiled frame (rather than just
# the animation frame number) we reproduce that exact same linspace here.
# If you've changed N_BURNIN/N_SAMPLES in mcmc_fit.py, override them below
# via --n-burnin / --n-samples (or pass --total-iterations directly).
N_BURNIN_DEFAULT  = 6000
N_SAMPLES_DEFAULT = 14000


def get_matching_harvests(csv_path: str) -> list[str]:
    harvests = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row["has_match"].strip() == "True":
                harvests.append(row["harvest"].strip())
    return harvests


def locate_gif_bytes(harvest: str) -> Optional[bytes]:
    """Find a mouse's gibbs_fitting_animation.gif and return its raw bytes."""
    # 1. Joint-fit mode
    p = os.path.join(RESULTS_JOINT, harvest, GIF_NAME)
    if os.path.isfile(p):
        with open(p, "rb") as fh:
            return fh.read()

    # 2. Zipped per-harvest results
    zpath = os.path.join(HARVEST_ZIPS, f"{harvest}.zip")
    if os.path.isfile(zpath):
        with zipfile.ZipFile(zpath) as zf:
            for name in zf.namelist():
                if os.path.basename(name) == GIF_NAME:
                    return zf.read(name)

    # 3. Unzipped per-harvest results dir
    p = os.path.join(f"{RESULTS_PREFIX}_{harvest}", GIF_NAME)
    if os.path.isfile(p):
        with open(p, "rb") as fh:
            return fh.read()

    return None


def load_frames(gif_bytes: bytes) -> tuple[list[Image.Image], list[int]]:
    """Return (frames as RGB images, per-frame durations in ms)."""
    im = Image.open(io.BytesIO(gif_bytes))
    frames, durations = [], []
    try:
        i = 0
        while True:
            im.seek(i)
            frames.append(im.convert("RGB"))
            durations.append(im.info.get("duration", 80))
            i += 1
    except EOFError:
        pass
    return frames, durations


def make_label_font(size: int):
    try:
        return ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", size
        )
    except OSError:
        return ImageFont.load_default()


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
    ap.add_argument("--no-label", dest="label", action="store_false", default=True,
                     help="Disable the harvest-name label overlay on each panel")
    ap.add_argument("--n-burnin", type=int, default=N_BURNIN_DEFAULT,
                     help="N_BURNIN from mcmc_fit.py (used to compute the true "
                          "iteration number shown in the header, and burn-in "
                          "vs sampling phase)")
    ap.add_argument("--n-samples", type=int, default=N_SAMPLES_DEFAULT,
                     help="N_SAMPLES from mcmc_fit.py (used to compute the true "
                          "iteration number shown in the header)")
    ap.add_argument("--total-iterations", type=int, default=None,
                     help="Override: total MCMC iterations "
                          "(default: --n-burnin + --n-samples)")
    ap.add_argument("--no-header", dest="header", action="store_false", default=True,
                     help="Disable the top iteration-count header banner")
    ap.add_argument("--header-height", type=int, default=70,
                     help="Height in px of the top iteration-count banner")
    ap.add_argument("--output", default=os.path.join(HARVEST_ZIPS, "gibbs_fitting_grid.gif"))
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

    # ── Load every mouse's GIF ───────────────────────────────────────────
    per_mouse_frames: list[Optional[list[Image.Image]]] = []
    durations_ref: Optional[list[int]] = None
    tile_size: Optional[tuple[int, int]] = None

    for harvest in mice:
        gif_bytes = locate_gif_bytes(harvest)
        if gif_bytes is None:
            print(f"  MISSING: {harvest} — no {GIF_NAME} found, leaving slot blank")
            per_mouse_frames.append(None)
            continue

        frames, durations = load_frames(gif_bytes)
        print(f"  Loaded {harvest}: {len(frames)} frames")

        if tile_size is None:
            w0, h0 = frames[0].size
            tile_h = int(round(args.tile_width * h0 / w0))
            tile_size = (args.tile_width, tile_h)
        if durations_ref is None:
            durations_ref = durations

        # Resize + optionally label each frame
        resized = []
        font = make_label_font(max(12, args.tile_width // 24)) if args.label else None
        for fr in frames:
            fr = fr.resize(tile_size, Image.LANCZOS)
            if args.label:
                fr = fr.copy()
                draw = ImageDraw.Draw(fr)
                text = harvest
                pad = 4
                bbox = draw.textbbox((0, 0), text, font=font)
                tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
                draw.rectangle([0, 0, tw + 2 * pad, th + 2 * pad], fill=(0, 0, 0))
                draw.text((pad, pad), text, fill=(255, 255, 255), font=font)
            resized.append(fr)
        per_mouse_frames.append(resized)

    if tile_size is None:
        raise RuntimeError("No mouse GIFs could be found — nothing to compile.")

    # ── Sync frame counts (mice are fit for the same # iterations, so
    #     counts should already match; pad on the rare mismatch by holding
    #     the last available frame) ────────────────────────────────────────
    max_frames = max(len(f) for f in per_mouse_frames if f is not None)
    blank_tile = Image.new("RGB", tile_size, (30, 30, 30))
    for i, frames in enumerate(per_mouse_frames):
        if frames is None:
            per_mouse_frames[i] = [blank_tile] * max_frames
        elif len(frames) < max_frames:
            frames.extend([frames[-1]] * (max_frames - len(frames)))

    durations = durations_ref if durations_ref else [80] * max_frames
    if len(durations) < max_frames:
        durations = durations + [durations[-1]] * (max_frames - len(durations))

    # ── Map each animation frame back to its true MCMC iteration number ────
    # (reproduces the np.linspace(0, total_iters-1, N_ANIM_FRAMES) sampling
    # done when each per-mouse gif was generated)
    total_iters = args.total_iterations or (args.n_burnin + args.n_samples)
    if max_frames > 1:
        frame_iters = np.linspace(0, total_iters - 1, max_frames, dtype=int)
    else:
        frame_iters = np.array([total_iters - 1])

    header_font = None
    if args.header:
        header_font = make_label_font(max(16, int(args.header_height * 0.42)))

    # ── Compose the grid, frame by frame ────────────────────────────────
    canvas_w = args.cols * tile_size[0]
    grid_h   = args.rows * tile_size[1]
    header_h = args.header_height if args.header else 0
    canvas_h = grid_h + header_h

    grid_frames = []
    for t in range(max_frames):
        canvas = Image.new("RGB", (canvas_w, canvas_h), (0, 0, 0))
        for idx in range(n_slots):
            row, col = divmod(idx, args.cols)
            tile = (per_mouse_frames[idx][t] if idx < len(per_mouse_frames)
                     else blank_tile)
            canvas.paste(tile, (col * tile_size[0], header_h + row * tile_size[1]))

        if args.header:
            draw = ImageDraw.Draw(canvas)
            it_num = int(frame_iters[t]) + 1
            phase = "burn-in" if it_num <= args.n_burnin else "sampling"
            text = f"Iteration {it_num:,} / {total_iters:,}  ({phase})"
            bbox = draw.textbbox((0, 0), text, font=header_font)
            tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
            draw.text(((canvas_w - tw) / 2, (header_h - th) / 2 - bbox[1]),
                       text, fill=(255, 255, 255), font=header_font)

            # progress bar underneath the text
            bar_y0, bar_y1 = header_h - 8, header_h - 4
            draw.rectangle([0, bar_y0, canvas_w, bar_y1], fill=(60, 60, 60))
            frac = it_num / total_iters
            draw.rectangle([0, bar_y0, canvas_w * frac, bar_y1], fill=(70, 160, 235))

        grid_frames.append(canvas)
        if (t + 1) % 25 == 0 or t == max_frames - 1:
            print(f"  Composed frame {t + 1}/{max_frames}")

    # ── Save ──────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    grid_frames[0].save(
        args.output,
        save_all=True,
        append_images=grid_frames[1:],
        duration=durations[:max_frames],
        loop=0,
        optimize=True,
    )
    print(f"\nSaved compilation GIF: {args.output}  "
          f"({canvas_w}x{canvas_h}px, {max_frames} frames, "
          f"{args.rows}x{args.cols} grid"
          f"{', with iteration header' if args.header else ''})")


if __name__ == "__main__":
    main()