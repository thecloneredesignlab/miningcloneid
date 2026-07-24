# Immutable Raster Policy

## Purpose

Manuscript figures should normally be regenerated from scripts and direct
inputs. This policy defines the narrow exception for a source image whose pixels
must be embedded unchanged, such as a microscopy field, scanned specimen,
instrument screenshot required as evidence, or externally produced image for
which no underlying renderable source exists.

The repository has exactly one canonical exception directory:

`user-approved-raster-figures/`

The directory is an input store, not an output directory.

## What qualifies

An asset may be added only when all of the following are true:

1. The panel cannot reasonably be regenerated from code, tabular data, vector
   source, or the original scientific image source.
2. The user has explicitly approved the exact displayed image for manuscript
   use.
3. The original bytes are retained without cropping, recoloring, resampling,
   compression, or other pixel modification.
4. The asset is entered in `user-approved-raster-figures/manifest.tsv` with its
   SHA-256 checksum and approval reference.
5. The raster validator passes before figure assembly.

Generated model plots, extracted HTML report plots, screenshots of regenerable
figures, and montage images do not qualify. Missing plotting code is a
provenance problem to resolve, not an automatic raster exception.

## Approval and intake

1. Inspect the candidate and its scientific source.
2. Record the user approval in a durable issue, pull request, commit message, or
   dated project note.
3. Copy the approved file into the canonical directory without modifying it.
4. Calculate its SHA-256 checksum.
5. Add one manifest row with a stable asset ID, relative path, byte count,
   approval reference, intended panel, source description, and notes.
6. Run `python3 scripts/validate_approved_rasters.py`.

If an image must be altered, make the alteration upstream, obtain approval for
the new exact image, and register it as a new asset. Never overwrite an approved
asset while retaining its old manifest identity.

## Use in figure code

- Read the registered file directly from the canonical directory.
- Do not crop, trim, recolor, resample, rewrite, or replace its pixels.
- Record the asset ID and checksum in the figure rebuild manifest.
- Layout scaling is allowed only as a non-destructive display operation that
  preserves aspect ratio; the source file must remain unchanged.
- If the intended panel cannot use the image without pixel editing, stop and
  request approval for a newly prepared source asset.

## Storage

The manifest and policy are tracked in normal Git. Large raster assets should
use a path-specific Git LFS rule before they are committed. Do not add a global
LFS rule that silently changes the handling of existing repository PNG files.

