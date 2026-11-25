"""Repository source and archived evidence locations for the joint workflow."""

from pathlib import Path

SOURCE = Path(__file__).resolve().parent
IN_VITRO = SOURCE.parents[1]
ROOT = next(parent for parent in IN_VITRO.parents if (parent / "oxygen").is_dir())
RESULTS = ROOT / "oxygen/results"
SYNTHETIC = RESULTS / "synthetic_identifiability"
STRUCTURAL = RESULTS / "structural_identifiability"
FROZEN_DESIGN = SYNTHETIC / "20261122_in_vitro_joint_structural_steps01_03_01"
ARCHIVED_RUN = SYNTHETIC / "20261122_in_vitro_joint_structural_steps04_07_01"
FIT = SYNTHETIC / "20261122_experimental_joint_fit_step02_01"
PROFILES = SYNTHETIC / "20261122_experimental_joint_identifiability_step03_01"
ALIGNMENT = SYNTHETIC / "20261119_experimental_karyotype_alignment_01"
FINITE_STATE = IN_VITRO / "02_structural_certificates/finite_state"


def output_directory(path):
    """Resolve an output path and protect sealed packages and the source tree."""
    path = Path(path).resolve()
    if ROOT.is_relative_to(path) or path.is_relative_to(ROOT / "oxygen/code") or path in (RESULTS, SYNTHETIC, STRUCTURAL):
        raise ValueError(f"Choose a dedicated run directory outside the source tree: {path}")
    archive = next((p for p in (path, *path.parents) if (p / "SHA256SUMS").is_file()), None)
    if archive is not None:
        raise ValueError(f"Choose an output directory outside the sealed package: {archive}")
    return path
