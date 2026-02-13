# Gemcitabine xenograft dataset: tumor burden time courses + terminal ploidy

This folder contains the **minimal data needed** to relate tumor-burden trajectories in mice treated with **gemcitabine** to **chromosome-length–weighted ploidy** measured at harvest.

Files included:

- `dt_Gem_VT_20260209_v5.xlsx` — tumor burden trajectories and metadata (wide format)
- `all_ploidy.tsv` — per-cell chromosome-length–weighted ploidy at harvest (long format)

> **Important:** Ploidy in `all_ploidy.tsv` is measured **only at harvest** (end of study). It should be interpreted as an *end-state readout* after both (i) early selection in the injection microenvironment and (ii) any treatment selection.

---

## Experimental design (conceptual)

### Cohorts / initial conditions (2N vs 4N)
Mice are injected with 1E6 of either:
- **2N** cells (near-diploid starting condition), or
- **4N** cells (near-tetraploid starting condition)

The cohort identity is encoded in the sample naming convention (typically visible in the `harvest` identifier; e.g., containing `2N` or `4N`).

**How this informs simulation initial conditions**
For simulations, the injection cohort provides the **initial ploidy state**:

- **2N cohort** → initialize the model with an initial distribution centered on **2N** (e.g., `N ≈ N_unit`, in the “pre-WGD” layer if your model distinguishes layers)
- **4N cohort** → initialize the model with an initial distribution centered on **4N** (e.g., `N ≈ 2*N_unit`, typically in the “post-WGD” layer if your model treats WGD explicitly)

### Treatment (gemcitabine)
Gemcitabine dosing is recorded per lesion/mouse in the tumor burden sheet:

- `Dose` encodes the assigned gemcitabine dose group.
  - `Dose = 0` represents **untreated / control**.
- `Day of 1st treatment` marks **treatment start** for that time series.

Tumor burden is measured at scheduled days (`Day_XX` columns) until the study endpoint. Once treatment has started each measurement day also corresponds to a treatment day. Missing values occur as **trailing NA only**, indicating measurement stops after sacrifice/endpoint.

---

## File: `dt_Gem_VT_20260209_v5.xlsx`

### Contents
This is the primary time-course table (one row per lesion/time series).

Typical c
