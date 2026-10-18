# Existing-output phase and spectrum diagnostics

This follow-up reads existing certified fixed-oxygen responses only. It performs
zero model evaluations, generates zero new scientific designs, and submits zero
cluster jobs. N/M configurations and five phase repeats remain separate.

## Scope

* seed25: N4097/M4 versus N8193/M4; `p_mis_base`, `mu_hp`, `O2_crit`, and the stable-ST `buffer_beta` control; all 201 oxygen points.
* seed138: the same N comparison; five key parameters; focus 1.5–2.8% O2.
* seed101: N4097/M4 versus N4097/M8; five key parameters; focus 3–5% O2.
* seed174: the same M comparison; five key parameters; focus 0–1% O2.

Five key parameters are `buffer_beta`, `buffer_smax`, `p_mis_base`, `mu_hp`,
and `O2_crit`. There are 190 focal trajectories, each varying all 14 parameters.
Both outputs and the entire 0–5% grid at 0.025% spacing are parsed for every
trajectory, so comparisons outside focus windows are available as controls.
The cached ST screen separately retains all 14 parameters for all four seeds
and both settings; it does not claim that all 14 raw spectra were collected.

## Reproduce

Run from the repository root. Analysis requires NumPy and Matplotlib; the
collection inventory requires only Python's standard library. The final PDF
audit additionally uses pypdf; visual QA uses Poppler-rendered pages. The scientific
model is not imported. The inventory was made on HPC against the original
frozen code and completion receipts.

```bash
CODE=oxygen/code/O2_supply_demand_MAP/analysis/eFAST
ROOT=oxygen/results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8/phase_spectrum_diagnostics
mkdir -p "$ROOT/raw_cache"
rsync -a --files-from="$ROOT/rsync_files.txt" -e ssh \
  4482173@red.moffitt.org:/share/lab_crd/taoli/Project/soft_couping_eFAST/ \
  "$ROOT/raw_cache/"
python3 "$CODE/collect_phase_spectrum.py" verify-local --out "$ROOT" --cache "$ROOT/raw_cache"
python3 "$CODE/phase_spectrum.py" process --root "$ROOT"
python3 "$CODE/phase_spectrum.py" figures --root "$ROOT"
python3 "$CODE/report_phase_spectrum.py" --root "$ROOT"
python3 "$CODE/validate_phase_spectrum.py" --root "$ROOT"
```

`raw_cache/` (5.386 GB) and `processed_cache/` are ignored by Git. All input
identities, byte sizes, full SHA-256 hashes, compact figure source tables and
numerical checks are versioned. Rebuild the inventory only if intentionally
creating a new diagnostic snapshot; do not replace this snapshot's provenance.
The validation command records visual review as pending by default; pass
`--visual-render-reviewed` only after personally inspecting all six rendered PDF
pages. This delivered snapshot includes that completed review. Mathematical
fixtures are in `test_phase_spectrum.py` and can be run with `python3 -m unittest
-v test_phase_spectrum.py` from the code directory.

## Method and caveats

FFT is along the original ordered FAST trajectory separately at each oxygen
point. The original N/M determines the focal frequency, S1 harmonics and ST
cutoff. There is no smoothing, windowing, tail removal or oxygen-axis FFT.
Nonzero Fourier power is doubled for an odd-length real sequence. Parseval's
identity and the cached formal S1/ST and output variances are checked for every
trajectory and output. Stored input coordinates independently recover the
original shared phase and frequency assignments. The full raw files match
their certified receipts; original row flags were validated by the original
analyzer, while this local numeric parser reads the first five numeric columns.

Formal formulas and sample construction follow the primary implementation:
[SALib 1.5.2 sampler](https://github.com/SALib/SALib/blob/v1.5.2/src/SALib/sample/fast_sampler.py)
and [analyzer](https://github.com/SALib/SALib/blob/v1.5.2/src/SALib/analyze/fast.py).

An auxiliary cutoff sweep diagnoses how much power lies near the formal
partition. Its alternative partitions are **not** valid new-M analyses. A
neighborhood sensitivity index remains conditional on the prescribed ranges
and independent marginal distributions; it is not a posterior probability.
Phase is matched across settings, but sample IDs do not denote matched 14D
parameter vectors. Higher moments and phase correlations describe the paths;
they do not establish causation or calibrated significance with five repeats.

Representative figure O2 values are selected by the largest test-setting ST
phase range within each focus window **among the fixed 23 cached candidate
values** listed in the script. Exact full-grid maxima are also reported in
`band_diagnostics.tsv`. This selection is for diagnosis, not an unbiased
population summary. The four selected seeds cannot represent all 500 endpoints.

## Files

* `interpretation_zh.md`: detailed findings, figure reading and recommendations.
* `figures/diagnostic_overview.pdf`: distributions, V/L/H and ploidy-versus-growth ST.
* `figures/seed*_phase_spectrum.pdf`: all focused parameters; conditional joint-path response, Fourier CDF, formal ST across O2, raw spectra.
* `figures/phase_response_sequences.pdf`: original sequence responses and response-versus-gap scatter.
* `phase_metrics.tsv.gz`, `phase_spread_by_oxygen.tsv.gz`, `band_diagnostics.tsv`: complete quantitative diagnostic tables.
* `representative_phase_evidence.tsv`: phase-specific powers, tails, gaps and exact ST ratio decomposition.
* `representative_spectra.tsv.gz`, `representative_response_bins.tsv`, `representative_histograms.tsv.gz`, `raw_representative_responses.tsv.gz`: figure source data.
* `influential_response_points.tsv`: extreme points and their actual normalized 14D coordinates; points are retained in formal indices.
* `auxiliary_cutoff_sweep.tsv.gz`, `top_spectral_peaks.tsv.gz`: frequency-partition diagnostics.
* `input_design_diagnostics.tsv`, `all_parameter_cached_ST_screen.tsv`: stored input geometry and full cached-parameter screen.
* `inventory.json`, manifests, collection receipt, `numerical_validation.json`, `evidence_validation.json`, and final artifact audit: reproducibility and QA evidence.

`exact_low_power_term` and `exact_variance_term` use phase-mean powers Lbar/Vbar:

```
ST_r - (1-Lbar/Vbar)
  = -(L_r-Lbar)/V_r + Lbar*(V_r-Vbar)/(V_r*Vbar)
```

The terms can cancel. They are an algebraic decomposition, not independent
causal contributions. Historical source manifests and their files are rechecked
without altering previous results.
