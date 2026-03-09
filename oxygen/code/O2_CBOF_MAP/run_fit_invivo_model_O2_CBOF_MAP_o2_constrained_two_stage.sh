#!/usr/bin/env bash
set -euo pipefail

# Two-stage constrained O2 fitting wrapper for O2_CBOF_MAP.
# Stage 1: constrain O2-decline terms and fix A_ang=0.
# Stage 2: keep constraints, release A_ang (bounded), warm-start each seed from Stage 1.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OXYGEN_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BASE_RUNNER="${SCRIPT_DIR}/run_fit_invivo_model_O2_CBOF_MAP.sh"
DEFAULT_CONFIG="${OXYGEN_DIR}/config/O2_CBOF.yaml"

CONFIG_FILE="${DEFAULT_CONFIG}"
RUN_PREFIX=""
OUT_ROOT=""
BASE_PARAM_TABLE=""
SEEDS_CSV=""
PASSTHROUGH_ARGS=()

# O2-constraint defaults
KDOWN_LOG10_MIN="8.0"
KDOWN_LOG10_MAX="9.3"
HDOWN_LOG10_MIN="0.0"
HDOWN_LOG10_MAX="0.45"
AANG_MAX="3.0"
M_ON_MIN="7.2"
M_ON_MAX="7.9"
DELTA_M_MIN="0.2"
DELTA_M_MAX="0.8"

STEP1_AUTO_VIZ="TRUE"
STEP2_AUTO_VIZ="TRUE"
STEP2_APPEND_TIMESTAMP="FALSE"

trim() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "$s"
}

resolve_path() {
  local p="$1"
  local base="$2"
  if [[ -z "$p" ]]; then
    printf '%s' "$p"
    return
  fi
  if [[ "$p" == /* ]]; then
    printf '%s' "$p"
    return
  fi
  if [[ "$p" == ~* ]]; then
    eval printf '%s' "$p"
    return
  fi
  printf '%s' "${base}/$p"
}

yaml_value() {
  local key="$1"
  local file="$2"
  awk -F: -v k="$key" '
    BEGIN { found=0 }
    /^[[:space:]]*#/ { next }
    {
      line=$0
      sub(/[[:space:]]*#.*/, "", line)
      if (line ~ "^[[:space:]]*"k"[[:space:]]*:") {
        sub("^[[:space:]]*"k"[[:space:]]*:[[:space:]]*", "", line)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
        gsub(/^"|"$/, "", line)
        gsub(/^'\''|'\''$/, "", line)
        print line
        found=1
        exit
      }
    }
    END { if (!found) exit 0 }
  ' "$file"
}

usage() {
  cat <<'USAGE'
Usage:
  bash run_fit_invivo_model_O2_CBOF_MAP_o2_constrained_two_stage.sh [--key=value ...]

Core:
  --config=PATH
  --run_prefix=NAME
  --out_root=DIR
  --parameter_table=PATH
  --seeds_csv=1,2,3

Constraint bounds:
  --kdown_log10_min=8.0      --kdown_log10_max=9.3
  --hdown_log10_min=0.0      --hdown_log10_max=0.45
  --aang_max=3.0
  --m_on_min=7.2             --m_on_max=7.9
  --delta_m_min=0.2          --delta_m_max=0.8

Viz switches:
  --step1_auto_viz=TRUE
  --step2_auto_viz=TRUE
  --step2_append_timestamp=FALSE

All unrecognized --key=value args are forwarded to the base runner.
USAGE
}

for arg in "$@"; do
  case "$arg" in
    --help|-h)
      usage
      exit 0
      ;;
    --config=*) CONFIG_FILE="${arg#*=}" ;;
    --run_prefix=*) RUN_PREFIX="${arg#*=}" ;;
    --out_root=*) OUT_ROOT="${arg#*=}" ;;
    --parameter_table=*) BASE_PARAM_TABLE="${arg#*=}" ;;
    --seeds_csv=*) SEEDS_CSV="${arg#*=}" ;;
    --kdown_log10_min=*) KDOWN_LOG10_MIN="${arg#*=}" ;;
    --kdown_log10_max=*) KDOWN_LOG10_MAX="${arg#*=}" ;;
    --hdown_log10_min=*) HDOWN_LOG10_MIN="${arg#*=}" ;;
    --hdown_log10_max=*) HDOWN_LOG10_MAX="${arg#*=}" ;;
    --aang_max=*) AANG_MAX="${arg#*=}" ;;
    --m_on_min=*) M_ON_MIN="${arg#*=}" ;;
    --m_on_max=*) M_ON_MAX="${arg#*=}" ;;
    --delta_m_min=*) DELTA_M_MIN="${arg#*=}" ;;
    --delta_m_max=*) DELTA_M_MAX="${arg#*=}" ;;
    --step1_auto_viz=*) STEP1_AUTO_VIZ="${arg#*=}" ;;
    --step2_auto_viz=*) STEP2_AUTO_VIZ="${arg#*=}" ;;
    --step2_append_timestamp=*) STEP2_APPEND_TIMESTAMP="${arg#*=}" ;;
    *)
      PASSTHROUGH_ARGS+=("$arg")
      ;;
  esac
done

[[ "$CONFIG_FILE" == /* ]] || CONFIG_FILE="$(resolve_path "$CONFIG_FILE" "$PWD")"
CONFIG_DIR="$(cd "$(dirname "$CONFIG_FILE")" && pwd)"
[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config not found: $CONFIG_FILE" >&2; exit 1; }
[[ -f "$BASE_RUNNER" ]] || { echo "ERROR: base runner not found: $BASE_RUNNER" >&2; exit 1; }

if [[ -z "$RUN_PREFIX" ]]; then
  RUN_PREFIX="$(yaml_value run_prefix "$CONFIG_FILE")"
fi
if [[ -z "$RUN_PREFIX" ]]; then
  RUN_PREFIX="fit_invivo_model_O2_CBOF_MAP"
fi
RUN_PREFIX="${RUN_PREFIX}_o2constrained2step"

if [[ -z "$OUT_ROOT" ]]; then
  OUT_ROOT="$(yaml_value out_root "$CONFIG_FILE")"
fi
OUT_ROOT="$(resolve_path "${OUT_ROOT:-../results}" "$CONFIG_DIR")"

if [[ -z "$BASE_PARAM_TABLE" ]]; then
  BASE_PARAM_TABLE="$(yaml_value parameter_table "$CONFIG_FILE")"
fi
BASE_PARAM_TABLE="$(resolve_path "${BASE_PARAM_TABLE:-../data/O2_CBOF/parameter_table.csv}" "$CONFIG_DIR")"
[[ -f "$BASE_PARAM_TABLE" ]] || { echo "ERROR: parameter_table not found: $BASE_PARAM_TABLE" >&2; exit 1; }

TMP_DIR="$(mktemp -d)"
STEP1_TABLE="${TMP_DIR}/parameter_table_step1.csv"
STEP2_TABLE="${TMP_DIR}/parameter_table_step2.csv"

python3 - "$BASE_PARAM_TABLE" "$STEP1_TABLE" "$STEP2_TABLE" \
  "$KDOWN_LOG10_MIN" "$KDOWN_LOG10_MAX" \
  "$HDOWN_LOG10_MIN" "$HDOWN_LOG10_MAX" \
  "$AANG_MAX" "$M_ON_MIN" "$M_ON_MAX" "$DELTA_M_MIN" "$DELTA_M_MAX" <<'PY'
import csv
import math
import sys

base, step1, step2 = sys.argv[1], sys.argv[2], sys.argv[3]
kdown_min, kdown_max = float(sys.argv[4]), float(sys.argv[5])
hdown_min, hdown_max = float(sys.argv[6]), float(sys.argv[7])
aang_max = float(sys.argv[8])
m_on_min, m_on_max = float(sys.argv[9]), float(sys.argv[10])
delta_m_min, delta_m_max = float(sys.argv[11]), float(sys.argv[12])

if delta_m_min <= 0 or delta_m_max <= 0:
    raise SystemExit("delta_m_min/max must be > 0")
if delta_m_min > delta_m_max:
    raise SystemExit("delta_m_min must be <= delta_m_max")

ldm_min = math.log10(delta_m_min)
ldm_max = math.log10(delta_m_max)

with open(base, newline="") as f:
    rows = list(csv.DictReader(f))
    fieldnames = list(rows[0].keys()) if rows else [
        "param_name","estimate","init_value","lower_bound","upper_bound","source","note"
    ]

def _clip(v, lo, hi):
    return max(lo, min(hi, v))

def apply_constraints(rows, stage):
    out = []
    for r in rows:
        rr = dict(r)
        name = rr.get("param_name", "").strip()
        if not name:
            out.append(rr)
            continue

        if name == "log10_K_down":
            rr["lower_bound"] = f"{kdown_min:.6f}"
            rr["upper_bound"] = f"{kdown_max:.6f}"
            rr["init_value"] = f"{_clip(float(rr['init_value']), kdown_min, kdown_max):.6f}"
        elif name == "log10_h_down":
            rr["lower_bound"] = f"{hdown_min:.6f}"
            rr["upper_bound"] = f"{hdown_max:.6f}"
            rr["init_value"] = f"{_clip(float(rr['init_value']), hdown_min, hdown_max):.6f}"
        elif name == "A_ang":
            if stage == 1:
                rr["lower_bound"] = "0.000000"
                rr["upper_bound"] = "0.000000"
                rr["init_value"] = "0.000000"
            else:
                rr["lower_bound"] = "0.000000"
                rr["upper_bound"] = f"{aang_max:.6f}"
                rr["init_value"] = f"{_clip(float(rr['init_value']), 0.0, aang_max):.6f}"
        elif name == "m_on":
            rr["lower_bound"] = f"{m_on_min:.6f}"
            rr["upper_bound"] = f"{m_on_max:.6f}"
            rr["init_value"] = f"{_clip(float(rr['init_value']), m_on_min, m_on_max):.6f}"
        elif name == "log10_delta_m":
            rr["lower_bound"] = f"{ldm_min:.6f}"
            rr["upper_bound"] = f"{ldm_max:.6f}"
            rr["init_value"] = f"{_clip(float(rr['init_value']), ldm_min, ldm_max):.6f}"
        out.append(rr)
    return out

rows1 = apply_constraints(rows, stage=1)
rows2 = apply_constraints(rows, stage=2)

for path, rows_use in [(step1, rows1), (step2, rows2)]:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows_use)
PY

mkdir -p "$OUT_ROOT"
STEP1_PREFIX="${RUN_PREFIX}_step1"
STEP2_PREFIX="${RUN_PREFIX}_step2"
STEP1_DIR="${OUT_ROOT}/${STEP1_PREFIX}"
STEP2_DIR="${OUT_ROOT}/${STEP2_PREFIX}"

echo "Two-stage constrained O2 fit (O2_CBOF_MAP)"
echo "  config: ${CONFIG_FILE}"
echo "  out_root: ${OUT_ROOT}"
echo "  step1_table: ${STEP1_TABLE}"
echo "  step2_table: ${STEP2_TABLE}"
echo "  constraints:"
echo "    log10_K_down in [${KDOWN_LOG10_MIN}, ${KDOWN_LOG10_MAX}]"
echo "    log10_h_down in [${HDOWN_LOG10_MIN}, ${HDOWN_LOG10_MAX}]"
echo "    A_ang stage1 fixed 0; stage2 in [0, ${AANG_MAX}]"
echo "    m_on in [${M_ON_MIN}, ${M_ON_MAX}]"
echo "    delta_m in [${DELTA_M_MIN}, ${DELTA_M_MAX}]"

STEP1_CMD=(
  bash "$BASE_RUNNER"
  "--config=${CONFIG_FILE}"
  "--run_prefix=${STEP1_PREFIX}"
  "--out_root=${OUT_ROOT}"
  "--append_run_prefix_timestamp=FALSE"
  "--parameter_table=${STEP1_TABLE}"
  "--auto_viz=${STEP1_AUTO_VIZ}"
)
if [[ -n "$SEEDS_CSV" ]]; then
  STEP1_CMD+=("--seeds_file=/dev/null" "--seeds_csv=${SEEDS_CSV}")
fi
if [[ ${#PASSTHROUGH_ARGS[@]} -gt 0 ]]; then
  STEP1_CMD+=("${PASSTHROUGH_ARGS[@]}")
fi

echo
echo "[Stage 1] Running base fit with A_ang fixed to 0 ..."
"${STEP1_CMD[@]}"

[[ -d "$STEP1_DIR" ]] || { echo "ERROR: Stage1 output dir not found: $STEP1_DIR" >&2; exit 1; }

seed_dirs=( "$STEP1_DIR"/seed* )
if [[ ! -d "${seed_dirs[0]}" ]]; then
  echo "ERROR: no seed directories found under stage1 output: $STEP1_DIR" >&2
  exit 1
fi

echo
echo "[Stage 2] Warm-start per seed from Stage 1 ..."
for sd in "${seed_dirs[@]}"; do
  [[ -d "$sd" ]] || continue
  seed="$(basename "$sd" | sed 's/^seed//')"
  init_tsv="${sd}/fit_parameter_stages.tsv"
  [[ -f "$init_tsv" ]] || { echo "ERROR: missing warm-start file: $init_tsv" >&2; exit 1; }

  STEP2_CMD=(
    bash "$BASE_RUNNER"
    "--config=${CONFIG_FILE}"
    "--run_prefix=${STEP2_PREFIX}"
    "--out_root=${OUT_ROOT}"
    "--append_run_prefix_timestamp=${STEP2_APPEND_TIMESTAMP}"
    "--parameter_table=${STEP2_TABLE}"
    "--seeds_file=/dev/null"
    "--seeds_csv=${seed}"
    "--init_params_tsv=${init_tsv}"
    "--auto_viz=${STEP2_AUTO_VIZ}"
  )
  if [[ ${#PASSTHROUGH_ARGS[@]} -gt 0 ]]; then
    STEP2_CMD+=("${PASSTHROUGH_ARGS[@]}")
  fi
  echo "  - seed=${seed}"
  "${STEP2_CMD[@]}"
done

echo
echo "Done."
echo "  Stage1: ${STEP1_DIR}"
echo "  Stage2: ${STEP2_DIR}"
echo "  Temp files: ${TMP_DIR}"
