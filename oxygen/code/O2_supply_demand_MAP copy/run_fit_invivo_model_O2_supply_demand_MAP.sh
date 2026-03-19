#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OXYGEN_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
FIT_SCRIPT="${SCRIPT_DIR}/fit_invivo_model_O2_supply_demand_MAP.R"
VIZ_SCRIPT="${SCRIPT_DIR}/viz_invivo_model_O2_supply_demand_MAP_results.R"

DEFAULT_CONFIG="${OXYGEN_DIR}/config/O2_supply_demand.yaml"
CONFIG_FILE="${DEFAULT_CONFIG}"

# Runner-level settings (YAML-overridable).
RUN_PREFIX=""
OUT_ROOT=""
DATA_DIR=""
SEEDS_FILE=""
SEEDS_CSV=""
APPEND_RUN_PREFIX_TIMESTAMP="FALSE"
RUN_PREFIX_TIMESTAMP_FORMAT="%Y%m%d_%H%M%S"
AUTO_VIZ="TRUE"
VIZ_REPORT_DT="1"
VIZ_TOP_N="6"

# Fit arguments loaded from YAML (and CLI passthrough).
FIT_CFG_ARGS=()
EXTRA_ARGS=()

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

append_fit_arg() {
  local key="$1"
  local val="$2"
  FIT_CFG_ARGS+=("--${key}=${val}")
}

has_fit_arg() {
  local key="$1"
  local a
  for a in "${FIT_CFG_ARGS[@]-}"; do
    [[ "$a" == --"${key}"=* ]] && return 0
  done
  return 1
}

get_fit_arg_value() {
  local key="$1"
  local a out=""
  for a in "${FIT_CFG_ARGS[@]-}"; do
    if [[ "$a" == --"${key}"=* ]]; then
      out="${a#*=}"
    fi
  done
  printf '%s' "$out"
}

peek_extra_arg_value() {
  local key="$1"
  local a out=""
  for a in "${EXTRA_ARGS[@]-}"; do
    if [[ "$a" == --"${key}"=* ]]; then
      out="${a#*=}"
    fi
  done
  printf '%s' "$out"
}

pop_extra_arg_value() {
  local key="$1"
  local a out=""
  local kept=()
  for a in "${EXTRA_ARGS[@]-}"; do
    if [[ "$a" == --"${key}"=* ]]; then
      out="${a#*=}"
    else
      kept+=("$a")
    fi
  done
  EXTRA_ARGS=("${kept[@]-}")
  printf '%s' "$out"
}

append_or_replace_fit_arg() {
  local key="$1"
  local val="$2"
  local a
  local kept=()
  for a in "${FIT_CFG_ARGS[@]-}"; do
    if [[ "$a" != --"${key}"=* ]]; then
      kept+=("$a")
    fi
  done
  FIT_CFG_ARGS=("${kept[@]-}")
  FIT_CFG_ARGS+=("--${key}=${val}")
}

remove_fit_arg_key() {
  local key="$1"
  local a
  local kept=()
  for a in "${FIT_CFG_ARGS[@]-}"; do
    if [[ "$a" != --"${key}"=* ]]; then
      kept+=("$a")
    fi
  done
  FIT_CFG_ARGS=("${kept[@]-}")
}

bool_true() {
  local v="$(echo "$1" | tr '[:upper:]' '[:lower:]')"
  [[ "$v" == "1" || "$v" == "true" || "$v" == "t" || "$v" == "yes" || "$v" == "y" ]]
}

log10_value() {
  local x="$1"
  awk -v x="$x" 'BEGIN { if (x+0<=0) { printf ""; exit 0 } printf "%.17g", log(x)/log(10) }'
}

lookup_param_table_natural() {
  local table_path="$1"
  local param_name="$2"
  local col_name="$3"
  [[ -n "$table_path" && -f "$table_path" ]] || { printf ''; return; }
  Rscript --vanilla - "$table_path" "$param_name" "$col_name" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
pname <- args[[2]]
col <- args[[3]]
if (!file.exists(path)) { cat(""); quit(save = "no", status = 0) }
tab <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
if (is.null(tab) || !all(c("param_name", "init", "lower", "upper") %in% names(tab))) {
  cat(""); quit(save = "no", status = 0)
}
i <- match(pname, tab$param_name)
if (is.na(i) || !(col %in% names(tab))) { cat(""); quit(save = "no", status = 0) }
tv <- suppressWarnings(as.numeric(tab[[col]][i]))
if (!is.finite(tv)) { cat(""); quit(save = "no", status = 0) }
nv <- 10^tv
if (!is.finite(nv) || nv <= 0) { cat(""); quit(save = "no", status = 0) }
cat(format(nv, scientific = TRUE, digits = 17, trim = TRUE))
RS
}

resolve_scalar_arg_value() {
  local key="$1"
  local yaml_group_key="$2"
  local allow_group="$3"
  local ptab_param="$4"
  local ptab_col="$5"
  local code_default="$6"
  local ptab_path="$7"

  local cli_val
  local yaml_base_val
  local yaml_group_val
  local ptab_val
  local resolved=""
  local source=""

  cli_val="$(pop_extra_arg_value "$key")"
  if [[ -n "$cli_val" ]]; then
    resolved="$cli_val"
    source="cli"
  else
    yaml_base_val="$(get_fit_arg_value "$key")"
    if [[ -n "$yaml_base_val" ]]; then
      resolved="$yaml_base_val"
      source="yaml_base"
    elif [[ "$allow_group" == "TRUE" && -n "$yaml_group_key" ]]; then
      yaml_group_val="$(get_fit_arg_value "$yaml_group_key")"
      if [[ -n "$yaml_group_val" ]]; then
        resolved="$yaml_group_val"
        source="yaml_group_default"
      fi
    fi
  fi

  if [[ -z "$resolved" && -n "$ptab_param" && -n "$ptab_col" ]]; then
    ptab_val="$(lookup_param_table_natural "$ptab_path" "$ptab_param" "$ptab_col")"
    if [[ -n "$ptab_val" ]]; then
      resolved="$ptab_val"
      source="parameter_table"
    fi
  fi

  if [[ -z "$resolved" ]]; then
    resolved="$code_default"
    source="code_default"
  fi

  append_or_replace_fit_arg "$key" "$resolved"
  append_or_replace_fit_arg "${key}_source" "$source"
}

resolve_death_args_precedence() {
  local ptab_path="$1"

  local death_cli death_yaml death_resolved death_source
  death_cli="$(pop_extra_arg_value "death")"
  if [[ -n "$death_cli" ]]; then
    death_resolved="$death_cli"
    death_source="cli"
  else
    death_yaml="$(get_fit_arg_value "death")"
    if [[ -n "$death_yaml" ]]; then
      death_resolved="$death_yaml"
      death_source="yaml_base"
    else
      death_resolved="TRUE"
      death_source="code_default"
    fi
  fi
  append_or_replace_fit_arg "death" "$death_resolved"
  append_or_replace_fit_arg "death_source" "$death_source"

  local use_group_defaults="FALSE"
  if bool_true "$death_resolved"; then
    use_group_defaults="TRUE"
  fi

  resolve_scalar_arg_value "mu_hp_init" "death_mu_hp_init" "$use_group_defaults" "log10_mu_hp" "init" "1e-3" "$ptab_path"
  resolve_scalar_arg_value "mu_hp_min" "death_mu_hp_min" "$use_group_defaults" "log10_mu_hp" "lower" "1e-8" "$ptab_path"
  resolve_scalar_arg_value "mu_hp_max" "death_mu_hp_max" "$use_group_defaults" "log10_mu_hp" "upper" "1.0" "$ptab_path"
  resolve_scalar_arg_value "k_clear_init" "death_k_clear_init" "$use_group_defaults" "log10_k_clear" "init" "1e-3" "$ptab_path"
  resolve_scalar_arg_value "k_clear_min" "death_k_clear_min" "$use_group_defaults" "log10_k_clear" "lower" "1e-8" "$ptab_path"
  resolve_scalar_arg_value "k_clear_max" "death_k_clear_max" "$use_group_defaults" "log10_k_clear" "upper" "1.0" "$ptab_path"

  local mu_hp_init_v k_clear_init_v prior_mu_center_default prior_k_center_default
  mu_hp_init_v="$(get_fit_arg_value "mu_hp_init")"
  k_clear_init_v="$(get_fit_arg_value "k_clear_init")"
  prior_mu_center_default="$(log10_value "$mu_hp_init_v")"
  prior_k_center_default="$(log10_value "$k_clear_init_v")"
  [[ -n "$prior_mu_center_default" ]] || prior_mu_center_default="-3"
  [[ -n "$prior_k_center_default" ]] || prior_k_center_default="-3"

  resolve_scalar_arg_value "prior_center_log10_mu_hp" "" "FALSE" "" "" "$prior_mu_center_default" ""
  resolve_scalar_arg_value "prior_sd_log10_mu_hp" "" "FALSE" "" "" "1.0" ""
  resolve_scalar_arg_value "prior_center_log10_k_clear" "" "FALSE" "" "" "$prior_k_center_default" ""
  resolve_scalar_arg_value "prior_sd_log10_k_clear" "" "FALSE" "" "" "1.0" ""

  append_or_replace_fit_arg "mu_hp_source" "$(get_fit_arg_value "mu_hp_init_source")"
  append_or_replace_fit_arg "k_clear_source" "$(get_fit_arg_value "k_clear_init_source")"
  append_or_replace_fit_arg "prior_mu_hp_source" "$(get_fit_arg_value "prior_center_log10_mu_hp_source")"
  append_or_replace_fit_arg "prior_k_clear_source" "$(get_fit_arg_value "prior_center_log10_k_clear_source")"

  # Keep resolved base keys in snapshot/CLI; drop raw death-group defaults to avoid ambiguity.
  remove_fit_arg_key "death_mu_hp_init"
  remove_fit_arg_key "death_mu_hp_min"
  remove_fit_arg_key "death_mu_hp_max"
  remove_fit_arg_key "death_k_clear_init"
  remove_fit_arg_key "death_k_clear_min"
  remove_fit_arg_key "death_k_clear_max"
}

resolve_parameter_table_path() {
  local a
  for (( idx=${#EXTRA_ARGS[@]}-1; idx>=0; idx-- )); do
    a="${EXTRA_ARGS[$idx]}"
    if [[ "$a" == --parameter_table=* ]]; then
      printf '%s' "${a#*=}"
      return
    fi
  done
  for (( idx=${#FIT_CFG_ARGS[@]}-1; idx>=0; idx-- )); do
    a="${FIT_CFG_ARGS[$idx]}"
    if [[ "$a" == --parameter_table=* ]]; then
      printf '%s' "${a#*=}"
      return
    fi
  done
  printf '%s' ""
}

yaml_escape() {
  local s="$1"
  s="${s//\\/\\\\}"
  s="${s//\"/\\\"}"
  printf '"%s"' "$s"
}

write_config_snapshots() {
  local run_dir="$1"
  local config_input="$2"
  local config_input_snapshot="${run_dir}/config.input.yaml"
  local config_resolved_snapshot="${run_dir}/config.resolved.yaml"
  local tmp_kv
  local tmp_resolved
  tmp_kv="$(mktemp)"
  tmp_resolved="$(mktemp)"

  cp "$config_input" "$config_input_snapshot"

  : > "$tmp_kv"
  for arg in "${FIT_CFG_ARGS[@]-}" "${EXTRA_ARGS[@]-}"; do
    case "$arg" in
      --*=*)
        local key="${arg%%=*}"
        local val="${arg#*=}"
        key="${key#--}"
        printf '%s\t%s\n' "$key" "$val" >> "$tmp_kv"
        ;;
    esac
  done

  if [[ -s "$tmp_kv" ]]; then
    awk -F '\t' '
      { idx[$1]=NR; key[NR]=$1; val[NR]=$2 }
      END {
        for (i=1; i<=NR; i++) {
          if (idx[key[i]] == i) print key[i] "\t" val[i]
        }
      }
    ' "$tmp_kv" > "$tmp_resolved"
  else
    : > "$tmp_resolved"
  fi

  {
    printf '# Resolved runtime config snapshot generated by runner\n'
    printf 'config_source: %s\n' "$(yaml_escape "$config_input")"
    printf 'run_prefix: %s\n' "$(yaml_escape "$RUN_PREFIX")"
    printf 'out_root: %s\n' "$(yaml_escape "$OUT_ROOT")"
    printf 'data_dir: %s\n' "$(yaml_escape "$DATA_DIR")"
    printf 'seeds_file: %s\n' "$(yaml_escape "${SEEDS_FILE}")"
    printf 'seeds_csv: %s\n' "$(yaml_escape "${SEEDS_CSV}")"
    printf 'seeds_use: %s\n' "$(yaml_escape "${SEEDS_USE}")"
    printf 'seed_source: %s\n' "$(yaml_escape "${SEED_SOURCE}")"
    printf 'append_run_prefix_timestamp: %s\n' "$(yaml_escape "${APPEND_RUN_PREFIX_TIMESTAMP}")"
    printf 'run_prefix_timestamp_format: %s\n' "$(yaml_escape "${RUN_PREFIX_TIMESTAMP_FORMAT}")"
    printf 'auto_viz: %s\n' "$(yaml_escape "${AUTO_VIZ}")"
    printf 'viz_report_dt: %s\n' "$(yaml_escape "${VIZ_REPORT_DT}")"
    printf 'viz_top_n: %s\n' "$(yaml_escape "${VIZ_TOP_N}")"
    printf 'run_dir: %s\n' "$(yaml_escape "${run_dir}")"
    printf 'fit_script: %s\n' "$(yaml_escape "${FIT_SCRIPT}")"
    printf 'viz_script: %s\n' "$(yaml_escape "${VIZ_SCRIPT}")"
    printf 'fit_args:\n'
    while IFS=$'\t' read -r k v || [[ -n "${k:-}" ]]; do
      [[ -z "${k:-}" ]] && continue
      printf '  %s: %s\n' "$k" "$(yaml_escape "$v")"
    done < "$tmp_resolved"
  } > "$config_resolved_snapshot"

  rm -f "$tmp_kv" "$tmp_resolved"
  printf '%s\t%s\n' "$config_input_snapshot" "$config_resolved_snapshot"
}

load_config() {
  local cfg_path="$1"
  local cfg_dir="$2"
  [[ -f "$cfg_path" ]] || { echo "ERROR: config not found: $cfg_path" >&2; exit 1; }

  while IFS= read -r raw || [[ -n "$raw" ]]; do
    local line="${raw%%#*}"
    line="$(trim "$line")"
    [[ -z "$line" ]] && continue
    [[ "$line" != *:* ]] && continue

    local key="${line%%:*}"
    local val="${line#*:}"
    key="$(trim "$key")"
    val="$(trim "$val")"

    if [[ "$val" == "\""*"\"" && "$val" == *"\"" ]]; then
      val="${val:1:${#val}-2}"
    elif [[ "$val" == "'"*"'" && "$val" == *"'" ]]; then
      val="${val:1:${#val}-2}"
    fi

    case "$key" in
      run_prefix) RUN_PREFIX="$val" ;;
      out_root) OUT_ROOT="$(resolve_path "$val" "$cfg_dir")" ;;
      data_dir) DATA_DIR="$(resolve_path "$val" "$cfg_dir")" ;;
      seeds_file) SEEDS_FILE="$(resolve_path "$val" "$cfg_dir")" ;;
      seeds_csv) SEEDS_CSV="$val" ;;
      append_run_prefix_timestamp) APPEND_RUN_PREFIX_TIMESTAMP="$val" ;;
      run_prefix_timestamp_format) RUN_PREFIX_TIMESTAMP_FORMAT="$val" ;;
      auto_viz) AUTO_VIZ="$val" ;;
      viz_report_dt) VIZ_REPORT_DT="$val" ;;
      viz_top_n) VIZ_TOP_N="$val" ;;
      parameter_table)
        append_fit_arg "$key" "$(resolve_path "$val" "$cfg_dir")"
        ;;
      *)
        append_fit_arg "$key" "$val"
        ;;
    esac
  done < "$cfg_path"
}

for arg in "$@"; do
  case "$arg" in
    --config=*) CONFIG_FILE="${arg#*=}" ;;
  esac
done

if [[ "$CONFIG_FILE" != /* ]]; then
  CONFIG_FILE="$(resolve_path "$CONFIG_FILE" "$PWD")"
fi
CONFIG_DIR="$(cd "$(dirname "$CONFIG_FILE")" && pwd)"
load_config "$CONFIG_FILE" "$CONFIG_DIR"

for arg in "$@"; do
  case "$arg" in
    --config=*) ;;
    --run_prefix=*) RUN_PREFIX="${arg#*=}" ;;
    --out_root=*) OUT_ROOT="$(resolve_path "${arg#*=}" "$PWD")" ;;
    --data_dir=*) DATA_DIR="$(resolve_path "${arg#*=}" "$PWD")" ;;
    --seeds_file=*) SEEDS_FILE="$(resolve_path "${arg#*=}" "$PWD")" ;;
    --seeds_csv=*) SEEDS_CSV="${arg#*=}" ;;
    --append_run_prefix_timestamp=*) APPEND_RUN_PREFIX_TIMESTAMP="${arg#*=}" ;;
    --run_prefix_timestamp_format=*) RUN_PREFIX_TIMESTAMP_FORMAT="${arg#*=}" ;;
    --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
    --viz_report_dt=*) VIZ_REPORT_DT="${arg#*=}" ;;
    --viz_top_n=*) VIZ_TOP_N="${arg#*=}" ;;
    --parameter_table=*)
      EXTRA_ARGS+=("--parameter_table=$(resolve_path "${arg#*=}" "$PWD")")
      ;;
    --*=*)
      EXTRA_ARGS+=("$arg")
      ;;
    *)
      EXTRA_ARGS+=("$arg")
      ;;
  esac
done

require_nonempty() {
  local name="$1"
  local value="$2"
  if [[ -z "${value}" ]]; then
    echo "ERROR: required value '${name}' is empty. Set it in YAML or CLI." >&2
    exit 1
  fi
}

require_nonempty "run_prefix" "${RUN_PREFIX}"
require_nonempty "out_root" "${OUT_ROOT}"
require_nonempty "data_dir" "${DATA_DIR}"
if [[ -z "${SEEDS_FILE}" && -z "${SEEDS_CSV}" ]]; then
  echo "ERROR: set seeds_file or seeds_csv in YAML/CLI." >&2
  exit 1
fi

if [[ "${APPEND_RUN_PREFIX_TIMESTAMP}" == "TRUE" || "${APPEND_RUN_PREFIX_TIMESTAMP}" == "true" || "${APPEND_RUN_PREFIX_TIMESTAMP}" == "1" ]]; then
  ts_suffix="$(date +"${RUN_PREFIX_TIMESTAMP_FORMAT}")"
  RUN_PREFIX="${RUN_PREFIX}_${ts_suffix}"
fi

read_seeds_from_file() {
  local f="$1"
  Rscript --vanilla - "$f" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
f <- args[[1]]
if (!file.exists(f)) { cat(""); quit(save = "no", status = 0) }
ln <- readLines(f, warn = FALSE)
ln <- gsub("\\r", "", ln)
ln <- ln[nzchar(trimws(ln))]
if (length(ln) == 0L) { cat(""); quit(save = "no", status = 0) }
raw <- paste(ln, collapse = ",")
parts <- trimws(unlist(strsplit(raw, "[,]", fixed = FALSE)))
parts <- parts[nzchar(parts)]
parts <- parts[!tolower(parts) %in% c("seed", "seeds")]
nums <- suppressWarnings(as.integer(parts))
nums <- nums[is.finite(nums)]
if (length(nums) == 0L) { cat(""); quit(save = "no", status = 0) }
cat(paste(nums, collapse = ","))
RS
}

SEEDS_FROM_FILE=""
if [[ -n "${SEEDS_FILE}" ]]; then
  SEEDS_FROM_FILE="$(read_seeds_from_file "${SEEDS_FILE}")"
fi
if [[ -n "${SEEDS_FROM_FILE}" ]]; then
  SEEDS_USE="${SEEDS_FROM_FILE}"
  SEED_SOURCE="file:${SEEDS_FILE}"
else
  SEEDS_USE="${SEEDS_CSV}"
  SEED_SOURCE="arg:--seeds_csv"
fi
if [[ -z "${SEEDS_USE}" ]]; then
  echo "ERROR: no seeds found. Provide ${SEEDS_FILE} or --seeds_csv." >&2
  exit 1
fi

# Resolve switch-gated death/prior arguments into a single effective set
# before writing snapshots and invoking fitter.
PARAM_TABLE_PATH_FOR_RESOLVE="$(resolve_parameter_table_path)"
resolve_death_args_precedence "${PARAM_TABLE_PATH_FOR_RESOLVE}"

mkdir -p "${OUT_ROOT}"
RUN_DIR="${OUT_ROOT}/${RUN_PREFIX}"
mkdir -p "${RUN_DIR}"
RUN_LOG="${RUN_DIR}/run_status.log"
SNAPSHOT_PATHS="$(write_config_snapshots "${RUN_DIR}" "${CONFIG_FILE}")"
CONFIG_INPUT_SNAPSHOT="${SNAPSHOT_PATHS%%$'\t'*}"
CONFIG_RESOLVED_SNAPSHOT="${SNAPSHOT_PATHS#*$'\t'}"
PARAM_TABLE_PATH="$(resolve_parameter_table_path)"
if [[ -n "${PARAM_TABLE_PATH}" && -f "${PARAM_TABLE_PATH}" ]]; then
  cp -f "${PARAM_TABLE_PATH}" "${RUN_DIR}/parameter_table.csv"
fi
touch "${RUN_LOG}"
exec > >(tee -a "${RUN_LOG}") 2>&1

echo "Running O2_supply_demand_MAP"
echo "  Config: ${CONFIG_FILE}"
echo "  Config input snapshot: ${CONFIG_INPUT_SNAPSHOT}"
echo "  Config resolved snapshot: ${CONFIG_RESOLVED_SNAPSHOT}"
echo "  Fit script: ${FIT_SCRIPT}"
echo "  Data dir: ${DATA_DIR}"
echo "  Seeds: ${SEEDS_USE} (${SEED_SOURCE})"
echo "  Run dir: ${RUN_DIR}"
echo "  Run log: ${RUN_LOG}"
if [[ -n "${PARAM_TABLE_PATH}" && -f "${PARAM_TABLE_PATH}" ]]; then
  echo "  Parameter table snapshot: ${RUN_DIR}/parameter_table.csv"
else
  echo "  Parameter table snapshot: (missing --parameter_table or file not found)"
fi
echo "  Run prefix timestamp suffix: ${APPEND_RUN_PREFIX_TIMESTAMP} (format=${RUN_PREFIX_TIMESTAMP_FORMAT})"
echo "  Auto viz: ${AUTO_VIZ} (report_dt=${VIZ_REPORT_DT}, top_n=${VIZ_TOP_N})"

IFS=',' read -r -a seed_arr <<< "${SEEDS_USE}"
for seed in "${seed_arr[@]}"; do
  seed="$(echo "$seed" | tr -d '[:space:]')"
  [[ -z "$seed" ]] && continue
  run_dir="${RUN_DIR}/seed${seed}"
  mkdir -p "${run_dir}"
  if [[ -n "${PARAM_TABLE_PATH}" && -f "${PARAM_TABLE_PATH}" ]]; then
    cp -f "${PARAM_TABLE_PATH}" "${run_dir}/parameter_table.csv"
  fi
  fit_log="${run_dir}/fit_status.log"
  viz_log="${run_dir}/viz_status.log"
  cmd=(
    Rscript "${FIT_SCRIPT}"
    "--seed=${seed}"
    "--out_dir=${run_dir}"
    "--data_dir=${DATA_DIR}"
  )
  if [[ ${#FIT_CFG_ARGS[@]} -gt 0 ]]; then
    cmd+=("${FIT_CFG_ARGS[@]}")
  fi
  if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    cmd+=("${EXTRA_ARGS[@]}")
  fi
  echo "[$(date '+%F %T')] seed=${seed}: start"
  echo "[$(date '+%F %T')] seed=${seed}: fit_log=${fit_log}"
  echo "Command: ${cmd[*]}"
  "${cmd[@]}" 2>&1 | tee "${fit_log}"
  echo "[$(date '+%F %T')] seed=${seed}: done"

  if [[ "${AUTO_VIZ}" == "TRUE" || "${AUTO_VIZ}" == "true" || "${AUTO_VIZ}" == "1" ]]; then
    viz_cmd=(
      Rscript "${VIZ_SCRIPT}"
      "--fit_dir=${run_dir}"
      "--data_dir=${DATA_DIR}"
      "--report_dt=${VIZ_REPORT_DT}"
      "--top_n=${VIZ_TOP_N}"
      "--n_cores=1"
    )
    echo "[$(date '+%F %T')] seed=${seed}: viz start"
    echo "[$(date '+%F %T')] seed=${seed}: viz_log=${viz_log}"
    echo "Viz command: ${viz_cmd[*]}"
    "${viz_cmd[@]}" 2>&1 | tee "${viz_log}"
    echo "[$(date '+%F %T')] seed=${seed}: viz done"
  fi
done

echo "All done. Run directory: ${RUN_DIR}"
