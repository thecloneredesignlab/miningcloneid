#!/usr/bin/env bash

# Deprecated compatibility loader. Canonical source-only helpers live in util/.

O2SD_PROV_COMPAT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
O2SD_SHELL_UTILS="$(cd "${O2SD_PROV_COMPAT_DIR}/../../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"
unset O2SD_PROV_COMPAT_DIR
