#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${ROOT_DIR}/../.." && pwd)"
EXE="${ROOT_DIR}/simulate_tc_scan"
PYTHON_BIN="${PYTHON_BIN:-${REPO_ROOT}/.venv/bin/python}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
    PYTHON_BIN="$(command -v python3)"
fi

MAX_JOBS="${MAX_JOBS:-$(nproc)}"
STEPS_PER_PERIOD="${STEPS_PER_PERIOD:-600}"
TRANSIENT_PERIODS="${TRANSIENT_PERIODS:-12}"
MEASURE_PERIODS="${MEASURE_PERIODS:-6}"

T_MIN="${T_MIN:-0.2}"
T_MAX="${T_MAX:-12.0}"
NT="${NT:-48}"
X0="${X0:-0.2}"
Y0="${Y0:-1e-6}"

export X0 Y0 T_MIN T_MAX NT

g++ -std=c++17 -O3 -march=native -pipe -o "${EXE}" "${SCRIPT_DIR}/simulate_tc_scan.cpp"

TASK_FILE="$(mktemp)"
trap 'rm -f "${TASK_FILE}"' EXIT

"${PYTHON_BIN}" - "$TASK_FILE" <<'PY'
import numpy as np
import os
import sys
from pathlib import Path

task_path = Path(sys.argv[1])
T_MIN = float(os.environ['T_MIN'])
T_MAX = float(os.environ['T_MAX'])
NT = int(os.environ['NT'])
X0 = float(os.environ['X0'])
Y0 = float(os.environ['Y0'])

t_vals = np.geomspace(T_MIN, T_MAX, NT)

with task_path.open('w', encoding='utf-8') as handle:
    task_id = 0
    for T in t_vals:
        handle.write(f"{task_id} {X0:.16g} {Y0:.16g} {T:.16g}\n")
        task_id += 1
PY

running=0
while read -r task_id x0 y0 T; do
    "${EXE}" "${x0}" "${y0}" "${T}" "${task_id}" "${STEPS_PER_PERIOD}" "${TRANSIENT_PERIODS}" "${MEASURE_PERIODS}" &
    ((running += 1))
    if (( running >= MAX_JOBS )); then
        wait -n
        ((running -= 1))
    fi
done < "${TASK_FILE}"

wait

echo "Finished Tc scan. Raw outputs are in ${ROOT_DIR}/outputs/tc_scan/raw"
