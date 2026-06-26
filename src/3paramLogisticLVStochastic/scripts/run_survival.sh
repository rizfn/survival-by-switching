#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${ROOT_DIR}/../.." && pwd)"
EXE="${ROOT_DIR}/simulate_survival"
PYTHON_BIN="${PYTHON_BIN:-${REPO_ROOT}/.venv/bin/python}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
    PYTHON_BIN="$(command -v python3)"
fi

# ── Switching rates to simulate ────────────────────────────────────────────────
# alpha_c ~ 1.5565 for the standard parameters; these straddle it.
ALPHAS=(
    1.24521
    1.55652
    1.86782
)

# ── Per-run knobs (override from the environment if you like) ───────────────────
T_MAX="${T_MAX:-6000}"
N_TRAJ="${N_TRAJ:-20000}"
DT="${DT:-0.01}"
Y_EXT="${Y_EXT:-1e-3}"
SEED="${SEED:-}"

g++ -std=c++17 -O3 -march=native -pipe -pthread -o "${EXE}" "${SCRIPT_DIR}/simulate_survival.cpp"

sim_id=0
for alpha in "${ALPHAS[@]}"; do
    "${EXE}" "${alpha}" "${sim_id}" "${T_MAX}" "${N_TRAJ}" "${DT}" "${Y_EXT}" ${SEED:+"${SEED}"}
    ((sim_id += 1))
done

echo "Finished survival runs. Raw outputs are in ${ROOT_DIR}/outputs/survivalProb"

# ── Plot ───────────────────────────────────────────────────────────────────────
"${PYTHON_BIN}" "${SCRIPT_DIR}/plot_survival_prob.py"
