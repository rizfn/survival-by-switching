#!/usr/bin/env python3
"""Plot numerical and analytical invasion-growth curves versus switching period.

The C++ scan produces a finite-time numerical estimate of the invasion growth
rate for a grid of switching periods. The analytical asymptotic curve is
computed from the periodic-orbit Floquet exponent and overlaid for direct
comparison.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.optimize import root_scalar


def env_params():
    return {
        'A': {'r': 1.0, 'K': 1.0, 'gamma': 1.0},
        'B': {'r': 4.0, 'K': 2.0, 'gamma': 2.2},
    }


def phi(t, x, K, r):
    expo = np.exp(r * t)
    return (K * x * expo) / (K + x * (expo - 1.0))


def periodic_map(T, x, params):
    KA = params['A']['K']
    KB = params['B']['K']
    rA = params['A']['r']
    rB = params['B']['r']
    xA = phi(T / 2.0, x, KA, rA)
    return phi(T / 2.0, xA, KB, rB)


def periodic_fixed_point(T, params):
    KA = params['A']['K']
    KB = params['B']['K']

    def H(x):
        return periodic_map(T, x, params) - x

    lo = KA
    hi = KB
    flo = H(lo)
    fhi = H(hi)

    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        sample = np.linspace(max(1e-12, 0.5 * KA), 0.999999 * KB, 256)
        prev_x = float(sample[0])
        prev_f = H(prev_x)
        for x in sample[1:]:
            fx = H(float(x))
            if prev_f * fx < 0.0:
                sol = root_scalar(H, bracket=(prev_x, float(x)), method='brentq', xtol=1e-10)
                return float(sol.root) if sol.converged else np.nan
            prev_x = float(x)
            prev_f = fx
        return np.nan

    sol = root_scalar(H, bracket=(lo, hi), method='brentq', xtol=1e-10)
    return float(sol.root) if sol.converged else np.nan


def asymptotic_lambda(T, params):
    KA = params['A']['K']
    KB = params['B']['K']
    rA = params['A']['r']
    rB = params['B']['r']
    gammaA = params['A']['gamma']
    gammaB = params['B']['gamma']

    x0 = periodic_fixed_point(T, params)
    if not np.isfinite(x0):
        return np.nan
    x1 = phi(T / 2.0, x0, KA, rA)
    return 0.5 * (KA + KB - gammaA - gammaB) - (1.0 / T) * (
        (KA / rA) - (KB / rB)
    ) * np.log(x1 / x0)


def load_records(raw_dir: Path):
    rows = []
    for path in sorted(raw_dir.glob('task_*.tsv')):
        with path.open('r', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for row in reader:
                rows.append({
                    'x0': float(row['x0']),
                    'y0': float(row['y0']),
                    'T': float(row['T']),
                    'growth_rate': float(row['growth_rate']),
                    'valid': int(row['valid']),
                })
    return rows


def reduce_by_T(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[round(row['T'], 12)].append(row['growth_rate'])

    t_values = np.array(sorted(grouped.keys()), dtype=float)
    numeric = np.array([float(np.median(np.asarray(grouped[key], dtype=float))) for key in sorted(grouped.keys())], dtype=float)
    return t_values, numeric


def numerical_zero_crossing(t_values, growth_values):
    finite = np.isfinite(t_values) & np.isfinite(growth_values)
    t = t_values[finite]
    g = growth_values[finite]
    if t.size < 2:
        return np.nan

    order = np.argsort(t)
    t = t[order]
    g = g[order]

    for idx in range(t.size - 1):
        g0 = g[idx]
        g1 = g[idx + 1]
        if g0 == 0.0:
            return float(t[idx])
        if g0 * g1 < 0.0:
            t0 = t[idx]
            t1 = t[idx + 1]
            return float(t0 + (0.0 - g0) * (t1 - t0) / (g1 - g0))
    return np.nan


def plot_comparison(t_values, numeric_values, analytic_values, outdir: Path, tc_numeric: float, tc_analytic: float):
    plt.rcParams.update({
        'font.size': 24,
        'axes.labelsize': 24,
        'xtick.labelsize': 20,
        'ytick.labelsize': 20,
        'legend.fontsize': 24,
    })

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor('none')
    ax.set_facecolor('none')

    ax.axhline(0.0, color='black', linewidth=1.2, linestyle='--', alpha=0.7)
    ax.plot(t_values, analytic_values, color='#9E1528', linewidth=2.5, label=r'Analytical $\lambda(T)$')
    ax.scatter(t_values, numeric_values, color='#5299a8', s=30, alpha=0.85, label='Numerical data')

    if np.isfinite(tc_analytic):
        ax.axvline(tc_analytic, color='#666666', linestyle=':', linewidth=2.0, alpha=0.9, label=fr'Analytical $T_c \approx {tc_analytic:.3f}$')
    # if np.isfinite(tc_numeric):
    #     ax.axvline(tc_numeric, color='#5299a8', linestyle='-.', linewidth=2.0, alpha=0.9, label=fr'Numerical $T_c \approx {tc_numeric:.3f}$')

    ax.set_xscale('log')
    ax.set_xlabel('Switching period $T$')
    ax.set_ylabel('Invasion growth rate')
    # ax.set_title('Numerical vs analytical invasion growth rate')
    ax.grid(True, linestyle='--', alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='best')

    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / 'tc_comparison.pdf', dpi=250, facecolor='none', bbox_inches='tight')
    fig.savefig(outdir / 'tc_comparison.svg', dpi=250, facecolor='none', bbox_inches='tight')
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parent.parent
    raw_dir = root / 'outputs' / 'tc_scan' / 'raw'
    outdir = root / 'plots' / 'tc_scan_numeric'

    rows = load_records(raw_dir)
    if not rows:
        raise SystemExit(f'No raw scan files found in {raw_dir}')

    t_values, numeric_values = reduce_by_T(rows)
    params = env_params()
    analytic_values = np.array([asymptotic_lambda(float(T), params) for T in t_values], dtype=float)

    tc_numeric = numerical_zero_crossing(t_values, numeric_values)
    tc_analytic = numerical_zero_crossing(t_values, analytic_values)

    plot_comparison(t_values, numeric_values, analytic_values, outdir, tc_numeric, tc_analytic)
    print(f'Saved numerical/analytical comparison to {outdir}')


if __name__ == '__main__':
    main()
