#!/usr/bin/env python3
"""Plot the analytical invasion-growth curve versus switching period.

This is the exact asymptotic Floquet exponent derived from the periodic
boundary orbit. The figure is meant to be visually comparable to the numerical
simulation plot.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar


def env_params():
    return {
        'A': {'r': 1.0, 'K': 1.0, 'gamma': 1.0},
        'B': {'r': 4.0, 'K': 2.0, 'gamma': 2.2},
    }


def phi(t, x, K, r):
    expo = math.exp(r * t)
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
    ) * math.log(x1 / x0)


def numerical_zero_crossing(t_values, g_values):
    finite = np.isfinite(t_values) & np.isfinite(g_values)
    t = t_values[finite]
    g = g_values[finite]
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


def make_curve():
    params = env_params()
    t_values = np.geomspace(1e-3, 20.0, 280)
    growth = np.array([asymptotic_lambda(float(T), params) for T in t_values], dtype=float)
    return t_values, growth, params


def plot_curve(t_values, growth, outdir: Path, tc_value: float):
    plt.rcParams.update({
        'font.size': 18,
        'axes.titlesize': 22,
        'axes.labelsize': 20,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'legend.fontsize': 16,
    })

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor('none')
    ax.set_facecolor('none')

    ax.axhline(0.0, color='black', linewidth=1.2, linestyle='--', alpha=0.7)
    ax.plot(t_values, growth, color='tab:blue', linewidth=2.6, label=r'Analytical $\lambda(T)$')
    if np.isfinite(tc_value):
        ax.axvline(tc_value, color='tab:blue', linestyle=':', linewidth=2.0, label=fr'$T_c \approx {tc_value:.3f}$')

    ax.set_xscale('log')
    ax.set_xlabel('Switching period $T$')
    ax.set_ylabel('Invasion growth rate')
    ax.set_title('Analytical asymptotic invasion growth rate')
    ax.grid(True, linestyle='--', alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='best')

    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / 'tc_analytic_curve.pdf', dpi=250, facecolor='none', bbox_inches='tight')
    fig.savefig(outdir / 'tc_analytic_curve.svg', dpi=250, facecolor='none', bbox_inches='tight')
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parent.parent
    outdir = root / 'plots' / 'tc_scan_analytic'
    t_values, growth, params = make_curve()
    tc_value = numerical_zero_crossing(t_values, growth)
    plot_curve(t_values, growth, outdir, tc_value)
    print(f'Saved analytical curve to {outdir}')


if __name__ == '__main__':
    main()
