#!/usr/bin/env python3
"""Plot y_max and y_mean on the limit cycle vs switching period T.

For each of several duty cycles p, simulate the 2D ODE across log-spaced
T values, plus a T→0 point (averaged ODE).

Uses numba @njit for the inner RK4 loops and multiprocessing across T values.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool, cpu_count
from numba import njit


# ── Numba-compiled RK4 for switching ODE ─────────────────────────────────────

@njit
def _rhs_switch(x, logy, t, T, tau_A, rA, KA, gA, rB, KB, gB):
    N     = len(x)
    dx    = np.empty(N)
    dlogy = np.empty(N)
    phase = t % T
    for j in range(N):
        in_A = phase < tau_A[j]
        r = rA if in_A else rB
        K = KA if in_A else KB
        g = gA if in_A else gB
        y = np.exp(min(max(logy[j], -500.0), 500.0))
        dx[j]    = r * x[j] * (1.0 - x[j] / K) - x[j] * y
        dlogy[j] = x[j] - g
    return dx, dlogy


@njit
def _rk4_switch(x, logy, t, T, tau_A, rA, KA, gA, rB, KB, gB, dt):
    k1x, k1l = _rhs_switch(x,           logy,           t,       T, tau_A, rA, KA, gA, rB, KB, gB)
    k2x, k2l = _rhs_switch(x+.5*dt*k1x, logy+.5*dt*k1l, t+.5*dt, T, tau_A, rA, KA, gA, rB, KB, gB)
    k3x, k3l = _rhs_switch(x+.5*dt*k2x, logy+.5*dt*k2l, t+.5*dt, T, tau_A, rA, KA, gA, rB, KB, gB)
    k4x, k4l = _rhs_switch(x+   dt*k3x, logy+   dt*k3l, t+dt,    T, tau_A, rA, KA, gA, rB, KB, gB)
    x_new    = np.empty_like(x)
    logy_new = np.empty_like(logy)
    for j in range(len(x)):
        x_new[j]    = max(x[j]    + (dt/6.0)*(k1x[j]+2*k2x[j]+2*k3x[j]+k4x[j]), 1e-14)
        logy_new[j] =     logy[j] + (dt/6.0)*(k1l[j]+2*k2l[j]+2*k3l[j]+k4l[j])
    return x_new, logy_new


@njit
def _simulate_T_core(T, p_arr, tau_A, steps_per_period,
                     transient_periods, measure_periods,
                     rA, KA, gA, rB, KB, gB):
    N  = len(p_arr)
    dt = T / steps_per_period

    x    = np.full(N, 0.8 * KA)
    logy = np.full(N, np.log(0.01))
    t    = 0.0

    for _ in range(transient_periods * steps_per_period):
        x, logy = _rk4_switch(x, logy, t, T, tau_A, rA, KA, gA, rB, KB, gB, dt)
        t += dt

    for j in range(N):
        logy[j] = max(logy[j], np.log(1e-14))

    y_max  = np.full(N, np.nan)
    y_sum  = np.zeros(N)
    n_meas = measure_periods * steps_per_period

    for _ in range(n_meas):
        phase      = t % T
        next_phase = (t + dt) % T
        y_now      = np.empty(N)
        for j in range(N):
            y_now[j] = np.exp(min(max(logy[j], -500.0), 500.0))
            crossing = (phase < tau_A[j]) and (
                (next_phase >= tau_A[j]) or (next_phase < phase))
            if crossing:
                y_max[j] = y_now[j]
            y_sum[j] += y_now[j]
        x, logy = _rk4_switch(x, logy, t, T, tau_A, rA, KA, gA, rB, KB, gB, dt)
        t += dt

    y_mean = y_sum / n_meas
    for j in range(N):
        if np.exp(min(max(logy[j], -500.0), 500.0)) < 1e-10:
            y_max[j]  = np.nan
            y_mean[j] = np.nan

    return y_max, y_mean


# ── Numba-compiled averaged ODE ───────────────────────────────────────────────

@njit
def _avg_stats_core(p_arr, dt, transient_steps, measure_steps,
                    rA, KA, gA, rB, KB, gB):
    N   = len(p_arr)
    r_b = np.empty(N); g_b = np.empty(N); c_b = np.empty(N)
    for j in range(N):
        r_b[j] = p_arr[j]*rA + (1-p_arr[j])*rB
        g_b[j] = p_arr[j]*gA + (1-p_arr[j])*gB
        c_b[j] = p_arr[j]*(rA/KA) + (1-p_arr[j])*(rB/KB)

    x    = np.full(N, 0.8*KA)
    logy = np.full(N, np.log(0.01))

    def step(x, logy):
        k1x = np.empty(N); k1l = np.empty(N)
        for j in range(N):
            y = np.exp(min(max(logy[j], -500.0), 500.0))
            k1x[j] = r_b[j]*x[j] - c_b[j]*x[j]**2 - x[j]*y
            k1l[j] = x[j] - g_b[j]
        k2x = np.empty(N); k2l = np.empty(N)
        for j in range(N):
            y  = np.exp(min(max(logy[j]+.5*dt*k1l[j], -500.0), 500.0))
            xj = x[j]+.5*dt*k1x[j]
            k2x[j] = r_b[j]*xj - c_b[j]*xj**2 - xj*y
            k2l[j] = xj - g_b[j]
        k3x = np.empty(N); k3l = np.empty(N)
        for j in range(N):
            y  = np.exp(min(max(logy[j]+.5*dt*k2l[j], -500.0), 500.0))
            xj = x[j]+.5*dt*k2x[j]
            k3x[j] = r_b[j]*xj - c_b[j]*xj**2 - xj*y
            k3l[j] = xj - g_b[j]
        k4x = np.empty(N); k4l = np.empty(N)
        for j in range(N):
            y  = np.exp(min(max(logy[j]+dt*k3l[j], -500.0), 500.0))
            xj = x[j]+dt*k3x[j]
            k4x[j] = r_b[j]*xj - c_b[j]*xj**2 - xj*y
            k4l[j] = xj - g_b[j]
        x_new = np.empty(N); logy_new = np.empty(N)
        for j in range(N):
            x_new[j]    = max(x[j]    + (dt/6)*(k1x[j]+2*k2x[j]+2*k3x[j]+k4x[j]), 1e-14)
            logy_new[j] =     logy[j] + (dt/6)*(k1l[j]+2*k2l[j]+2*k3l[j]+k4l[j])
        return x_new, logy_new

    for _ in range(transient_steps):
        x, logy = step(x, logy)

    y_max_acc = np.zeros(N)
    y_sum_acc = np.zeros(N)
    for _ in range(measure_steps):
        x, logy = step(x, logy)
        for j in range(N):
            y = np.exp(min(max(logy[j], -500.0), 500.0))
            if y > y_max_acc[j]:
                y_max_acc[j] = y
            y_sum_acc[j] += y

    y_mean_acc = y_sum_acc / measure_steps
    for j in range(N):
        if np.exp(min(max(logy[j], -500.0), 500.0)) < 1e-10:
            y_max_acc[j]  = np.nan
            y_mean_acc[j] = np.nan

    return y_max_acc, y_mean_acc


# ── Python wrappers (picklable for multiprocessing) ───────────────────────────

def simulate_T(args):
    T, p_values, steps_per_period, transient_periods, measure_periods, params = args
    rA, KA, gA, rB, KB, gB = params
    p_arr = np.asarray(p_values)
    tau_A = p_arr * T
    # Ensure at least MIN_TRANSIENT_TIME physical time units regardless of T.
    # For small T, transient_periods would only cover transient_periods*T time,
    # which can be far too short to converge from the initial condition.
    MIN_TRANSIENT_TIME = 500.0
    effective_transient = max(transient_periods, int(np.ceil(MIN_TRANSIENT_TIME / T)))
    y_max, y_mean = _simulate_T_core(
        T, p_arr, tau_A, steps_per_period, effective_transient, measure_periods,
        rA, KA, gA, rB, KB, gB)
    print(f"  T={T:.2e} (transient_periods={effective_transient}) done", flush=True)
    return T, y_max, y_mean


def avg_stats(p_values, rA, KA, gA, rB, KB, gB,
              dt=1e-5, transient_time=300.0, measure_time=100.0):
    return _avg_stats_core(
        np.asarray(p_values), dt,
        int(transient_time / dt), int(measure_time / dt),
        rA, KA, gA, rB, KB, gB)


# ── Main ──────────────────────────────────────────────────────────────────────

def main(rA, KA, gA, rB, KB, gB, p_values, T_scan,
         steps_per_period, transient_periods, measure_periods, outdir):

    os.makedirs(outdir, exist_ok=True)
    params = (rA, KA, gA, rB, KB, gB)
    p_arr  = np.asarray(p_values)

    print("Warming up numba JIT...")
    _simulate_T_core(1.0, p_arr, p_arr.copy(), 2, 2, 2, *params)
    _avg_stats_core(p_arr, 0.1, 5, 5, *params)
    print("  done")

    print("Computing T→0 (averaged ODE)...")
    avg_ymax, avg_ymean = avg_stats(p_arr, *params)
    print("  done")

    args_list = [
        (T, p_arr, steps_per_period, transient_periods, measure_periods, params)
        for T in T_scan
    ]
    n_workers = max(1, cpu_count() - 1)
    print(f"Computing {len(T_scan)} T values using {n_workers} workers...")
    with Pool(n_workers) as pool:
        results = pool.map(simulate_T, args_list)

    results.sort(key=lambda r: r[0])
    T_arr     = np.array([r[0] for r in results])
    ymax_mat  = np.array([r[1] for r in results])   # (N_T, N_p)
    ymean_mat = np.array([r[2] for r in results])

    # ── Plot ──────────────────────────────────────────────────────────────────
    plt.rcParams.update({
        'font.size': 20, 'axes.labelsize': 20,
        'xtick.labelsize': 16, 'ytick.labelsize': 16, 'legend.fontsize': 16,
    })
    cmap   = plt.get_cmap('cool_r')
    p_norm = (p_arr - p_arr.min()) / (p_arr.max() - p_arr.min())
    colors = [cmap(n) for n in p_norm]

    fname_base = f'r_{rA}-{rB}_K_{KA}-{KB}_g_{gA}-{gB}'

    for mat, avg_vals, ylabel, suffix in [
        (ymax_mat,  avg_ymax,  r'$y_{\max}$ on limit cycle',          'ymax_vs_T'),
        (ymean_mat, avg_ymean, r'$\langle y \rangle$ on limit cycle',  'ymean_vs_T'),
    ]:
        fig, ax = plt.subplots(figsize=(9, 6))
        fig.patch.set_facecolor('none')
        ax.set_facecolor('none')

        for j, p in enumerate(p_arr):
            ax.plot(T_arr, mat[:, j], color=colors[j], linewidth=2.0,
                    label=f'$p={p}$')
            ax.scatter([T_scan[0] * 0.3], [avg_vals[j]],
                       color=colors[j], marker='*', s=150, zorder=5)

        ax.scatter([], [], color='grey', marker='*', s=150,
                   label=r'$T\to 0$ (averaged ODE)')
        ax.set_xscale('log')
        ax.set_xlabel('Switching period $T$')
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle='--', alpha=0.25)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, loc='best')
        fig.tight_layout()

        fname = f'{fname_base}_{suffix}'
        fig.savefig(f'{outdir}/{fname}.pdf', dpi=250,
                    facecolor='none', bbox_inches='tight')
        fig.savefig(f'{outdir}/{fname}.svg', dpi=250,
                    facecolor='none', bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {fname}.pdf / .svg")


if __name__ == '__main__':
    main(
        rA=1.0, KA=1.0, gA=1.01,
        rB=4.0, KB=2.0, gB=2.2,
        p_values=np.array([0.3, 0.5, 0.667, 0.8, 0.9]),
        T_scan=np.logspace(-4, 1, 100),
        steps_per_period=300,
        transient_periods=500,
        measure_periods=50,
        outdir='src/3paramLogisticLV/plots/ystar_vs_T',
    )