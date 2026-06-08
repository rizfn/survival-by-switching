#!/usr/bin/env python3
"""Plot three quantities vs duty cycle p. Everything is numerical simulation.

The invasion exponent is computed via a two-phase approach:
  Phase 1 — settle x on the prey-only orbit (y=0 exactly) for each (T,p).
  Phase 2 — introduce y=eps and measure d(ln y)/dt over several periods.
This avoids the spurious lambda=0 that occurs when x starts near the
coexistence fixed point instead of the prey-only orbit.

Panel 1 — y_max(T, p):  peak predator on the 2D limit cycle (batched RK4).
                         T→0 limit: steady-state y of the averaged ODE.
Panel 2 — λ(T, p):      invasion exponent (two-phase RK4).
                         T→0 limit: invasion exponent of the averaged ODE.
Panel 3 — T_c(p):       zero-crossing of λ(T, p).
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import os


def main(rA, KA, gA, rB, KB, gB, T_FAMILY, N_P):

    cA, cB = rA / KA, rB / KB
    p_arr  = np.linspace(0.01, 1.0, N_P)

    def r_bar(p): return p * rA + (1 - p) * rB
    def g_bar(p): return p * gA + (1 - p) * gB
    def c_bar(p): return p * cA + (1 - p) * cB

    # ── Switching RK4 kernel ──────────────────────────────────────────────────

    def _make_switch(T_arr, p_arr, steps_per_period):
        T_arr   = np.asarray(T_arr,  dtype=float)
        p_arr   = np.asarray(p_arr,  dtype=float)
        T_grid  = T_arr[:, None]
        p_grid  = p_arr[None, :]
        dt_grid = (T_arr / steps_per_period)[:, None]
        tau_A   = p_grid * T_grid

        def rhs(x, logy, t_arr, y_zero=False):
            phase = t_arr % T_grid
            in_A  = phase < tau_A
            r = np.where(in_A, rA, rB)
            K = np.where(in_A, KA, KB)
            g = np.where(in_A, gA, gB)
            if y_zero:
                y = np.zeros_like(x)
            else:
                y = np.exp(np.clip(logy, -500, 500))
            return r * x * (1 - x / K) - x * y, x - g

        def rk4(x, logy, t_arr, y_zero=False):
            k1x, k1l = rhs(x,                  logy,                  t_arr, y_zero)
            k2x, k2l = rhs(x + .5*dt_grid*k1x, logy + .5*dt_grid*k1l, t_arr + .5*dt_grid, y_zero)
            k3x, k3l = rhs(x + .5*dt_grid*k2x, logy + .5*dt_grid*k2l, t_arr + .5*dt_grid, y_zero)
            k4x, k4l = rhs(x +    dt_grid*k3x, logy +    dt_grid*k3l, t_arr +    dt_grid,  y_zero)
            x_new    = np.maximum(x    + (dt_grid/6)*(k1x+2*k2x+2*k3x+k4x), 1e-14)
            logy_new =            logy + (dt_grid/6)*(k1l+2*k2l+2*k3l+k4l)
            return x_new, logy_new

        return rk4, dt_grid, tau_A, T_grid

    # ── Averaged ODE kernel ───────────────────────────────────────────────────

    def _make_avg(p_arr, dt):
        p_arr = np.asarray(p_arr, dtype=float)
        r_b, g_b, c_b = r_bar(p_arr), g_bar(p_arr), c_bar(p_arr)

        def rhs_avg(x, logy, y_zero=False):
            y = np.zeros_like(x) if y_zero else np.exp(np.clip(logy, -500, 500))
            return r_b*x - c_b*x**2 - x*y, x - g_b

        def rk4_avg(x, logy, y_zero=False):
            k1x, k1l = rhs_avg(x,             logy,            y_zero)
            k2x, k2l = rhs_avg(x+.5*dt*k1x,   logy+.5*dt*k1l, y_zero)
            k3x, k3l = rhs_avg(x+.5*dt*k2x,   logy+.5*dt*k2l, y_zero)
            k4x, k4l = rhs_avg(x+   dt*k3x,   logy+   dt*k3l, y_zero)
            x_new    = np.maximum(x    + (dt/6)*(k1x+2*k2x+2*k3x+k4x), 1e-14)
            logy_new =            logy + (dt/6)*(k1l+2*k2l+2*k3l+k4l)
            return x_new, logy_new

        return rk4_avg

    # ── T→0: averaged ODE steady-state y ─────────────────────────────────────

    def compute_avg_ymax(p_arr, dt=0.01, transient_time=500.0):
        rk4_avg = _make_avg(p_arr, dt)
        x    = np.full(len(p_arr), 0.8 * KA)
        logy = np.full(len(p_arr), np.log(0.01))
        for _ in range(int(transient_time / dt)):
            x, logy = rk4_avg(x, logy)
        return np.exp(np.clip(logy, -500, 500))

    def compute_avg_lambda(p_arr, eps=1e-10, dt=0.01,
                           settle_time=500.0, measure_time=100.0):
        """
        Two-phase invasion exponent for averaged ODE.
        Phase 1: settle x on prey-only equilibrium with y=0.
        Phase 2: introduce y=eps, measure d(ln y)/dt.
        """
        rk4_avg = _make_avg(p_arr, dt)
        # Phase 1: y=0
        x    = np.full(len(p_arr), 0.5*(KA+KB))
        logy = np.full(len(p_arr), np.log(eps))
        for _ in range(int(settle_time / dt)):
            x, logy = rk4_avg(x, logy, y_zero=True)
        # Phase 2: y=eps, measure
        logy = np.full(len(p_arr), np.log(eps))
        logy0 = logy.copy()
        for _ in range(int(measure_time / dt)):
            x, logy = rk4_avg(x, logy)
        return (logy - logy0) / measure_time

    # ── Panel 1: y_max (full 2D limit cycle) ─────────────────────────────────

    def compute_ymax_all(T_family, p_arr,
                         steps_per_period=200, transient_periods=500):
        T_arr = np.asarray(T_family, dtype=float)
        N_T, N_p = len(T_arr), len(p_arr)
        rk4, dt_grid, tau_A, T_grid = _make_switch(T_arr, p_arr, steps_per_period)

        x     = np.full((N_T, N_p), 0.8*KA)
        logy  = np.full((N_T, N_p), np.log(0.01))
        t_arr = np.zeros((N_T, 1))

        for _ in range(transient_periods * steps_per_period):
            x, logy = rk4(x, logy, t_arr)
            t_arr  += dt_grid

        logy = np.maximum(logy, np.log(1e-14))

        y_max = np.full((N_T, N_p), np.nan)
        for _ in range(steps_per_period):
            phase      = t_arr % T_grid
            next_phase = (t_arr + dt_grid) % T_grid
            crossing   = (phase < tau_A) & ((next_phase >= tau_A) | (next_phase < phase))
            y_now      = np.exp(np.clip(logy, -500, 500))
            y_max      = np.where(crossing, y_now, y_max)
            x, logy    = rk4(x, logy, t_arr)
            t_arr     += dt_grid

        extinct = np.exp(np.clip(logy, -500, 500)) < 1e-10
        y_max[extinct] = 0.0
        return y_max

    # ── Panel 2 & 3: invasion exponent (two-phase) ───────────────────────────

    def compute_lambda_2phase(T_arr, p_arr, eps=1e-10,
                              steps_per_period=400,
                              settle_periods=300,
                              measure_periods=50):
        """
        Two-phase invasion exponent for switching ODE.
        Phase 1: integrate with y=0 to settle x on the prey-only periodic orbit.
        Phase 2: set y=eps and measure d(ln y)/dt over measure_periods.
        """
        T_arr = np.asarray(T_arr, dtype=float)
        N_T, N_p = len(T_arr), len(p_arr)
        rk4, dt_grid, tau_A, T_grid = _make_switch(T_arr, p_arr, steps_per_period)

        # Phase 1: y=0, settle x
        x     = np.full((N_T, N_p), 0.5*(KA+KB))
        logy  = np.full((N_T, N_p), np.log(eps))
        t_arr = np.zeros((N_T, 1))
        for _ in range(settle_periods * steps_per_period):
            x, logy = rk4(x, logy, t_arr, y_zero=True)
            t_arr  += dt_grid

        # Phase 2: y=eps, measure
        logy  = np.full((N_T, N_p), np.log(eps))
        logy0 = logy.copy()
        for _ in range(measure_periods * steps_per_period):
            x, logy = rk4(x, logy, t_arr)
            t_arr  += dt_grid

        return (logy - logy0) / (measure_periods * T_grid)

    def compute_tc(lam_grid, T_scan, p_arr):
        N_T, N_p = lam_grid.shape
        tc = np.full(N_p, np.nan)
        for j in range(N_p):
            lams = lam_grid[:, j]
            for i in range(N_T - 1):
                f0, f1 = lams[i], lams[i+1]
                if np.isfinite(f0) and np.isfinite(f1) and f0*f1 < 0:
                    tc[j] = T_scan[i] + (T_scan[i+1]-T_scan[i])*(-f0)/(f1-f0)
                    break
        return tc

    # ── Compute ───────────────────────────────────────────────────────────────

    print("Computing y_max(T→0) [averaged ODE]...")
    ymax_fast = compute_avg_ymax(p_arr)

    print("Computing λ(T→0) [averaged ODE, two-phase]...")
    lam_fast = compute_avg_lambda(p_arr)

    print("Computing y_max(T, p) [2D limit cycle]...")
    ymax_grid = compute_ymax_all(T_FAMILY, p_arr,
                                 steps_per_period=200, transient_periods=500)
    print(f"  done — shape {ymax_grid.shape}")

    print("Computing λ(T, p) family [two-phase, per T]...")
    lam_family = {}
    for T_val in T_FAMILY:
        row = compute_lambda_2phase(
            [T_val], p_arr,
            eps=1e-10, steps_per_period=400,
            settle_periods=300, measure_periods=50,
        )
        lam_family[T_val] = row[0]
        print(f"  T={T_val} done")

    print("Computing T_c(p) [two-phase, per-decade T scan]...")
    lam_parts, T_parts = [], []
    for lo, hi, n in [(-1, 0, 25), (0, 1, 25), (1, 2, 25), (2, 3, 15)]:
        T_part = np.logspace(lo, hi, n)
        lam_part = compute_lambda_2phase(
            T_part, p_arr,
            eps=1e-10, steps_per_period=400,
            settle_periods=200, measure_periods=30,
        )
        T_parts.append(T_part)
        lam_parts.append(lam_part)
        print(f"  decade [{10**lo:.1f}, {10**hi:.0f}] done")

    T_scan   = np.concatenate(T_parts)
    lam_scan = np.concatenate(lam_parts, axis=0)
    _, uid   = np.unique(T_scan, return_index=True)
    T_scan, lam_scan = T_scan[uid], lam_scan[uid]
    tc_arr = compute_tc(lam_scan, T_scan, p_arr)

    # Optima
    p_opt1     = p_arr[np.nanargmax(np.where(np.isfinite(ymax_fast), ymax_fast, np.nan))]
    p_opt2     = p_arr[np.nanargmax(np.where(np.isfinite(lam_fast),  lam_fast,  np.nan))]
    valid_tc   = np.isfinite(tc_arr)
    p_opt3_idx = np.nanargmax(tc_arr[valid_tc]) if np.any(valid_tc) else None
    p_opt3     = p_arr[valid_tc][p_opt3_idx]    if p_opt3_idx is not None else np.nan
    tc_opt3    = tc_arr[valid_tc][p_opt3_idx]   if p_opt3_idx is not None else np.nan

    print(f"\np_opt (y*, T→0)  = {p_opt1:.4f}")
    print(f"p_opt (λ, T→0)   = {p_opt2:.4f}")
    print(f"p_opt (T_c max)  = {p_opt3:.4f}")

    # ── Colormap ──────────────────────────────────────────────────────────────
    log_T    = np.log10(T_FAMILY)
    t_norm   = (log_T - log_T.min()) / (log_T.max() - log_T.min())
    cmap     = plt.get_cmap('cool_r')
    T_colors = [cmap(n) for n in t_norm]
    C_FAST, C_OPT = "#910a2c", "#2a2a2a"

    # ── Plot ──────────────────────────────────────────────────────────────────
    plt.rcParams.update({
        'font.size': 22, 'axes.labelsize': 22,
        'xtick.labelsize': 18, 'ytick.labelsize': 18, 'legend.fontsize': 16,
    })
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), layout='constrained')
    fig.patch.set_facecolor('none')
    for ax in axes: ax.set_facecolor('none')

    def mark_opt(ax, p_opt, y_opt, label):
        ax.axvline(p_opt, color=C_OPT, linestyle=':', linewidth=2.0, alpha=0.9, zorder=7)
        ax.scatter([p_opt], [y_opt], color=C_OPT, zorder=8, s=80, label=label)

    # Panel 1
    ax = axes[0]
    ax.axhline(0.0, color='black', linewidth=1.0, linestyle='--', alpha=0.45)
    ax.plot(p_arr, ymax_fast, color=C_FAST, linewidth=2.8, zorder=5,
            label=r'$T\to 0$ (simulation)')
    for i, T_val in enumerate(T_FAMILY):
        ym = np.where(ymax_grid[i] == 0.0, np.nan, ymax_grid[i])
        ax.plot(p_arr, ym, color=T_colors[i], linewidth=1.6, alpha=0.85,
                label=f'$T={T_val}$')
    mark_opt(ax, p_opt1, np.interp(p_opt1, p_arr, ymax_fast),
             fr'$p_{{\rm opt}} \approx {p_opt1:.3f}$')
    ax.set_xlabel('Duty cycle $p$')
    ax.set_ylabel(r'Peak predator $y_{\max}$ on limit cycle')
    ax.set_xlim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.25); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='lower right', ncol=2,
              columnspacing=0.8, handlelength=1.4)

    # Panel 2
    ax = axes[1]
    ax.axhline(0.0, color='black', linewidth=1.0, linestyle='--', alpha=0.45)
    ax.plot(p_arr, lam_fast, color=C_FAST, linewidth=2.8, zorder=5,
            label=r'$T\to 0$ (simulation)')
    for i, T_val in enumerate(T_FAMILY):
        ax.plot(p_arr, lam_family[T_val], color=T_colors[i], linewidth=1.6,
                alpha=0.80, label=f'$T={T_val}$')
    mark_opt(ax, p_opt2, np.interp(p_opt2, p_arr, lam_fast),
             fr'$p_{{\rm opt}}^\lambda \approx {p_opt2:.3f}$')
    ax.set_xlabel('Duty cycle $p$')
    ax.set_ylabel(r'Invasion exponent $\lambda$')
    ax.set_xlim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.25); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='lower right', ncol=2,
              columnspacing=0.8, handlelength=1.4)

    # Panel 3
    ax = axes[2]
    ax.plot(p_arr, tc_arr, color=C_FAST, linewidth=2.5, label=r'$T_c(p)$')
    if np.isfinite(p_opt3):
        mark_opt(ax, p_opt3, tc_opt3,
                 fr'$p_{{\rm opt}}^{{T_c}} \approx {p_opt3:.3f}$')
    ax.set_xlabel('Duty cycle $p$')
    ax.set_ylabel(r'Critical period $T_c$')
    ax.set_xlim(0, 1); ax.set_yscale('log')
    ax.grid(True, linestyle='--', alpha=0.25); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='best')

    # Save
    outdir = 'src/3paramLogisticLV/plots/duty_cycle_sim'
    os.makedirs(outdir, exist_ok=True)
    fname = f'r_{rA}-{rB}_K_{KA}-{KB}_g_{gA}-{gB}'
    fig.savefig(f'{outdir}/{fname}.pdf', dpi=250, facecolor='none', bbox_inches='tight')
    fig.savefig(f'{outdir}/{fname}.svg', dpi=250, facecolor='none', bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved {outdir}/{fname}.pdf / .svg")


if __name__ == "__main__":
    rA, KA, gA = 1.0, 1.0, 1.1
    rB, KB, gB = 4.0, 2.0, 2.2
    T_FAMILY   = [0.5, 1.0, 2.0, 4.0, 10.0, 20.0]
    N_P        = 250
    main(rA, KA, gA, rB, KB, gB, T_FAMILY, N_P)