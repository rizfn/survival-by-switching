#!/usr/bin/env python3
"""Plot three quantities vs duty cycle p.

Switching is characterised by its rate alpha (as in the main text): the full
period is 2/alpha, of which env A is held a fraction p (a dwell 2p/alpha) and
env B a fraction (1-p). At p=1/2 this is the equal dwell 1/alpha per environment.

Panel 1 — y_max(alpha, p):  peak predator density on the 2D limit cycle
                            (at the A→B switch moment), computed via batched
                            vectorised RK4 over all (alpha, p) simultaneously.
                            The alpha→∞ limit is y*(p) analytically.
Panel 2 — λ(alpha, p):      invasion exponent family + analytical alpha→∞ limit.
Panel 3 — alpha_c(p):       critical switching rate (the predator persists for
                            alpha > alpha_c; the most robust p minimises alpha_c).

Parameters:
  Env A: r=1, K=1, γ=1.1
  Env B: r=4, K=2, γ=2.1
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.optimize import brentq
import os


def label_ends_only(ax):
    """Keep every default tick mark, but only write the tick label for the
    smallest and largest tick that falls within the current view limits.
    (e.g. ticks [0,1,2,3,4,5] with the view ending at 5.5 -> only 0 and 5 are
    labelled.) Done with a draw-time formatter so it tracks the auto ticks,
    rather than rounding the data range by hand."""
    def make_fmt(get_ticks, get_lim):
        def fmt(v, pos):
            lo, hi = sorted(get_lim())
            inview = [t for t in get_ticks() if lo - 1e-9 <= t <= hi + 1e-9]
            if not inview:
                return ''
            if abs(v - min(inview)) < 1e-9 or abs(v - max(inview)) < 1e-9:
                return f'{v:g}'
            return ''
        return fmt
    ax.xaxis.set_major_formatter(FuncFormatter(make_fmt(ax.get_xticks, ax.get_xlim)))
    ax.yaxis.set_major_formatter(FuncFormatter(make_fmt(ax.get_yticks, ax.get_ylim)))


def main(rA, KA, gA, rB, KB, gB, ALPHA_FAMILY, N_P):

    # ── Parameters ───────────────────────────────────────────────────────────
    cA, cB = rA / KA, rB / KB
    p_arr  = np.linspace(0.01, 1.0, N_P)

    # ── Averaged parameters (fast-switching) ─────────────────────────────────
    def r_bar(p): return p * rA + (1 - p) * rB
    def g_bar(p): return p * gA + (1 - p) * gB
    def c_bar(p): return p * cA + (1 - p) * cB

    # ── Panel 1: y_max via batched vectorised RK4 ─────────────────────────────

    def y_star(p):
        """Fast-switching (alpha→∞) predator steady state.
        The analytical expression can go negative, which is unphysical
        (the predator is extinct), so we cap it at 0."""
        return max(r_bar(p) - c_bar(p) * g_bar(p), 0.0)

    def p_opt_ystar():
        """Analytical p that maximises y*(p)."""
        Dr, Dg, Dc = rA - rB, gA - gB, cA - cB
        denom = 2 * Dc * Dg
        return np.nan if abs(denom) < 1e-14 else (Dr - cB * Dg - gB * Dc) / denom

    def compute_ymax_all(T_family, p_arr,
                         steps_per_period: int = 200,
                         transient_periods: int = 500) -> np.ndarray:
        """Batched RK4 over the full (N_T × N_p) grid simultaneously."""
        T_arr = np.asarray(T_family, dtype=float)
        N_T, N_p = len(T_arr), len(p_arr)

        T_grid  = T_arr[:, None]
        p_grid  = p_arr[None, :]
        dt_grid = (T_arr / steps_per_period)[:, None]
        tau_A   = p_grid * T_grid

        def rhs(x, logy, t_arr):
            phase = t_arr % T_grid
            in_A  = phase < tau_A
            r = np.where(in_A, rA, rB)
            K = np.where(in_A, KA, KB)
            g = np.where(in_A, gA, gB)
            y = np.exp(np.clip(logy, -500, 500))
            return r * x * (1 - x / K) - x * y, x - g

        def rk4(x, logy, t_arr):
            k1x, k1l = rhs(x,                  logy,                  t_arr)
            k2x, k2l = rhs(x + .5*dt_grid*k1x, logy + .5*dt_grid*k1l, t_arr + .5*dt_grid)
            k3x, k3l = rhs(x + .5*dt_grid*k2x, logy + .5*dt_grid*k2l, t_arr + .5*dt_grid)
            k4x, k4l = rhs(x +    dt_grid*k3x, logy +    dt_grid*k3l, t_arr +    dt_grid)
            x_new    = np.maximum(x    + (dt_grid/6)*(k1x+2*k2x+2*k3x+k4x), 1e-14)
            logy_new = np.maximum(logy + (dt_grid/6)*(k1l+2*k2l+2*k3l+k4l), np.log(1e-14))
            return x_new, logy_new

        x     = np.full((N_T, N_p), 0.8 * KA)
        logy  = np.full((N_T, N_p), np.log(0.01))
        t_arr = np.zeros((N_T, 1))

        for _ in range(transient_periods * steps_per_period):
            x, logy = rk4(x, logy, t_arr)
            t_arr  += dt_grid

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

    # ── Panel 2: invasion exponent ────────────────────────────────────────────

    def lambda_fast(p):
        c = c_bar(p)
        return np.nan if abs(c) < 1e-14 else r_bar(p) / c - g_bar(p)

    def p_opt_lambda():
        Dr, Dg, Dc = rA - rB, gA - gB, cA - cB
        if abs(Dc) < 1e-14 or abs(Dg) < 1e-14: return np.nan
        inner = (cB * Dr - rB * Dc) / Dg
        return np.nan if inner < 0 else (np.sqrt(inner) - cB) / Dc

    def phi(t, x_in, K, r):
        if t == 0: return x_in
        rt = r * t
        if rt > 500: return K * (1.0 - x_in / K * np.exp(-rt))
        expo = np.exp(rt)
        return K * x_in * expo / (K + x_in * (expo - 1.0))

    def periodic_fixed_points(tau_A, tau_B):
        def residual(xB):
            return phi(tau_B, phi(tau_A, xB, KA, rA), KB, rB) - xB
        lo, hi = min(KA, KB) * 1e-4, max(KA, KB) * (1 - 1e-6)
        try:
            flo, fhi = residual(lo), residual(hi)
            if flo * fhi > 0:
                xs = np.linspace(lo, hi, 512)
                for i in range(len(xs) - 1):
                    if residual(xs[i]) * residual(xs[i + 1]) < 0:
                        lo, hi = xs[i], xs[i + 1]; break
                else: return np.nan, np.nan
            xB_star = brentq(residual, lo, hi, xtol=1e-10)
            return phi(tau_A, xB_star, KA, rA), xB_star
        except Exception:
            return np.nan, np.nan

    def invasion_exponent(alpha, p):
        # Switching rate alpha: env A held 2p/alpha, env B held 2(1-p)/alpha, so
        # the full period is 2/alpha (equal dwell 1/alpha at p=1/2). The old 1/T
        # prefactor of the period form becomes alpha/2.
        tau_A, tau_B = 2.0 * p / alpha, 2.0 * (1.0 - p) / alpha
        xA, xB = periodic_fixed_points(tau_A, tau_B)
        if not (np.isfinite(xA) and np.isfinite(xB) and xA > 0 and xB > 0):
            return np.nan
        drift    = (KA - gA) * p + (KB - gB) * (1 - p)
        log_term = (KA / rA) * np.log(xA / xB) + (KB / rB) * np.log(xB / xA)
        return drift - (alpha / 2.0) * log_term

    # ── Panel 3: alpha_c(p) ──────────────────────────────────────────────────

    def alpha_c(p, a_lo=4e-3, a_hi=2e2):
        # Critical switching rate: the root of lambda(alpha, p) = 0. Below it the
        # predator dies (slow switching), above it it persists. The bracket
        # [4e-3, 2e2] mirrors the old period bracket [1e-2, 500] via alpha = 2/T.
        try:
            a_grid = np.logspace(np.log10(a_lo), np.log10(a_hi), 200)
            lam    = np.array([invasion_exponent(a, p) for a in a_grid])
            for i in range(len(a_grid) - 1):
                f0, f1 = lam[i], lam[i + 1]
                if np.isfinite(f0) and np.isfinite(f1) and f0 * f1 < 0:
                    return brentq(lambda a: invasion_exponent(a, p),
                                  a_grid[i], a_grid[i + 1], xtol=1e-6, maxiter=200)
            return np.nan
        except Exception:
            return np.nan

    # ── Compute all curves ────────────────────────────────────────────────────

    # Switching-rate family -> internal periods (T = 2/alpha) for the RK4 sim.
    periods = 2.0 / np.asarray(ALPHA_FAMILY, dtype=float)

    print("Computing y_max(alpha, p) via batched RK4...")
    ymax_grid = compute_ymax_all(periods, p_arr, steps_per_period=200, transient_periods=500)
    print(f"  done — shape {ymax_grid.shape}")

    print("Computing y*(p) [alpha→∞ analytical limit]...")
    ymax_fast = np.array([y_star(p) for p in p_arr])

    print("Computing λ_fast(p)...")
    lam_fast = np.array([lambda_fast(p) for p in p_arr])

    print("Computing λ(alpha, p) family...")
    lam_family = {}
    for a_val in ALPHA_FAMILY:
        lam_family[a_val] = np.array([invasion_exponent(a_val, p) for p in p_arr])
        print(f"  alpha={a_val} done")

    print("Computing alpha_c(p)...")
    ac_arr = np.array([alpha_c(p) for p in p_arr])

    # Optima
    p_opt1 = p_opt_ystar()
    p_opt2 = p_opt_lambda()
    # Most robust duty cycle = smallest critical rate (survives down to the
    # slowest switching), i.e. argmin of alpha_c -- the mirror of the old max-T_c.
    valid_ac   = np.isfinite(ac_arr)
    p_opt3_idx = np.nanargmin(ac_arr[valid_ac].reshape(-1)) if np.any(valid_ac) else None
    p_opt3     = p_arr[valid_ac][p_opt3_idx] if p_opt3_idx is not None else np.nan
    ac_opt3    = ac_arr[valid_ac][p_opt3_idx] if p_opt3_idx is not None else np.nan

    print(f"\np_opt (y*)         = {p_opt1:.4f}")
    print(f"p_opt (λ_fast)     = {p_opt2:.4f}")
    print(f"p_opt (α_c min)    = {p_opt3:.4f}")

    # ── Colormap shared across all panels ─────────────────────────────────────
    # alpha-family curves: cool_r (blue=fast/large-alpha, pink=slow/small-alpha)
    # alpha→∞  limit: bold accent colour  (used in panels 1 & 2)
    # alpha_c(p): accent colour, panel 3

    log_a    = np.log10(np.asarray(ALPHA_FAMILY, dtype=float))
    # Keep blue = fast: the largest alpha maps to 0 (the 'fast' end of cool_r),
    # the smallest alpha to 1 (pink, slow).
    a_norm   = (log_a.max() - log_a) / (log_a.max() - log_a.min())
    cmap     = plt.get_cmap('cool_r')
    A_colors = [cmap(n) for n in a_norm]

    C_FAST = "#910a2c"   # accent for the alpha→∞ analytical limit (panels 1 & 2)
    C_OPT  = "#2a2a2a"   # near-black for optimal p markers (all panels)

    # ── Plot ──────────────────────────────────────────────────────────────────

    plt.rcParams.update({
        'font.size': 24, 'axes.labelsize': 24,
        'xtick.labelsize': 20, 'ytick.labelsize': 20, 'legend.fontsize': 18,
    })

    fig, axes = plt.subplots(1, 3, figsize=(20, 6), layout='constrained')
    fig.patch.set_facecolor('none')
    for ax in axes:
        ax.set_facecolor('none')

    def mark_opt(ax, p_opt, y_opt, label):
        """Shared helper: vertical dotted line + scatter dot + legend entry."""
        ax.axvline(p_opt, color=C_OPT, linestyle=':', linewidth=2.0, alpha=0.9, zorder=7)
        ax.scatter([p_opt], [y_opt], color=C_OPT, zorder=8, s=80, label=label)

    # ── Panel 1: y_max ────────────────────────────────────────────────────────
    ax = axes[0]
    ax.axhline(0.0, color='black', linewidth=1.0, linestyle='--', alpha=0.45)

    # alpha→∞ analytical limit
    ax.plot(p_arr, ymax_fast, color=C_FAST, linewidth=2.8, zorder=5,
            label=r'$\alpha\to\infty$ (analytical)')

    for i, a_val in enumerate(ALPHA_FAMILY):
        ym = np.where(ymax_grid[i] == 0.0, np.nan, ymax_grid[i])
        ax.plot(p_arr, ym, color=A_colors[i], linewidth=1.6, alpha=0.85,
                label=fr'$\alpha={a_val}$')

    # p_opt marker on the alpha→∞ curve
    if np.isfinite(p_opt1):
        y_at1 = np.interp(p_opt1, p_arr, ymax_fast)
        mark_opt(ax, p_opt1, y_at1, fr'$p_{{\rm opt}} \approx {p_opt1:.3f}$')

    ax.set_xlabel('Duty cycle $p$', labelpad=-18)
    ax.set_ylabel(r'Peak predator $y_{\max}$', labelpad=-28)
    ax.set_xlim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='lower right', ncol=2, columnspacing=0.8, handlelength=1.4)

    # ── Panel 2: λ ────────────────────────────────────────────────────────────
    ax = axes[1]
    ax.axhline(0.0, color='black', linewidth=1.0, linestyle='--', alpha=0.45)

    ax.plot(p_arr, lam_fast, color=C_FAST, linewidth=2.8, zorder=5,
            label=r'$\alpha\to\infty$ (analytical)')

    for i, a_val in enumerate(ALPHA_FAMILY):
        ax.plot(p_arr, lam_family[a_val], color=A_colors[i], linewidth=1.6, alpha=0.80,
                label=fr'$\alpha={a_val}$')

    if np.isfinite(p_opt2):
        lam_at = np.interp(p_opt2, p_arr, lam_fast)
        mark_opt(ax, p_opt2, lam_at, fr'$p_{{\rm opt}}^\lambda \approx {p_opt2:.3f}$')

    ax.set_xlabel('Duty cycle $p$', labelpad=-18)
    ax.set_ylabel(r'Invasion exponent $\lambda$', labelpad=-36)
    ax.set_xlim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='lower right', ncol=2, columnspacing=0.8, handlelength=1.4)

    # ── Panel 3: alpha_c ─────────────────────────────────────────────────────
    ax = axes[2]
    ax.plot(p_arr, ac_arr, color=C_FAST, linewidth=2.5, label=r'$\alpha_c(p)$')

    if np.isfinite(p_opt3):
        mark_opt(ax, p_opt3, ac_opt3, fr'$p_{{\rm opt}}^{{\alpha_c}} \approx {p_opt3:.3f}$')

    ax.set_xlabel('Duty cycle $p$', labelpad=-18)
    ax.set_ylabel(r'Critical rate $\alpha_c$', labelpad=-12)
    ax.set_xlim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='best')

    # ── Keep all ticks, but label only the extreme ticks on every axis ──────────
    for ax in axes:
        label_ends_only(ax)

    # ── Save ──────────────────────────────────────────────────────────────────
    outdir = 'src/3paramLogisticLV/plots/duty_cycle'
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(f'{outdir}/r_{rA}-{rB}_K_{KA}-{KB}_g_{gA}-{gB}.pdf', dpi=250, facecolor='none', bbox_inches='tight')
    fig.savefig(f'{outdir}/r_{rA}-{rB}_K_{KA}-{KB}_g_{gA}-{gB}.svg', dpi=250, facecolor='none', bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved r_{rA}-{rB}_K_{KA}-{KB}_g_{gA}-{gB}.pdf / .svg")


if __name__ == "__main__":
    rA, KA, gA   = 1.0, 1.0, 1.1
    rB, KB, gB   = 4.0, 2.0, 2.1
    # Switching-rate family (fast → slow); internal periods are 2/alpha = [0.5,1,2,4,10,20]
    ALPHA_FAMILY = [4.0, 2.0, 1.0, 0.5, 0.2, 0.1]
    N_P          = 500
    main(rA, KA, gA, rB, KB, gB, ALPHA_FAMILY, N_P)