"""
Two-panel figure showing how multiple quantities all vanish at the same critical threshold.

Left:  Base LV model — quantities vs gamma_2  (threshold gamma_2* = 3)
Right: Shielding model — quantities vs alpha_2 (threshold alpha_2* = ln 2 ≈ 0.693)

Fixed params:  r1=K1=gamma1=1, r2=4, K2=2.

Quantities plotted (raw, un-normalised, with twin axes):
  1. Integral of R(x)+ (positive part) over (K1, K2)   — amber / gold (left y-axis)
  2. Optimal y* (max over p in (0,1))                   — blue  (right y-axis)
  3. Max invasion exponent (max_p lambda_inv)            — crimson (right y-axis)
"""

from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os

# ─────────────────────────────────────────────────────────────────────────────
r1, K1, g1 = 1.0, 1.0, 1.0
r2, K2     = 4.0, 2.0

AMBER_FILL = '#CBA810'
BLUE       = '#2d4f73'
CRIMSON    = '#901A1E'

# ─────────────────────────────────────────────────────────────────────────────
# BASE LV
# ─────────────────────────────────────────────────────────────────────────────

def R_base(x, g2):
    h2_f2 = (x - g2) / (2 * x * (2 - x))
    h1_f1 = -1.0 / x
    return h2_f2 - h1_f1

def base_R_zero(g2):
    xz = 4.0 - g2
    if K1 < xz < K2:
        return xz
    return None

def base_quantities(g2_vals):
    areas, ystars, lambdas = [], [], []

    for g2 in g2_vals:
        # ── Integral of R(x)+ over full (K1, K2) ─────────────────────────
        xs = np.linspace(K1 + 1e-6, K2 - 1e-6, 4000)
        r_vals = np.array([R_base(x, g2) for x in xs])
        area = np.trapezoid(np.maximum(r_vals, 0.0), xs)
        areas.append(area)

        # ── y* maximised over p ───────────────────────────────────────────
        # y* = (1-p)[4 - p - g2(2-p)], optimum at p_opt = (5-3g2)/(2-2g2)
        if abs(g2 - 1.0) < 1e-9:
            p_opt = 0.5
        else:
            p_opt = (5 - 3*g2) / (2 - 2*g2)

        def ystar(p):
            return max((1 - p) * (4 - p - g2 * (2 - p)), 0.0)

        p_clamped = max(1e-4, min(0.9999, p_opt))
        ystars.append(max(ystar(p_clamped), ystar(1e-4), ystar(0.9999)))

        # ── Invasion exponent maximised over p ────────────────────────────
        def x0_base(p):
            num = p*r1 + (1-p)*r2
            den = p*r1/K1 + (1-p)*r2/K2
            return num / den

        def lambda_inv_base(p):
            x0 = x0_base(p)
            return p*(x0 - g1) + (1-p)*(x0 - g2)

        lam_vals = [lambda_inv_base(p) for p in np.linspace(0.01, 0.99, 500)]
        lambdas.append(max(lam_vals))

    return np.array(areas), np.array(ystars), np.array(lambdas)


# ─────────────────────────────────────────────────────────────────────────────
# SHIELDING MODEL
# ─────────────────────────────────────────────────────────────────────────────

def R_shield(x, aB, g=1.0, KB=2.0):
    KA = 1.0
    h_A_f_A = -1.0 / x
    h_B = x * np.exp(-aB * x) - g
    f_B = x * (1 - x / KB)
    return h_B / f_B - h_A_f_A

def shield_quantities(aB_vals):
    KA, KB, g = 1.0, 2.0, 1.0
    areas, ystars, lambdas = [], [], []

    for aB in aB_vals:
        # ── Integral of R(x)+ over full (KA, KB) ─────────────────────────
        xs = np.linspace(KA + 1e-6, KB - 1e-6, 4000)
        r_vals = np.array([R_shield(x, aB) for x in xs])
        area = np.trapezoid(np.maximum(r_vals, 0.0), xs)
        areas.append(area)

        # ── y* maximised over p (parametrise by x* ∈ (KA, KB)) ──────────
        def quantities_at_xstar(xstar):
            s = np.exp(-aB * xstar)
            denom_p = xstar * (1 - s)
            if denom_p < 1e-12:
                return 0.0
            p = (g - xstar * s) / denom_p
            if not (0 < p < 1):
                return 0.0
            num = (g - xstar*s)*(1 - xstar/KA) + (xstar - g)*(1 - xstar/KB)
            den = g * (1 - s)
            return max(num / den, 0.0) if den > 0 else 0.0

        xstar_grid = np.linspace(KA + 1e-5, KB - 1e-5, 1000)
        ystars.append(max(quantities_at_xstar(xs) for xs in xstar_grid))

        # ── Invasion exponent maximised over p ────────────────────────────
        def lambda_inv_shield(p):
            x0 = 1.0 / (p/KA + (1-p)/KB)
            return p*(x0 - g) + (1-p)*(x0*np.exp(-aB*x0) - g)

        lam_vals = [lambda_inv_shield(p) for p in np.linspace(0.01, 0.99, 500)]
        lambdas.append(max(lam_vals))

    return np.array(areas), np.array(ystars), np.array(lambdas)


# ─────────────────────────────────────────────────────────────────────────────
# PANEL DRAWING
# ─────────────────────────────────────────────────────────────────────────────

def draw_panel(ax, xvals, areas, ystars, lambdas, x_crit, xlabel, 
               inset_xlim, show_ylabel_left=True, show_ylabel_right=True):
    """
    Main axes:  left y-axis  → amber R(x)+ area
                right y-axis → blue y*, crimson Lambda_inv
    Inset: zoom near threshold showing all three vanishing together.
    """

    # ── Twin axes ─────────────────────────────────────────────────────────
    ax_r = ax.twinx()

    # Fill + line for R area on left axis
    ax.fill_between(xvals, areas, color='#f5e6a3', alpha=0.75, zorder=0)
    l1, = ax.plot(xvals, areas, color=AMBER_FILL, lw=2.6, zorder=4,
                  label=r'$\int_{K_1}^{K_2} [R(x)]_+\,dx$')

    # y* and Lambda on right axis
    l2, = ax_r.plot(xvals, ystars,  color=BLUE,    lw=2.4, zorder=5,
                    label=r'$\max_p\, y^*$')
    l3, = ax_r.plot(xvals, lambdas, color=CRIMSON, lw=2.4, zorder=5, ls='-',
                    label=r'$\max_p\, \Lambda_{\rm inv}$')

    # Threshold line
    ax.axvline(x_crit, color='#888888', lw=1.2, ls=':', zorder=2)

    # Axis formatting
    ax.set_xlabel(xlabel, fontsize=20, labelpad=8)
    ax.spines['top'].set_visible(False)
    ax_r.spines['top'].set_visible(False)
    ax.tick_params(labelsize=15)
    ax_r.tick_params(labelsize=15)

    if show_ylabel_left:
        ax.set_ylabel(r'$\int [R(x)]_+\,dx$', fontsize=20, labelpad=8,
                      color=AMBER_FILL)
    ax.tick_params(axis='y', colors=AMBER_FILL)
    ax.spines['left'].set_edgecolor(AMBER_FILL)

    if show_ylabel_right:
        ax_r.set_ylabel(r'$\max_p\, y^*\ $ and$\ \max_p\, \Lambda_{\rm inv}$',
                        fontsize=20, labelpad=8, color=BLUE)
    ax_r.tick_params(axis='y', colors=BLUE)
    ax_r.spines['right'].set_edgecolor(BLUE)

    # Tidy y limits — floor at 0
    ax.set_ylim(bottom=0)
    ax_r.set_ylim(bottom=0)

    # Shared legend
    lines = [l1, l2, l3]
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, fontsize=20, frameon=True, framealpha=0.92,
              edgecolor='#dddddd', loc='upper right')

    # ── Inset: zoom near threshold, twin-axes matching the main panel ─────
    inset = ax.inset_axes([0.54, 0.06, 0.44, 0.44])  # [x0, y0, w, h]
    inset_r = inset.twinx()

    mask = (xvals >= inset_xlim[0]) & (xvals <= inset_xlim[1])
    xv = xvals[mask]

    inset.fill_between(xv, areas[mask], color='#f5e6a3', alpha=0.65, zorder=0)
    inset.plot(xv, areas[mask],   color=AMBER_FILL, lw=1.8, zorder=4)
    inset_r.plot(xv, ystars[mask],  color=BLUE,       lw=1.8, zorder=5)
    inset_r.plot(xv, lambdas[mask], color=CRIMSON,    lw=1.8, zorder=5)
    inset.axvline(x_crit, color='#888888', lw=1.0, ls=':', zorder=2)

    inset.set_xlim(inset_xlim)
    inset.set_ylim(bottom=0)
    inset_r.set_ylim(bottom=0)

    # Set tight ylims from data in window
    a_max = areas[mask].max()
    yr_max = max(ystars[mask].max(), lambdas[mask].max())
    inset.set_ylim(-a_max * 0.08, a_max * 1.15)
    inset_r.set_ylim(-yr_max * 0.08, yr_max * 1.15)

    inset.set_xticks([])
    inset.set_yticks([])
    inset_r.set_yticks([])

    ax.indicate_inset_zoom(inset, edgecolor='#aaaaaa', lw=0.8, alpha=0.6)

    return ax_r


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main(outpath='src/general/plots/threshold/threshold_quantities.pdf'):
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    matplotlib.rcParams.update({
        'font.family':       'sans-serif',
        'font.sans-serif':   ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size':         18,
        'axes.linewidth':    1.2,
        'xtick.major.width': 1.1,
        'ytick.major.width': 1.1,
        'xtick.major.size':  5,
        'ytick.major.size':  5,
    })

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8), layout='constrained')
    fig.patch.set_facecolor('white')

    # ── Left: Base LV vs gamma_2 ─────────────────────────────────────────────
    g2_crit = 3.0
    g2_vals = np.linspace(2.0, 3.5, 1000)
    areas, ystars, lambdas = base_quantities(g2_vals)

    draw_panel(
        axes[0], g2_vals, areas, ystars, lambdas,
        x_crit   = g2_crit,
        xlabel   = r'Predator death rate  $\gamma_2$',
        inset_xlim = (g2_crit * 0.9, g2_crit * 1.05),
        show_ylabel_left=True, show_ylabel_right=True,
    )

    # ── Right: Shielding vs alpha_B ──────────────────────────────────────────
    aB_crit = np.log(2)
    aB_vals = np.linspace(1/np.e, 1.4, 1000)
    areas2, ystars2, lambdas2 = shield_quantities(aB_vals)

    draw_panel(
        axes[1], aB_vals, areas2, ystars2, lambdas2,
        x_crit   = aB_crit,
        xlabel   = r'Shielding strength  $\alpha_B$',
        inset_xlim = (aB_crit * 0.9, aB_crit * 1.05),
        show_ylabel_left=True, show_ylabel_right=True,
    )

    fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
    svg_path = outpath.replace('.pdf', '.svg')
    fig.savefig(svg_path, dpi=300, facecolor='white', edgecolor='none')
    print(f"Saved {outpath} and {svg_path}")


if __name__ == '__main__':
    main()