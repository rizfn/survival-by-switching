"""
Two-panel figure showing how multiple quantities all vanish at the same critical threshold.

Left:  Base LV model — quantities vs gamma_2  (threshold gamma_2* = 3)
Right: Shielding model — quantities vs alpha_2 (threshold alpha_2* = ln 2 ≈ 0.693)

Fixed params:  r1=K1=gamma1=1, r2=4, K2=2.

Main panel:
  - Amber vertical shading where R(x) > 0 exists in (K1, K2)  [binary: yes/no]
  - Blue  y* on left y-axis
  - Crimson Lambda_inv on right y-axis (shared with blue via twinx)

Inset: same design, zoomed near threshold.
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
AMBER_BG   = '#f5e6a3'
BLUE       = '#2d4f73'
CRIMSON    = '#901A1E'

# ─────────────────────────────────────────────────────────────────────────────
# BASE LV
# ─────────────────────────────────────────────────────────────────────────────

def R_base(x, g2):
    h2_f2 = (x - g2) / (2 * x * (2 - x))
    h1_f1 = -1.0 / x
    return h2_f2 - h1_f1

def base_R_positive(g2):
    """True if R(x) > 0 anywhere in (K1, K2)."""
    xz = 4.0 - g2          # analytic zero of R(x)
    return xz > K1          # R > 0 on (K1, xz) when xz > K1

def base_quantities(g2_vals):
    ystars, lambdas, r_positive = [], [], []

    for g2 in g2_vals:
        r_positive.append(base_R_positive(g2))

        # y* maximised over p
        if abs(g2 - 1.0) < 1e-9:
            p_opt = 0.5
        else:
            p_opt = (5 - 3*g2) / (2 - 2*g2)

        def ystar(p):
            return max((1 - p) * (4 - p - g2 * (2 - p)), 0.0)

        p_clamped = max(1e-4, min(0.9999, p_opt))
        ystars.append(max(ystar(p_clamped), ystar(1e-4), ystar(0.9999)))

        # Invasion exponent maximised over p
        def x0_base(p):
            return (p*r1 + (1-p)*r2) / (p*r1/K1 + (1-p)*r2/K2)

        def lambda_inv_base(p):
            x0 = x0_base(p)
            return p*(x0 - g1) + (1-p)*(x0 - g2)

        lam_vals = [lambda_inv_base(p) for p in np.linspace(0.01, 0.99, 500)]
        lambdas.append(max(lam_vals))

    return np.array(ystars), np.array(lambdas), np.array(r_positive)


# ─────────────────────────────────────────────────────────────────────────────
# SHIELDING MODEL
# ─────────────────────────────────────────────────────────────────────────────

def R_shield(x, aB, g=1.0, KB=2.0):
    KA = 1.0
    h_A_f_A = -1.0 / x
    h_B = x * np.exp(-aB * x) - g
    f_B = x * (1 - x / KB)
    return h_B / f_B - h_A_f_A

def shield_R_positive(aB, KA=1.0, KB=2.0):
    """True if R(x) > 0 anywhere in (KA, KB)."""
    if aB < 1e-9:
        return True
    xz = np.log(KB) / aB   # analytic zero
    return xz > KA

def shield_quantities(aB_vals):
    KA, KB, g = 1.0, 2.0, 1.0
    ystars, lambdas, r_positive = [], [], []

    for aB in aB_vals:
        r_positive.append(shield_R_positive(aB))

        # y* parametrised by x* in (KA, KB)
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

        # Invasion exponent maximised over p
        def lambda_inv_shield(p):
            x0 = 1.0 / (p/KA + (1-p)/KB)
            return p*(x0 - g) + (1-p)*(x0*np.exp(-aB*x0) - g)

        lam_vals = [lambda_inv_shield(p) for p in np.linspace(0.01, 0.99, 500)]
        lambdas.append(max(lam_vals))

    return np.array(ystars), np.array(lambdas), np.array(r_positive)


# ─────────────────────────────────────────────────────────────────────────────
# SHARED SHADING HELPER
# ─────────────────────────────────────────────────────────────────────────────

def shade_R_positive(ax, xvals, r_positive, ymin=0, ymax=1,
                     color=AMBER_BG, alpha=0.75, zorder=0):
    """
    Fill vertical bands where r_positive is True.
    Finds contiguous True-runs and calls axvspan for each.
    """
    dx = xvals[1] - xvals[0]
    in_band = False
    x_start = None
    for i, (x, pos) in enumerate(zip(xvals, r_positive)):
        if pos and not in_band:
            x_start = x - dx/2
            in_band = True
        elif not pos and in_band:
            ax.axvspan(x_start, x - dx/2, ymin=ymin, ymax=ymax,
                       color=color, alpha=alpha, zorder=zorder, lw=0)
            in_band = False
    if in_band:
        ax.axvspan(x_start, xvals[-1] + dx/2, ymin=ymin, ymax=ymax,
                   color=color, alpha=alpha, zorder=zorder, lw=0)

    # Draw a single amber edge line at the right boundary (threshold)
    # Find rightmost transition True->False
    transitions = [xvals[i] for i in range(1, len(r_positive))
                   if r_positive[i-1] and not r_positive[i]]
    for xt in transitions:
        ax.axvline(xt, color=AMBER_FILL, lw=1.6, ls=(0,(5,3)),
                   zorder=zorder+1, alpha=0.9)


# ─────────────────────────────────────────────────────────────────────────────
# PANEL DRAWING
# ─────────────────────────────────────────────────────────────────────────────

def draw_panel(ax, xvals, ystars, lambdas, r_positive,
               x_crit, xlabel, inset_xlim,
               show_ylabel_left=True, show_ylabel_right=True):

    # ── Twin axes: blue y* left, crimson Lambda right ─────────────────────
    ax_r = ax.twinx()

    # Amber binary shading first (background)
    shade_R_positive(ax, xvals, r_positive, zorder=0)

    # Curves
    l1, = ax.plot(xvals, ystars,  color=BLUE,    lw=2.4, zorder=5,
                  label=r'$\max_p\, y^*$')
    l2, = ax_r.plot(xvals, lambdas, color=CRIMSON, lw=2.4, zorder=5,
                    label=r'$\max_p\, \Lambda_{\rm inv}$')

    # Threshold line
    ax.axvline(x_crit, color='#888888', lw=1.2, ls=':', zorder=3)

    # Axis formatting
    ax.set_xlabel(xlabel, fontsize=20, labelpad=8)
    ax.spines['top'].set_visible(False)
    ax_r.spines['top'].set_visible(False)
    ax.tick_params(labelsize=15)
    ax_r.tick_params(labelsize=15)

    if show_ylabel_left:
        ax.set_ylabel(r'$\max_p\, y^*$', fontsize=20, labelpad=8, color=BLUE)
    ax.tick_params(axis='y', colors=BLUE)
    ax.spines['left'].set_edgecolor(BLUE)

    if show_ylabel_right:
        ax_r.set_ylabel(r'$\max_p\, \Lambda_{\rm inv}$',
                        fontsize=20, labelpad=8, color=CRIMSON)
    ax_r.tick_params(axis='y', colors=CRIMSON)
    ax_r.spines['right'].set_edgecolor(CRIMSON)

    ax.set_ylim(bottom=0)
    ax_r.set_ylim(bottom=0)

    # Legend: add amber patch manually
    from matplotlib.patches import Patch
    amber_patch = Patch(facecolor=AMBER_BG, edgecolor=AMBER_FILL,
                        label=r'$R(x) > 0$ in $(K_1, K_2)$',
                        linewidth=1.4, linestyle=(0,(5,3)))
    lines = [amber_patch, l1, l2]
    labels = [h.get_label() for h in lines]
    ax.legend(lines, labels, fontsize=16, frameon=True, framealpha=0.92,
              edgecolor='#dddddd', loc='upper right')

    # ── Inset ─────────────────────────────────────────────────────────────
    inset = ax.inset_axes([0.54, 0.08, 0.43, 0.43])
    inset_r = inset.twinx()

    mask = (xvals >= inset_xlim[0]) & (xvals <= inset_xlim[1])
    xv   = xvals[mask]

    shade_R_positive(inset, xv, r_positive[mask], zorder=0)
    inset.plot(xv, ystars[mask],   color=BLUE,    lw=1.8, zorder=5)
    inset_r.plot(xv, lambdas[mask], color=CRIMSON, lw=1.8, zorder=5)
    inset.axvline(x_crit, color='#888888', lw=1.0, ls=':', zorder=3)

    inset.set_xlim(inset_xlim)

    y_max  = ystars[mask].max()
    lr_max = lambdas[mask].max()
    inset.set_ylim(-y_max  * 0.08, y_max  * 1.18)
    inset_r.set_ylim(-lr_max * 0.08, lr_max * 1.18)

    inset.tick_params(labelsize=9, colors=BLUE)
    inset_r.tick_params(labelsize=9, colors=CRIMSON)
    inset.spines['top'].set_visible(False)
    inset_r.spines['top'].set_visible(False)
    inset.spines['left'].set_edgecolor(BLUE)
    inset_r.spines['right'].set_edgecolor(CRIMSON)
    inset.yaxis.set_major_locator(MaxNLocator(nbins=3, prune='upper'))
    inset_r.yaxis.set_major_locator(MaxNLocator(nbins=3, prune='upper'))
    inset.xaxis.set_major_locator(MaxNLocator(nbins=3))
    inset.set_title('near threshold', fontsize=9, pad=3, color='#444444')

    ax.indicate_inset_zoom(inset, edgecolor='#aaaaaa', lw=0.8, alpha=0.6)

    return ax_r


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main(outpath='src/general/plots/threshold/discreteR_threshold_quantities.pdf'):
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
    ystars, lambdas, r_pos = base_quantities(g2_vals)

    draw_panel(
        axes[0], g2_vals, ystars, lambdas, r_pos,
        x_crit     = g2_crit,
        xlabel     = r'Predator death rate  $\gamma_2$',
        inset_xlim = (g2_crit * 0.9, g2_crit * 1.05),
        show_ylabel_left=True, show_ylabel_right=True,
    )

    # ── Right: Shielding vs alpha_B ──────────────────────────────────────────
    aB_crit = np.log(2)
    aB_vals = np.linspace(1/np.e, 1.4, 1000)
    ystars2, lambdas2, r_pos2 = shield_quantities(aB_vals)

    draw_panel(
        axes[1], aB_vals, ystars2, lambdas2, r_pos2,
        x_crit     = aB_crit,
        xlabel     = r'Shielding strength  $\alpha_B$',
        inset_xlim = (aB_crit * 0.9, aB_crit * 1.05),
        show_ylabel_left=True, show_ylabel_right=True,
    )

    fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
    svg_path = outpath.replace('.pdf', '.svg')
    fig.savefig(svg_path, dpi=300, facecolor='white', edgecolor='none')
    print(f"Saved {outpath} and {svg_path}")


if __name__ == '__main__':
    main()