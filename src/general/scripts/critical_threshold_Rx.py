"""
Three-panel figure showing how multiple quantities all vanish at the same critical threshold.

Left:   Base LV model   — quantities vs gamma_2  (threshold gamma_2* ≈ 2.31)
Middle: Shielding model — quantities vs alpha_2  (threshold alpha_2* = ln 2 ≈ 0.693)
Right:  Viral model     — quantities vs delta_2  (threshold delta_2* ≈ 2.70)

Base LV fixed params:  r1=K1=1, gamma1=1.1, r2=4, K2=2.

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
import os

# ─────────────────────────────────────────────────────────────────────────────
# Base LV parameters
r1, K1, g1 = 1.0, 1.0, 1.1
r2, K2     = 4.0, 2.0
c1, c2     = r1 / K1, r2 / K2

AMBER_LINE = '#E0C233'   # gold line (lightened)
AMBER_FILL = '#faefc4'   # pale gold fill (lightened)
BLUE       = '#547AA5'   # steel blue (lightened, matches limit_cycles_Rx.py)
CRIMSON    = '#b84042'   # crimson (lightened, matches limit_cycles_Rx.py)
THRESH_COL = '#666666'

# Canonical figure geometry (shared with limit_cycles_Rx.py so the two figures
# stack cleanly at the same \linewidth).
FIGSIZE   = (21, 5.6)
FS_TITLE  = 20
FS_LABEL  = 20
FS_TICK   = 15
FS_LEGEND = 20


# ─────────────────────────────────────────────────────────────────────────────
# BASE LV
# ─────────────────────────────────────────────────────────────────────────────

def R_base(x, g2):
    f1 = r1 * x * (1 - x / K1)
    f2 = r2 * x * (1 - x / K2)
    return (x - g2) / f2 - (x - g1) / f1

def base_quantities(g2_vals):
    areas, ystars, lambdas = [], [], []
    p = np.linspace(0.01, 0.99, 500)
    for g2 in g2_vals:
        xs = np.linspace(K1 + 1e-6, K2 - 1e-6, 4000)
        areas.append(np.trapezoid(np.maximum(R_base(xs, g2), 0.0), xs))

        r_bar = p * r1 + (1 - p) * r2
        c_bar = p * c1 + (1 - p) * c2
        g_bar = p * g1 + (1 - p) * g2
        ystars.append(max(np.max(r_bar - c_bar * g_bar), 0.0))   # max_p y*, capped
        lambdas.append(np.max(r_bar / c_bar - g_bar))             # max_p lambda

    return np.array(areas), np.array(ystars), np.array(lambdas)


# ─────────────────────────────────────────────────────────────────────────────
# SHIELDING MODEL  (predator death fixed at 1; env A has mu=0, env B mu=alpha_B)
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
        xs = np.linspace(KA + 1e-6, KB - 1e-6, 4000)
        r_vals = np.array([R_shield(x, aB) for x in xs])
        area = np.trapezoid(np.maximum(r_vals, 0.0), xs)
        areas.append(area)

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

        def lambda_inv_shield(p):
            x0 = 1.0 / (p/KA + (1-p)/KB)
            return p*(x0 - g) + (1-p)*(x0*np.exp(-aB*x0) - g)

        lam_vals = [lambda_inv_shield(p) for p in np.linspace(0.01, 0.99, 500)]
        lambdas.append(max(lam_vals))

    return np.array(areas), np.array(ystars), np.array(lambdas)


# ─────────────────────────────────────────────────────────────────────────────
# VIRAL MODEL  (large-c quasi-steady limit; x = T uninfected target cells)
#   dT/dt = lam - d T - a T I ,   dI/dt = a T I - delta I ,   a = beta p / c
#   so on the I=0 boundary  f(T) = lam - d T   and   h(T) = a T - delta,
#   with carrying capacity  K = lam / d.
# env1 (slow): d1=lam1=1 -> K1=1 ;  env2 (fast): d2=4, lam2=8 -> K2=2.
# ─────────────────────────────────────────────────────────────────────────────

vd1, vlam1, va1, vde1 = 1.0, 1.0, 1.0, 1.2
vd2, vlam2, va2       = 4.0, 8.0, 1.0          # delta_2 is swept
vK1, vK2 = vlam1 / vd1, vlam2 / vd2            # = 1, 2

def R_viral(x, de2):
    f1 = vlam1 - vd1 * x
    f2 = vlam2 - vd2 * x
    return (va2 * x - de2) / f2 - (va1 * x - vde1) / f1

def viral_quantities(de2_vals):
    areas, ystars, lambdas = [], [], []
    p = np.linspace(0.01, 0.99, 500)
    for de2 in de2_vals:
        Ts = np.linspace(vK1 + 1e-6, vK2 - 1e-6, 4000)
        areas.append(np.trapezoid(np.maximum(R_viral(Ts, de2), 0.0), Ts))

        a_bar  = p * va1   + (1 - p) * va2
        l_bar  = p * vlam1 + (1 - p) * vlam2
        d_bar  = p * vd1   + (1 - p) * vd2
        e_bar  = p * vde1  + (1 - p) * de2
        Istar  = l_bar / e_bar - d_bar / a_bar          # infected steady state
        ystars.append(max(np.max(np.clip(Istar, 0.0, None)), 0.0))
        lambdas.append(np.max(a_bar * l_bar / d_bar - e_bar))   # invasion exponent

    return np.array(areas), np.array(ystars), np.array(lambdas)


# ─────────────────────────────────────────────────────────────────────────────
# PANEL DRAWING
# ─────────────────────────────────────────────────────────────────────────────

def draw_panel(ax, xvals, areas, ystars, lambdas, x_crit, xlabel, title,
               inset_xlim, thresh_label,
               show_ylabel_left=True, show_ylabel_right=True, show_legend=False):
    """
    Main axes:  left y-axis  → amber R(x)+ area
                right y-axis → blue y*, crimson Lambda_inv
    Inset: same twin-axis layout, same left/right axis ratio as main panel.
    """

    # ── Twin axes ─────────────────────────────────────────────────────────
    ax_r = ax.twinx()

    ax.fill_between(xvals, areas, color=AMBER_FILL, alpha=0.75, zorder=0)
    l1, = ax.plot(xvals, areas, color=AMBER_LINE, lw=2.6, zorder=4,
                  label=r'$\int_{K_1}^{K_2} [R(x)]_+\,dx$')
    l2, = ax_r.plot(xvals, ystars,  color=BLUE,    lw=2.4, zorder=5,
                    label=r'$\max_p\, y^*$')
    l3, = ax_r.plot(xvals, lambdas, color=CRIMSON, lw=2.4, zorder=5, ls='-',
                    label=r'$\max_p\, \Lambda_{\rm inv}$')

    # 1. Solid threshold line (was dotted)
    ax.axvline(x_crit, color=THRESH_COL, lw=1.4, ls='-', zorder=2)

    ax.set_title(title, fontsize=FS_TITLE, pad=10)
    ax.set_xlabel(xlabel, fontsize=FS_LABEL, labelpad=8)
    ax.spines['top'].set_visible(False)
    ax_r.spines['top'].set_visible(False)
    ax.tick_params(labelsize=FS_TICK)
    ax_r.tick_params(labelsize=FS_TICK)

    if show_ylabel_left:
        ax.set_ylabel(r'$\int [R(x)]_+\,dx$', fontsize=FS_LABEL, labelpad=8,
                      color=AMBER_LINE)
    ax.tick_params(axis='y', colors=AMBER_LINE)
    ax.spines['left'].set_edgecolor(AMBER_LINE)

    if show_ylabel_right:
        ax_r.set_ylabel(r'$\max_p\, y^*\ $ and$\ \max_p\, \Lambda_{\rm inv}$',
                        fontsize=FS_LABEL, labelpad=8, color=BLUE)
    ax_r.tick_params(axis='y', colors=BLUE)
    ax_r.spines['right'].set_edgecolor(BLUE)

    ax.set_ylim(bottom=0)
    ax_r.set_ylim(bottom=0)

    # Record the axis limits that matplotlib chose for the full data
    ax.figure.canvas.draw()
    main_left_top  = ax.get_ylim()[1]
    main_right_top = ax_r.get_ylim()[1]
    ratio = main_right_top / main_left_top

    # 2. Threshold label on main panel — placed just to the left of the line,
    ax.text(x_crit - (xvals[-1] - xvals[0]) * 0.015, main_left_top * 0.5,
            thresh_label, ha='right', va='top', fontsize=18,
            color=THRESH_COL, zorder=6)

    if show_legend:
        lines = [l1, l2, l3]
        labels = [l.get_label() for l in lines]
        ax.legend(lines, labels, fontsize=FS_LEGEND, frameon=True, framealpha=0.92,
                  edgecolor='#dddddd', loc='upper right')

    # ── Inset ─────────────────────────────────────────────────────────────
    inset = ax.inset_axes([0.54, 0.06, 0.44, 0.44])
    inset_r = inset.twinx()

    mask = (xvals >= inset_xlim[0]) & (xvals <= inset_xlim[1])
    xv = xvals[mask]

    inset.fill_between(xv, areas[mask], color=AMBER_FILL, alpha=0.65, zorder=0)
    inset.plot(xv, areas[mask],   color=AMBER_LINE, lw=1.8, zorder=4)
    inset_r.plot(xv, ystars[mask],  color=BLUE,    lw=1.8, zorder=5)
    inset_r.plot(xv, lambdas[mask], color=CRIMSON, lw=1.8, zorder=5)

    # 1. Solid threshold line in inset
    inset.axvline(x_crit, color=THRESH_COL, lw=1.0, ls='-', zorder=2)

    inset.set_xlim(inset_xlim)

    a_max_inset  = areas[mask].max()
    left_top  = a_max_inset * 1.15
    right_top = left_top * ratio

    inset.set_ylim(-0.05 * left_top, left_top)
    inset_r.set_ylim(-0.05 * right_top, right_top)

    inset.set_xticks([])
    inset.set_yticks([])
    inset_r.set_yticks([])

    # 1. Add top spine to inset to close the box
    inset.spines['top'].set_visible(True)
    inset.spines['top'].set_color('black')
    inset.spines['top'].set_linewidth(0.8)
    inset_r.spines['top'].set_visible(False)   # twin's top would double-draw

    # # 2. Threshold label inside inset — placed just left of line near top
    # inset_x_span = inset_xlim[1] - inset_xlim[0]
    # inset.text(x_crit - inset_x_span * 0.03, left_top * 0.95,
    #            thresh_label, ha='right', va='top', fontsize=16,
    #            color=THRESH_COL, zorder=6)

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

    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE, layout='constrained')
    fig.patch.set_facecolor('white')

    # ── Left: Base LV vs gamma_2 ─────────────────────────────────────────────
    # Sweep only the physical regime gamma_2 >= K2 (both environments exclude
    # the predator). Below K2 the R(x)+ integral is improper (R ~ 1/(K2-x) at
    # the boundary) and diverges, which is unphysical here.
    g2_crit = 2.306
    g2_vals = np.linspace(2.0, 2.6, 1000)
    areas, ystars, lambdas = base_quantities(g2_vals)

    draw_panel(
        axes[0], g2_vals, areas, ystars, lambdas,
        x_crit          = g2_crit,
        xlabel          = r'Predator death rate  $\gamma_2$',
        title           = 'Lotka-Volterra',
        inset_xlim      = (g2_crit * 0.96, g2_crit * 1.02),
        thresh_label    = r'$\gamma_2\!\approx\!2.31$',
        show_ylabel_left=True, show_ylabel_right=False, show_legend=True,
    )

    # ── Middle: Shielding vs alpha_B ─────────────────────────────────────────
    aB_crit = np.log(2)
    aB_vals = np.linspace(1/np.e, 1.2, 1000)
    areas2, ystars2, lambdas2 = shield_quantities(aB_vals)

    draw_panel(
        axes[1], aB_vals, areas2, ystars2, lambdas2,
        x_crit          = aB_crit,
        xlabel          = r'Shielding strength  $\mu_2$',
        title           = 'Shielding',
        inset_xlim      = (aB_crit * 0.9, aB_crit * 1.05),
        thresh_label    = r'$\mu_2\!=\!\ln 2$',
        show_ylabel_left=False, show_ylabel_right=False,
    )

    # ── Right: Viral model vs delta_2 ────────────────────────────────────────
    de2_crit = 2.701
    de2_vals = np.linspace(2.0, 3.4, 1000)   # physical regime delta_2 >= a2*K2
    areas3, ystars3, lambdas3 = viral_quantities(de2_vals)

    draw_panel(
        axes[2], de2_vals, areas3, ystars3, lambdas3,
        x_crit          = de2_crit,
        xlabel          = r'Infected-cell death rate  $\delta_2$',
        title           = 'Viral infection',
        inset_xlim      = (de2_crit * 0.96, de2_crit * 1.02),
        thresh_label    = r'$\delta_2\!\approx\!2.70$',
        show_ylabel_left=False, show_ylabel_right=True,
    )

    fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
    svg_path = outpath.replace('.pdf', '.svg')
    fig.savefig(svg_path, dpi=300, facecolor='white', edgecolor='none')
    print(f"Saved {outpath} and {svg_path}")


if __name__ == '__main__':
    main()
