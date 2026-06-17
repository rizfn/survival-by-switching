"""
Two-panel proof figure for R(x) theory (base LV model).

Parameters: r1=K1=g1=1, r2=4, K2=2, g2=2.1  (threshold g2*=3)

Panel A — Heatmap of integral_{xa}^{xb} R(x) dx over the (xa, xb) plane.
  Green  = positive net log gain (cycle spanning that range can sustain y).
  Red    = negative net log gain (y dies out).
  Overlay: turning points (xa^T, xb^T) of simulated limit cycles for
           several T values, coloured by T.  Points in the green zone
           persist; points in the red zone go extinct.

Panel B — Two curves vs gamma2, both crossing zero at gamma2*=3:
  (i)  max_p y*  (coexistence fixed-point height of averaged ODE)
  (ii) x_{R=0}  position of the R(x)=0 crossing (right axis, or normalised)
  Vertical line at gamma2*=3.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import TwoSlopeNorm
import os

# ── Parameters ────────────────────────────────────────────────────────────────
r1, K1, g1 = 1.0, 1.0, 1.0
r2, K2     = 4.0, 2.0
G2_PLOT    = 2.1          # case shown in heatmap
G2_STAR    = 3.0          # analytical threshold

T_VALUES   = [0.3, 0.5, 1.0, 2.0, 4.0, 8.0]   # T values for overlay

# ── R(x) and its integral ─────────────────────────────────────────────────────

def R_func(x, g2):
    """R(x) for base LV.  h1/f1 = -1/x (exact, removable singularity at K1)."""
    f2 = r2 * x * (1.0 - x / K2)
    with np.errstate(divide='ignore', invalid='ignore'):
        term2 = np.where(np.abs(f2) > 1e-12, (x - g2) / f2, np.nan)
    return term2 + 1.0 / x   # = h2/f2 - h1/f1,  since -h1/f1 = +1/x


def R_integral(xa, xb, g2, n=1000):
    if xa >= xb:
        return np.nan
    x = np.linspace(xa, xb, n)
    v = R_func(x, g2)
    m = np.isfinite(v)
    if m.sum() < 2:
        return np.nan
    return np.trapezoid(v[m], x[m])


# ── ODE + RK4 ─────────────────────────────────────────────────────────────────

def rhs(x, y, env, g2):
    if env == 'A':
        return r1*x*(1-x/K1) - x*y,  x*y - g1*y
    else:
        return r2*x*(1-x/K2) - x*y,  x*y - g2*y


def rk4_step(x, y, t, T, p, dt, g2):
    e = lambda tt: 'A' if (tt % T) < p*T else 'B'
    k1x, k1y = rhs(x,             y,             e(t),       g2)
    k2x, k2y = rhs(x+.5*dt*k1x,   y+.5*dt*k1y,   e(t+.5*dt), g2)
    k3x, k3y = rhs(x+.5*dt*k2x,   y+.5*dt*k2y,   e(t+.5*dt), g2)
    k4x, k4y = rhs(x+   dt*k3x,   y+   dt*k3y,   e(t+dt),    g2)
    return (max(x + (dt/6)*(k1x+2*k2x+2*k3x+k4x), 1e-14),
            max(y + (dt/6)*(k1y+2*k2y+2*k3y+k4y), 1e-14))


def turning_points(T, g2, p=0.5, dt=None, transient=3000.0, measure_periods=12):
    if dt is None:
        dt = min(5e-4, T / 30)
    x, y, t = 1.5, 0.1, 0.0
    for _ in range(int(transient / dt)):
        x, y = rk4_step(x, y, t, T, p, dt, g2)
        t += dt
    xs, ys = [], []
    n_meas = int(measure_periods * T / dt)
    for _ in range(n_meas):
        x, y = rk4_step(x, y, t, T, p, dt, g2)
        t += dt
        xs.append(x); ys.append(y)
    xs, ys = np.array(xs), np.array(ys)
    persists = np.max(ys) > 1e-4
    return np.min(xs), np.max(xs), np.max(ys), persists


# ── Panel B quantities ────────────────────────────────────────────────────────

def max_ystar(g2):
    """max over p of y* = (1-p)[4-p-g2(2-p)], with x* in (K1,K2)."""
    from scipy.optimize import minimize_scalar
    p_min = max((g2 - K2) / (g2 - 1), 0.0) if g2 > K2 else 0.0
    p_min = max(p_min, 1e-4)
    def neg_ys(p):
        return -(1-p) * (4 - p - g2*(2-p))
    res = minimize_scalar(neg_ys, bounds=(p_min, 1-1e-4), method='bounded')
    return -res.fun


# ── Build figure ──────────────────────────────────────────────────────────────

def make_figure(outpath):
    # ── 1. Heatmap grid ───────────────────────────────────────────────────────
    n = 140
    xa_arr = np.linspace(K1 + 0.01, K2 - 0.12, n)
    xb_arr = np.linspace(K1 + 0.12, K2 - 0.01, n)
    Z = np.full((n, n), np.nan)
    for i, xb in enumerate(xb_arr):
        for j, xa in enumerate(xa_arr):
            if xa < xb - 0.05:
                Z[i, j] = R_integral(xa, xb, G2_PLOT)

    # ── 2. Limit cycle turning points ────────────────────────────────────────
    print("Simulating limit cycles...")
    tp_data = []
    for T in T_VALUES:
        xa_t, xb_t, ymax, pers = turning_points(T, G2_PLOT)
        tp_data.append((T, xa_t, xb_t, ymax, pers))
        print(f"  T={T}: xa={xa_t:.4f}, xb={xb_t:.4f}, ymax={ymax:.5f}, persists={pers}")

    # ── 3. Panel B data ───────────────────────────────────────────────────────
    g2_arr    = np.linspace(1.6, 3.6, 300)
    my_arr    = np.array([max_ystar(g) for g in g2_arr])
    Rzero_arr = 4.0 - g2_arr    # analytical: R(x)=0 at x=4-g2

    # ── Plotting ──────────────────────────────────────────────────────────────
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9.5,
    })

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.patch.set_facecolor('white')

    # ═══════════════════════════════════════════════════════════════
    # Panel A: heatmap
    # ═══════════════════════════════════════════════════════════════
    ax = axes[0]

    vmax = np.nanpercentile(np.abs(Z), 97)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.pcolormesh(xa_arr, xb_arr, Z,
                       cmap='RdYlGn', norm=norm,
                       shading='auto', rasterized=True)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r'$\int_{x_a}^{x_b} R(x)\,\mathrm{d}x$', fontsize=11)
    cb.ax.tick_params(labelsize=9)

    # diagonal xa=xb guide
    lim_lo, lim_hi = K1 + 0.01, K2 - 0.01
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi],
            color='#555555', lw=0.8, ls=':', zorder=2)

    # R(x)=0 crossing: xb = 4-g2 line (constant for fixed g2)
    R_zero_x = 4.0 - G2_PLOT   # = 1.9
    ax.axhline(R_zero_x, color='#222222', lw=1.2, ls='--', zorder=3,
               label=rf'$x_{{R=0}} = {R_zero_x:.1f}$')
    ax.axvline(R_zero_x, color='#222222', lw=1.2, ls='--', zorder=3)

    # Overlay turning points
    cmap_T = plt.get_cmap('plasma')
    T_log  = np.log10(T_VALUES)
    T_norm = (T_log - T_log.min()) / (T_log.max() - T_log.min())

    for k, (T, xa_t, xb_t, ymax, pers) in enumerate(tp_data):
        colour = cmap_T(0.1 + 0.8 * T_norm[k])
        marker = 'o' if pers else 'X'
        ms     = 70  if pers else 80
        ec     = 'black'
        ax.scatter([xa_t], [xb_t], c=[colour], s=ms, marker=marker,
                   edgecolors=ec, linewidths=0.8, zorder=5,
                   label=f'$T={T}$  ' + ('(persist)' if pers else '(extinct)'))

    ax.set_xlabel(r'$x_a$  (lower turning point)', fontsize=12)
    ax.set_ylabel(r'$x_b$  (upper turning point)', fontsize=12)
    ax.set_title(rf'(a)  Net log gain $\int R\,\mathrm{{d}}x$ for $\gamma_2={G2_PLOT}$',
                 fontsize=12, pad=7)
    ax.set_xlim(K1 + 0.01, K2 - 0.01)
    ax.set_ylim(K1 + 0.01, K2 - 0.01)

    # Annotate zones
    ax.text(1.03, 1.94, 'net gain\n(survive)', color='#1a6e1a',
            fontsize=8.5, va='top', style='italic')
    ax.text(1.03, 1.65, 'net loss\n(extinct)', color='#8b0000',
            fontsize=8.5, va='top', style='italic')

    ax.legend(loc='lower right', frameon=True, framealpha=0.9,
              title=r'$T$ value', title_fontsize=8.5)
    ax.tick_params(labelsize=10)

    # ═══════════════════════════════════════════════════════════════
    # Panel B: threshold curves
    # ═══════════════════════════════════════════════════════════════
    ax2 = axes[1]

    # max_p y* (left y-axis, black)
    ax2.plot(g2_arr, my_arr, color='#1f4e79', lw=2.2,
             label=r'$\max_p\, y^*$  (averaged FP)')
    ax2.axhline(0, color='#333333', lw=0.8, ls='-')

    # x_{R=0} shifted to cross 0 at same place: plot (x_{R=0} - K1)
    ax2.plot(g2_arr, Rzero_arr - K1, color='#8b0000', lw=2.2, ls='--',
             label=r'$x_{R=0} - K_1$  (zero of $R$ relative to $K_1$)')

    # Vertical threshold line
    ax2.axvline(G2_STAR, color='#555555', lw=1.2, ls=':')
    ax2.text(G2_STAR + 0.04, ax2.get_ylim()[1] if False else 0.85,
             r'$\gamma_2^* = 3$', fontsize=10, color='#333333', va='top')

    # Mark gamma2=2.1 (the heatmap case)
    y_at_21_ys   = max_ystar(G2_PLOT)
    y_at_21_Rz   = (4.0 - G2_PLOT) - K1
    ax2.scatter([G2_PLOT], [y_at_21_ys], color='#1f4e79', s=60, zorder=5)
    ax2.scatter([G2_PLOT], [y_at_21_Rz], color='#8b0000', s=60, zorder=5,
                marker='s')
    ax2.axvline(G2_PLOT, color='#888', lw=0.8, ls='--', alpha=0.6)
    ax2.text(G2_PLOT + 0.03, -0.07,
             rf'$\gamma_2={G2_PLOT}$' + '\n(panel a)', fontsize=8.5,
             color='#555', va='top')

    ax2.set_xlabel(r'$\gamma_2$', fontsize=13)
    ax2.set_ylabel('Value', fontsize=12)
    ax2.set_title(r'(b)  Both criteria vanish at $\gamma_2^* = 3$', fontsize=12, pad=7)
    ax2.set_xlim(1.6, 3.6)
    ax2.set_ylim(-0.25, 1.05)
    ax2.legend(loc='upper right', frameon=True, framealpha=0.9)
    ax2.tick_params(labelsize=10)
    ax2.grid(True, ls='--', alpha=0.25)

    fig.tight_layout(w_pad=3.0)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=220, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {outpath}")


if __name__ == '__main__':
    make_figure('src/3paramLogisticLV/plots/critical_gamma2/proof_Rx.pdf')