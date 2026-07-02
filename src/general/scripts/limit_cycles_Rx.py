"""
Three-panel phase-portrait figure demonstrating R(x) theory.
Base LV: r1=K1=1, g1=1.1, r2=4, K2=2.

Switching is characterised by its rate alpha (period 2/alpha; equal dwell
1/alpha per environment at p=1/2). Three panels for three gamma2 values
straddling the rescue threshold (gamma2* ≈ 2.31). Plotting order: blues
alpha=0.5..4, then reds alpha=0.5..4 (small-alpha / slow curves underneath).
Extinct attractors shown as lines at y=0. Shared y-scale derived from
simulated data. Constrained layout.

Geometry/fontsizes match critical_threshold_Rx.py so the two figures stack
cleanly at the same \\linewidth.
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

# ── Parameters ────────────────────────────────────────────────────────────────
r1, K1, g1 = 1.0, 1.0, 1.1
r2, K2     = 4.0, 2.0
# Switching rates alpha (fast → slow); the simulation uses the period T = 2/alpha.
ALPHA_VALUES = [4.0, 2.0, 1.0, 0.5]            # index 0=alpha=4 (fast) .. 3=alpha=0.5 (slow)
T_VALUES     = [2.0 / a for a in ALPHA_VALUES]  # periods = [0.5, 1.0, 2.0, 4.0]

# Canonical figure geometry (shared with critical_threshold_Rx.py)
FIGSIZE   = (21, 5.6)
FS_TITLE  = 20
FS_LABEL  = 20
FS_TICK   = 15
FS_LEGEND = 15

# ── Palette ───────────────────────────────────────────────────────────────────
# Steel blues anchored at #547AA5, lightening for larger alpha (faster)
BLUES = ['#b8cfe0', '#7aa3c0', '#547AA5', '#2d4f73']
# Deep crimsons anchored at #901A1E, lightening for larger alpha (faster)
REDS  = ['#d9888a', '#b84042', '#901A1E', '#5a0f11']

# R(x)>0 shading: golden amber (lightened), very light fill
SHADE_FC = '#faefc4'   # pale amber fill (lightened)
SHADE_EC = '#E0C233'   # amber boundary line (lightened)

# ── R(x) and its positive interval ──────────────────────────────────────────────
def R_base(x, g2):
    f1 = r1 * x * (1 - x / K1)
    f2 = r2 * x * (1 - x / K2)
    return (x - g2) / f2 - (x - g1) / f1

def Rpos_interval(g2):
    """Interval within (K1, K2) where R(x)>0 (now interior, since g1 != K1).
    Returns (x_lo, x_hi) or None if R<=0 everywhere."""
    xs = np.linspace(K1 + 1e-6, K2 - 1e-6, 40000)
    pos = xs[R_base(xs, g2) > 0]
    if len(pos) == 0:
        return None
    return pos.min(), pos.max()

# ── ODE ───────────────────────────────────────────────────────────────────────
def rhs(x, y, env, g2):
    if env == 'A': return r1*x*(1-x/K1) - x*y,  x*y - g1*y
    else:          return r2*x*(1-x/K2) - x*y,  x*y - g2*y

def rk4_step(x, y, t, T, p, dt, g2):
    e = lambda tt: 'A' if (tt % T) < p*T else 'B'
    k1x, k1y = rhs(x,            y,            e(t),        g2)
    k2x, k2y = rhs(x+.5*dt*k1x,  y+.5*dt*k1y,  e(t+.5*dt),  g2)
    k3x, k3y = rhs(x+.5*dt*k2x,  y+.5*dt*k2y,  e(t+.5*dt),  g2)
    k4x, k4y = rhs(x+   dt*k3x,  y+   dt*k3y,  e(t+dt),     g2)
    return (max(x + (dt/6)*(k1x+2*k2x+2*k3x+k4x), 1e-14),
            max(y + (dt/6)*(k1y+2*k2y+2*k3y+k4y), 1e-14))

def simulate_attractor(T, p, g2, transient=3000., measure_periods=14, dt=None):
    """Return only the limit cycle attractor (transient discarded)."""
    if dt is None:
        dt = min(4e-4, T / 40)
    x0 = p*g1 + (1-p)*g2
    y0 = max((1-p)*(4 - p - g2*(2-p)), 0.01)
    x, y, t = x0, y0, 0.0
    for _ in range(int(transient / dt)):
        x, y = rk4_step(x, y, t, T, p, dt, g2); t += dt
    xs, ys = [], []
    for _ in range(int(measure_periods * T / dt)):
        x, y = rk4_step(x, y, t, T, p, dt, g2); t += dt
        xs.append(x); ys.append(y)
    xs, ys = np.array(xs), np.array(ys)
    return xs, ys, ys.max() > 1e-4

# ── Single panel ──────────────────────────────────────────────────────────────
def plot_panel(ax, g2, label, show_ylabel, data, y_hi, show_legend=False):
    """data: dict keyed (T, p) -> (xs, ys, persists)"""
    all_xs = [K1, K2]

    # ── R(x)>0 shading (interior interval) ──────────────────────────────────────
    interval = Rpos_interval(g2)
    if interval is not None:
        xs_lo, xs_hi = interval
        ax.axvspan(xs_lo, xs_hi, color=SHADE_FC, alpha=0.85, zorder=0, linewidth=0)
        ax.axvline(xs_lo, color=SHADE_EC, lw=1.6, ls=(0, (5, 3)), zorder=2)
        ax.axvline(xs_hi, color=SHADE_EC, lw=1.6, ls=(0, (5, 3)), zorder=2)

    # ── Draw order: blues T=4..0.5, then reds T=4..0.5 ──────────────────────
    for p, colors in [(0.5, BLUES), (0.7, REDS)]:
        for k, T in reversed(list(enumerate(T_VALUES))):
            xs, ys, pers = data[(T, p)]
            col = colors[k]
            if pers:
                ax.plot(xs, ys, color=col, lw=1.8, alpha=0.92, zorder=4)
            else:
                ax.plot([xs.min(), xs.max()], [0, 0],
                        color=col, lw=2.4, alpha=0.90, zorder=3,
                        solid_capstyle='round')
            all_xs += [xs.min(), xs.max()]

    # ── K markers (dots only, no labels) ──────────────────────────────────────
    for xk in [K1, K2]:
        ax.axvline(xk, color='#cccccc', lw=0.9, ls=':', zorder=1)
        ax.scatter([xk], [0], marker='o', s=40, color='#222222',
                   zorder=7, clip_on=False)

    # ── Axes limits ───────────────────────────────────────────────────────────
    x_lo = max(min(all_xs) - 0.04, 0.88)
    x_hi = min(max(all_xs) + 0.06, K2 + 0.12)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-0.004, y_hi)

    # ── Axes decoration ───────────────────────────────────────────────────────
    ax.set_title(label, fontsize=FS_TITLE, pad=10, fontweight='normal')
    ax.set_xlabel('Prey  $x$', fontsize=FS_LABEL, labelpad=8)
    if show_ylabel:
        ax.set_ylabel('Predator  $y$', fontsize=FS_LABEL, labelpad=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=FS_TICK)

    if show_legend:
        ax.legend(handles=make_legend_handles(), loc='upper right',
                  ncol=2, frameon=True, framealpha=0.95, edgecolor='#dddddd',
                  fontsize=FS_LEGEND, handlelength=1.8, columnspacing=1.0,
                  handletextpad=0.5, labelspacing=0.3)

# ── Legend handles (two columns: blues, reds; T=4..0.5 in each) ──────────────
def make_legend_handles():
    # Blues column: alpha=0.5 (dark, index 3) .. alpha=4 (light, index 0)
    blue_handles = [
        Line2D([0],[0], color=BLUES[k], lw=2.5,
               label=rf'$\alpha={ALPHA_VALUES[k]},\ p=0.5$')
        for k, T in reversed(list(enumerate(T_VALUES)))
    ]
    # Reds column: same alpha order
    red_handles = [
        Line2D([0],[0], color=REDS[k], lw=2.5,
               label=rf'$\alpha={ALPHA_VALUES[k]},\ p=0.7$')
        for k, T in reversed(list(enumerate(T_VALUES)))
    ]
    interleaved_handles = []
    for bh, rh in zip(blue_handles, red_handles):
        interleaved_handles += [bh, rh]
    return interleaved_handles

# ── Main ──────────────────────────────────────────────────────────────────────
def main(outpath='src/general/plots/phase_space/base_LV.pdf'):
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

    G2_VALUES = [2.1, 2.25, 2.4]

    # ── Simulate everything once, store in nested dict ────────────────────────
    print("Simulating all limit cycles...")
    all_data  = {}   # all_data[g2][(T, p)] = (xs, ys, persists)
    global_ymax = 0.0
    for g2 in G2_VALUES:
        all_data[g2] = {}
        for T in T_VALUES:
            for p in [0.5, 0.7]:
                xs, ys, pers = simulate_attractor(T, p, g2)
                all_data[g2][(T, p)] = (xs, ys, pers)
                if pers:
                    global_ymax = max(global_ymax, ys.max())
                print(f"  g2={g2}, T={T}, p={p}: "
                      f"x=[{xs.min():.3f},{xs.max():.3f}], "
                      f"ymax={ys.max():.4f}, persist={pers}")

    y_hi = global_ymax * 1.20
    print(f"Shared y_hi = {y_hi:.4f}")

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE, layout='constrained')
    fig.patch.set_facecolor('white')

    configs = [
        (G2_VALUES[0], rf'(a)  $\gamma_2={G2_VALUES[0]}$'),
        (G2_VALUES[1], rf'(b)  $\gamma_2={G2_VALUES[1]}$'),
        (G2_VALUES[2], rf'(c)  $\gamma_2={G2_VALUES[2]}$'),
    ]
    for i, (ax, (g2, label), show_y) in enumerate(zip(axes, configs, [True, False, False])):
        plot_panel(ax, g2=g2, label=label, show_ylabel=show_y,
                   data=all_data[g2], y_hi=y_hi, show_legend=(i == 2))

    fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
    fig.savefig(outpath.replace('.pdf', '.svg'), dpi=300, facecolor='none', edgecolor='none')
    print(f"Saved {outpath}")


if __name__ == '__main__':
    main()