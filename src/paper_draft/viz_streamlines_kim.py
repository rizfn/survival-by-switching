import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

# 'Arial' is requested for the panel labels; on systems without it, matplotlib
# silently falls back to the metric-compatible Liberation Sans.  Silence the
# (harmless) per-call fallback warnings.
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)


# ── Ecology functions ─────────────────────────────────────────────────────────

def env_params():
    return {
        'A': {'r': 1.0, 'K': 1.0, 'gamma': 1.1},
        'B': {'r': 4.0, 'K': 2.0, 'gamma': 2.1},
    }

def vector_field(x, y, env):
    r, K, g = env['r'], env['K'], env['gamma']
    return r * x * (1 - x / K) - x * y, y * (x - g)

def averaged_vector_field(x, y, envA, envB):
    dx1, dy1 = vector_field(x, y, envA)
    dx2, dy2 = vector_field(x, y, envB)
    return 0.5*(dx1+dx2), 0.5*(dy1+dy2)

def averaged_params(envA, envB):
    return {'r': 0.5*(envA['r']+envB['r']),
            'gamma': 0.5*(envA['gamma']+envB['gamma'])}

def find_fixed_points(env):
    r, K, g = env['r'], env['K'], env['gamma']
    pts = [(0., 0.), (K, 0.)]
    ystar = r * (1 - g / K)
    if ystar > 0:
        pts.append((g, ystar))
    return pts

def averaged_fixed_points(envA, envB):
    rA, KA, gA = envA['r'], envA['K'], envA['gamma']
    rB, KB, gB = envB['r'], envB['K'], envB['gamma']
    pts = [(0., 0.)]
    denom = rA/KA + rB/KB
    if denom > 0:
        pts.append(((rA+rB)/denom, 0.))
    xs = 0.5*(gA+gB)
    ys = 0.5*(rA*(1-xs/KA) + rB*(1-xs/KB))
    if ys > 0:
        pts.append((xs, ys))
    return pts

def jacobian(env, x, y):
    r, K, g = env['r'], env['K'], env['gamma']
    return np.array([[r*(1-2*x/K)-y, -x],[y, x-g]])

def averaged_jacobian(envA, envB, x, y):
    return 0.5*(jacobian(envA,x,y)+jacobian(envB,x,y))

def is_stable(jac, tol=1e-8):
    return np.all(np.real(np.linalg.eigvals(jac)) <= tol)


# ── Panel label ─────────────────────────────────────────────────────────────
# A large capital letter (Arial, not bold) in a small white box nestled snug
# into the top-left corner, sitting just inside the axes so the spine stays
# visible around it.

PANEL_LABEL_FONT = ['Arial', 'Liberation Sans', 'sans-serif']

def add_panel_label(ax, letter):
    ax.text(0.010, 0.974, letter, transform=ax.transAxes,
            ha='left', va='top', fontsize=32, fontweight='normal',
            family=PANEL_LABEL_FONT, color='#000000',
            bbox=dict(boxstyle='square,pad=0.12', facecolor='white',
                      edgecolor='none'), zorder=10)


# ── Fixed-point marker ────────────────────────────────────────────────────────
# A black core with a white halo ring, so it reads clearly against dark
# streamlines.  Deliberately the inverse of the initial-condition marker in
# panel D (white face / black rim), so the two are not confused.

def mark_fixed_point(ax, fx, fy):
    ax.plot(fx, fy, 'o', color='#666666', markersize=14, zorder=5)   # halo
    ax.plot(fx, fy, 'o', color='#000000', markersize=11, zorder=6) # core


# ── Switching trajectory integration ───────────────────────────────────────────
# Periodic switching: first half-period in env A, second half in env B.
# Integrated in (x, log y) coordinates with RK4 so that y stays strictly positive
# and can collapse smoothly towards extinction.

def _rhs(x, logy, env):
    y = np.exp(logy)
    dx = env['r'] * x * (1 - x / env['K']) - x * y
    dlogy = x - env['gamma']
    return dx, dlogy

def _active_env(t, T, envA, envB):
    phase = t % T
    return envA if phase < 0.5 * T else envB

def integrate_switching(x0, y0, T, n_periods, envA, envB, steps_per_period=400):
    dt = T / steps_per_period
    x, logy, t = x0, np.log(y0), 0.0
    xs, ys = [x], [np.exp(logy)]
    for _ in range(n_periods * steps_per_period):
        e1 = _active_env(t, T, envA, envB)
        e2 = _active_env(t + 0.5*dt, T, envA, envB)
        e4 = _active_env(t + dt, T, envA, envB)
        k1x, k1l = _rhs(x, logy, e1)
        k2x, k2l = _rhs(x + 0.5*dt*k1x, logy + 0.5*dt*k1l, e2)
        k3x, k3l = _rhs(x + 0.5*dt*k2x, logy + 0.5*dt*k2l, e2)
        k4x, k4l = _rhs(x + dt*k3x, logy + dt*k3l, e4)
        x += dt/6.0 * (k1x + 2*k2x + 2*k3x + k4x)
        logy += dt/6.0 * (k1l + 2*k2l + 2*k3l + k4l)
        x = max(x, 1e-14)
        t += dt
        xs.append(x)
        ys.append(np.exp(logy))
    return np.array(xs), np.array(ys)


def _fade_cmap(light, dark):
    """Colormap fading from a light (transient) to a dark (final) shade."""
    return LinearSegmentedColormap.from_list('fade', [light, dark])


# y-label pad: negative so the (rotated) label tucks into the unlabelled middle
# of the axis, between the top and bottom tick numbers.
YLABEL_PAD = -14


def plot_switching_pair(ax, trajs, x_lo, x_hi, y_lo, y_hi):
    """Plot several switching trajectories on one phase plane.

    Each trajectory is coloured by progression in time, fading from a light
    (transient) shade to a dark (final-state) shade.  Slow switching collapses
    onto the extinction axis y = 0 (dead); fast switching settles onto a
    positive-y limit cycle (alive).
    """
    for tr in trajs:
        xs, ys = tr['xs'], tr['ys']
        pts = np.column_stack([xs, ys]).reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        t = np.linspace(0.0, 1.0, len(segs))
        lc = LineCollection(segs, cmap=tr['cmap'], norm=plt.Normalize(0, 1),
                            linewidth=1.9, zorder=3)
        lc.set_array(t)
        ax.add_collection(lc)

    # shared start marker (all trajectories use the same initial condition):
    # white face with a black rim — distinct from the fixed-point marker.
    x0, y0 = trajs[0]['xs'][0], trajs[0]['ys'][0]
    ax.plot(x0, y0, 'o', color='#1a1a1a', markersize=9,
            markerfacecolor='white', markeredgewidth=1.8, zorder=8)

    # annotate each trajectory in its own dark shade, in DATA coordinates
    for tr in trajs:
        ax.text(tr['label_xy'][0], tr['label_xy'][1], tr['label'],
                color=tr['label_color'], ha=tr['label_ha'], va=tr['label_va'],
                fontsize=18, zorder=7)

    ax.set_xlabel('Prey  $x$', labelpad=4)
    ax.set_ylabel('Predator  $y$', labelpad=YLABEL_PAD*1.5)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.grid(False)
    ax.patch.set_facecolor('none')


# ── Main plot ─────────────────────────────────────────────────────────────────

def _param_tag(envA, envB, alpha_fast, alpha_slow):
    return (f"K1-{envA['K']:g}_r1-{envA['r']:g}_g1-{envA['gamma']:g}"
            f"_K2-{envB['K']:g}_r2-{envB['r']:g}_g2-{envB['gamma']:g}"
            f"_af-{alpha_fast:g}_as-{alpha_slow:g}")


def plot_phase_portrait(params=None, outdir=None,
                        alpha_fast=2.0, alpha_slow=0.2,
                        x0=0.3, y0=0.048,
                        n_periods_fast=40, n_periods_slow=14):
    if params is None:
        params = env_params()

    plt.rcParams.update({
        'font.size':        20,
        'axes.labelsize':   20,
        'xtick.labelsize':  16,
        'ytick.labelsize':  16,
        'axes.linewidth':   1.1,
        'xtick.direction':  'in',
        'ytick.direction':  'in',
        'xtick.major.size': 5,
        'ytick.major.size': 5,
        'pdf.fonttype':     42,
        'svg.fonttype':     'none',
    })

    envA   = params['A']
    envB   = params['B']
    envAvg = averaged_params(envA, envB)

    fps_A   = find_fixed_points(envA)
    fps_B   = find_fixed_points(envB)
    fps_avg = averaged_fixed_points(envA, envB)

    # ── Axis limits ───────────────────────────────────────────────────────────
    all_y      = [fy for (_, fy) in fps_A + fps_B + fps_avg]
    fp_y_max   = max(all_y)
    x_data_max = max(envA['K'], envB['K'],
                     max(fx for (fx, _) in fps_A + fps_B + fps_avg))
    x_lo, x_hi = 0.0, 2.5
    # reduced y-limit: common ceiling across every panel
    y_lo, y_hi = 0.0, 0.2

    # ── Grid for streamplot ───────────────────────────────────────────────────
    Ng = 400
    xg = np.linspace(x_lo + 1e-3, x_hi, Ng)
    yg = np.linspace(y_lo + 1e-3, y_hi, Ng)
    Xg, Yg = np.meshgrid(xg, yg)

    def make_UV(vf_fn):
        U, V = vf_fn(Xg, Yg)
        return U, V, np.hypot(U, V)

    UA,  VA,  spA  = make_UV(lambda x, y: vector_field(x, y, envA))
    UB,  VB,  spB  = make_UV(lambda x, y: vector_field(x, y, envB))
    Uav, Vav, spAv = make_UV(lambda x, y: averaged_vector_field(x, y, envA, envB))

    sp_max = max(spA.max(), spB.max(), spAv.max())
    norm   = plt.Normalize(vmin=0.0, vmax=sp_max)

    # ── Switching trajectories ─────────────────────────────────────────────────
    # The switching rate alpha sets the dwell time 1/alpha in each environment,
    # so one full A→B period lasts 2/alpha.
    steps_per_period = 400
    period_fast = 2.0 / alpha_fast
    period_slow = 2.0 / alpha_slow
    xs_fast, ys_fast = integrate_switching(x0, y0, period_fast, n_periods_fast,
                                           envA, envB, steps_per_period)
    xs_slow, ys_slow = integrate_switching(x0, y0, period_slow, n_periods_slow,
                                           envA, envB, steps_per_period)

    # the switching panel gets its own (tighter) ceiling: the dynamics — start,
    # limit cycle, decay to extinction — all live in a low band, so zoom in on
    # it.  This panel does NOT share the y-axis with the streamline panels.
    y_lo_traj = -0.004
    y_hi_traj = max(ys_fast.max(), ys_slow.max(), y0) * 1.05

    # ── Style ─────────────────────────────────────────────────────────────────
    STREAM_CMAP = 'plasma'   # dark-purple → magenta → bright orange; no harsh yellow

    # fast switching → red (light transient → dark final limit cycle)
    # slow switching → blue (light transient → dark final dead state)
    CMAP_FAST = _fade_cmap('#f5a3a3', '#8b0000')
    CMAP_SLOW = _fade_cmap('#9ecae1', '#08306b')
    C_FAST_DARK = '#8b0000'
    C_SLOW_DARK = '#08306b'

    # ── Layout:  4 × 1  (single stacked column) ─────────────────────────────────
    #   Env1 (A) / Env2 (B) / Avg (C) / Switching (D)
    fig = plt.figure(figsize=(8.0, 12.0))
    fig.patch.set_facecolor('none')
    gs = GridSpec(4, 1, figure=fig,
                  hspace=0.14,
                  left=0.14, right=0.965, top=0.985, bottom=0.055)

    axA = fig.add_subplot(gs[0, 0])   # Environment 1
    axB = fig.add_subplot(gs[1, 0])   # Environment 2
    axC = fig.add_subplot(gs[2, 0])   # Time-averaged field
    axS = fig.add_subplot(gs[3, 0])   # Switching trajectories

    # streamline panels share the same x-range; only the bottom-most (axS) shows
    # the x-tick labels.
    panels = [
        dict(ax=axA, U=UA,  V=VA,  sp=spA,  fps=fps_A, letter='A',
             show_xlabel=False, jac_fn=lambda x, y: jacobian(envA, x, y)),
        dict(ax=axB, U=UB,  V=VB,  sp=spB,  fps=fps_B, letter='B',
             show_xlabel=False, jac_fn=lambda x, y: jacobian(envB, x, y)),
        dict(ax=axC, U=Uav, V=Vav, sp=spAv, fps=fps_avg, letter='C',
             show_xlabel=False, jac_fn=lambda x, y: averaged_jacobian(envA, envB, x, y)),
    ]

    for p in panels:
        ax = p['ax']
        ax.streamplot(
            xg, yg, p['U'], p['V'],
            color=p['sp'], cmap=STREAM_CMAP, norm=norm,
            linewidth=1.1, density=1.1, arrowsize=1.4,
            arrowstyle='->', minlength=0.05, zorder=1,
        )
        for (fx, fy) in p['fps']:
            if is_stable(p['jac_fn'](fx, fy)):
                mark_fixed_point(ax, fx, fy)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        # keep tick MARKERS at the intermediate values, but only LABEL the ends
        ax.set_yticks([0.0, 0.05, 0.10, 0.15, 0.20])
        ax.set_yticklabels(['0.0', '', '', '', '0.2'])
        if p['show_xlabel']:
            ax.set_xlabel('Prey  $x$', labelpad=4)
        else:
            ax.tick_params(labelbottom=False)
        ax.set_ylabel('Predator  $y$', labelpad=YLABEL_PAD)
        ax.grid(False)
        ax.patch.set_facecolor('none')
        add_panel_label(ax, p['letter'])

    # ── Switching panel: fast (red, survives) + slow (blue, dies) combined ──────
    trajs = [
        dict(xs=xs_fast, ys=ys_fast, cmap=CMAP_FAST,
             label='Fast switching', label_color=C_FAST_DARK,
             label_xy=(0.7, 0.03), label_ha='left', label_va='center'),
        dict(xs=xs_slow, ys=ys_slow, cmap=CMAP_SLOW,
             label='Slow switching', label_color=C_SLOW_DARK,
             label_xy=(0.28, 0.0), label_ha='left', label_va='bottom'),
    ]
    plot_switching_pair(axS, trajs, x_lo, x_hi, y_lo_traj, y_hi_traj)
    # keep tick MARKERS at 0.01–0.04, but only LABEL 0.00 and 0.05
    axS.set_yticks([0.00, 0.01, 0.02, 0.03, 0.04, 0.05])
    axS.set_yticklabels(['0.00', '', '', '', '', '0.05'])
    add_panel_label(axS, 'D')

    # ── Save ──────────────────────────────────────────────────────────────────
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f'streamlines.svg'
    plt.savefig(outpath, dpi=300, format='svg',
                bbox_inches='tight', facecolor='none', transparent=True)
    plt.savefig(outpath.with_suffix('.pdf'), dpi=300, format='pdf',
                bbox_inches='tight', facecolor='none', transparent=True)
    plt.savefig(outpath.with_suffix('.png'), dpi=200, format='png',
                bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {outpath}")
    return outpath


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    plot_phase_portrait(env_params(), outdir=here / 'plots' / 'streamlines_kim')