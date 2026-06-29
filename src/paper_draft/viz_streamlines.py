import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from pathlib import Path


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


def plot_trajectory(ax, xs, ys, x_lo, x_hi, y_lo, y_hi, title,
                    interior_fp=None, show_ylabel=True, cmap='inferno'):
    """Plot a phase-plane trajectory coloured by progression in time.

    Early times are dark, late times bright (yellow in 'inferno').  Because
    later segments are drawn on top, the region the trajectory keeps re-tracing
    — the attractor — is overwritten in bright yellow: a positive-y limit cycle
    when the population survives, or the extinction axis y = 0 when it dies.
    """
    pts = np.column_stack([xs, ys]).reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    t = np.linspace(0.0, 1.0, len(segs))
    lc = LineCollection(segs, cmap=cmap, norm=plt.Normalize(0, 1),
                        linewidth=1.7, zorder=3)
    lc.set_array(t)
    ax.add_collection(lc)

    # start marker
    ax.plot(xs[0], ys[0], 'o', color='#1a1a1a', markersize=8,
            markerfacecolor='white', markeredgewidth=1.7, zorder=6)

    # averaged interior fixed point (fast-switching attractor) for reference
    if interior_fp is not None:
        ax.plot(*interior_fp, marker='*', color='#1a1a1a',
                markersize=16, zorder=5, linestyle='none')

    ax.set_title(title, pad=7)
    ax.set_xlabel('Prey  $x$', labelpad=4)
    if show_ylabel:
        ax.set_ylabel('Predator  $y$', labelpad=4)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.grid(False)
    ax.patch.set_facecolor('none')


# ── Main plot ─────────────────────────────────────────────────────────────────

def plot_phase_portrait(params=None, outdir=None,
                        alpha_fast=2.0, alpha_slow=0.2,
                        x0=1.0, y0=0.05,
                        n_periods_fast=25, n_periods_slow=12):
    if params is None:
        params = env_params()

    plt.rcParams.update({
        'font.size':        20,
        'axes.titlesize':   20,
        'axes.titleweight': 'normal',
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

    # interior coexistence fixed point of the averaged field (if it exists)
    interior_fp = next(((fx, fy) for (fx, fy) in fps_avg if fy > 1e-9), None)

    # ── Axis limits ───────────────────────────────────────────────────────────
    all_y      = [fy for (_, fy) in fps_A + fps_B + fps_avg]
    fp_y_max   = max(all_y)
    x_data_max = max(envA['K'], envB['K'],
                     max(fx for (fx, _) in fps_A + fps_B + fps_avg))
    x_lo, x_hi = 0.0, x_data_max * 1.18
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

    # trajectory panels get a tighter ceiling: the dynamics (start, limit cycle,
    # decay to extinction) all live in a low band, so zoom in on it.
    y_hi_traj = max(ys_fast.max(), ys_slow.max(), y0) * 1.18

    # ── Style ─────────────────────────────────────────────────────────────────
    STREAM_CMAP = 'plasma'   # dark-purple → magenta → bright orange; no harsh yellow
    C_FP_STABLE = '#1a1a1a'

    # ── Layout:  A  B  /  C  /  D  E  ─────────────────────────────────────────
    #   3 rows × 4 cols; centre panel C spans the middle two columns so it sits
    #   centred with the same width as the other panels.
    fig = plt.figure(figsize=(11, 13.5))
    fig.patch.set_facecolor('none')
    gs = GridSpec(3, 4, figure=fig,
                  wspace=0.55, hspace=0.42,
                  left=0.10, right=0.97, top=0.95, bottom=0.06)

    axA = fig.add_subplot(gs[0, 0:2])
    axB = fig.add_subplot(gs[0, 2:4])
    axC = fig.add_subplot(gs[1, 1:3])
    axD = fig.add_subplot(gs[2, 0:2])
    axE = fig.add_subplot(gs[2, 2:4])

    panels = [
        dict(ax=axA, U=UA,  V=VA,  sp=spA,  fps=fps_A,
             title='Environment 1', show_ylabel=True,
             jac_fn=lambda x, y: jacobian(envA, x, y)),
        dict(ax=axB, U=UB,  V=VB,  sp=spB,  fps=fps_B,
             title='Environment 2', show_ylabel=False,
             jac_fn=lambda x, y: jacobian(envB, x, y)),
        dict(ax=axC, U=Uav, V=Vav, sp=spAv, fps=fps_avg,
             title='Time-averaged field ($\\alpha\\to\\infty$)', show_ylabel=True,
             jac_fn=lambda x, y: averaged_jacobian(envA, envB, x, y)),
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
                ax.plot(fx, fy, 'o', color=C_FP_STABLE,
                        markeredgecolor=C_FP_STABLE, markeredgewidth=1.6,
                        markersize=12, zorder=5)
        ax.set_title(p['title'], pad=7)
        ax.set_xlabel('Prey  $x$', labelpad=4)
        if p['show_ylabel']:
            ax.set_ylabel('Predator  $y$', labelpad=4)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.grid(False)
        ax.patch.set_facecolor('none')

    # ── Trajectory panels D (fast → limit cycle) and E (slow → extinction) ──────
    plot_trajectory(axD, xs_fast, ys_fast, x_lo, x_hi, y_lo, y_hi_traj,
                    rf'Fast switching ($\alpha={alpha_fast:g}$)',
                    interior_fp=None, show_ylabel=True, cmap='inferno')
    plot_trajectory(axE, xs_slow, ys_slow, x_lo, x_hi, y_lo, y_hi_traj,
                    rf'Slow switching ($\alpha={alpha_slow:g}$)',
                    interior_fp=None, show_ylabel=False, cmap='inferno')

    # ── Colorbar for the streamline speed, placed to the right of C ─────────────
    #   positioned from C's own bbox so C stays centred and full-width.
    sm = plt.cm.ScalarMappable(cmap=STREAM_CMAP, norm=norm)
    sm.set_array([])
    posC = axC.get_position()
    cax = fig.add_axes([posC.x1 + 0.020, posC.y0, 0.018, posC.height])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(r'Speed', labelpad=8)
    cax.patch.set_facecolor('none')

    # # ── Panel letters ───────────────────────────────────────────────────────────
    # for ax, letter in [(axA, 'A'), (axB, 'B'), (axC, 'C'),
    #                    (axD, 'D'), (axE, 'E')]:
    #     ax.text(-0.18, 1.04, letter, transform=ax.transAxes,
    #             fontsize=24, fontweight='bold', va='bottom', ha='right')

    # ── Save ──────────────────────────────────────────────────────────────────
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / 'streamlines.svg'
    plt.savefig(outpath, dpi=300, format='svg',
                bbox_inches='tight', facecolor='none', transparent=True)
    plt.savefig(outpath.with_suffix('.pdf'), dpi=300, format='pdf',
                bbox_inches='tight', facecolor='none', transparent=True)
    plt.savefig(outpath.with_suffix('.png'), dpi=200, format='png',
                bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    plot_phase_portrait(env_params(), outdir=here / 'plots' / 'streamlines')
