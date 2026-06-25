import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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


# ── Main plot ─────────────────────────────────────────────────────────────────

def plot_phase_portrait(params=None, outdir=None):
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

    # ── Axis limits ───────────────────────────────────────────────────────────
    all_y      = [fy for (_, fy) in fps_A + fps_B + fps_avg]
    fp_y_max   = max(all_y)
    x_data_max = max(envA['K'], envB['K'],
                     max(fx for (fx, _) in fps_A + fps_B + fps_avg))
    x_lo, x_hi = 0.0, x_data_max * 1.18
    y_lo, y_hi = 0.0, max(fp_y_max * 3.2, x_data_max * 0.25)

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

    # ── Style ─────────────────────────────────────────────────────────────────
    STREAM_CMAP = 'plasma'   # dark-purple → magenta → bright orange; no harsh yellow
    C_FP_STABLE   = '#1a1a1a'   # filled: near-black
    C_FP_UNSTABLE = 'white'     # open circle

    # ── Layout: three panels, no colorbar ────────────────────────────────────
    fig = plt.figure(figsize=(14, 4.6))
    fig.patch.set_facecolor('none')
    gs = GridSpec(1, 3, wspace=0.42, left=0.07, right=0.98,
                  top=0.91, bottom=0.16)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

    panels = [
        dict(ax=axes[0], U=UA,  V=VA,  sp=spA,  fps=fps_A,
             title='Environment A',
             jac_fn=lambda x, y: jacobian(envA, x, y)),
        dict(ax=axes[1], U=UB,  V=VB,  sp=spB,  fps=fps_B,
             title='Environment B',
             jac_fn=lambda x, y: jacobian(envB, x, y)),
        dict(ax=axes[2], U=Uav, V=Vav, sp=spAv, fps=fps_avg,
             title='Time-averaged field',
             jac_fn=lambda x, y: averaged_jacobian(envA, envB, x, y)),
    ]

    for p in panels:
        ax = p['ax']

        # streamlines
        ax.streamplot(
            xg, yg, p['U'], p['V'],
            color=p['sp'],
            cmap=STREAM_CMAP,
            norm=norm,
            linewidth=1.1,
            density=1.1,
            arrowsize=1.4,
            arrowstyle='->',
            minlength=0.05,
            zorder=1,
        )

        # fixed points
        for (fx, fy) in p['fps']:
            stab = is_stable(p['jac_fn'](fx, fy))
            # fc = C_FP_STABLE if stab else C_FP_UNSTABLE
            # ax.plot(fx, fy, 'o',
            #         color=fc,
            #         markeredgecolor=C_FP_STABLE,
            #         markeredgewidth=1.6,
            #         markersize=12,
            #         zorder=5)
            if stab:
                ax.plot(fx, fy, 'o', 
                        color=C_FP_STABLE,
                        markeredgecolor=C_FP_STABLE,
                        markeredgewidth=1.6,
                        markersize=12,
                        zorder=5)

        ax.set_title(p['title'], pad=7)
        ax.set_xlabel('Prey  $x$', labelpad=4)
        ax.set_ylabel('Predator  $y$', labelpad=4)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.grid(False)
        ax.patch.set_facecolor('none')

    # ── Save ──────────────────────────────────────────────────────────────────
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / 'streamlines.svg'
    plt.savefig(outpath, dpi=300, format='svg',
                bbox_inches='tight', facecolor='none', transparent=True)
    plt.savefig(outpath.with_suffix('.pdf'), dpi=300, format='pdf',
                bbox_inches='tight', facecolor='none', transparent=True)
    plt.close(fig)
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    plot_phase_portrait(env_params(), outdir='src/3paramLogisticLV/plots/vector_fields')