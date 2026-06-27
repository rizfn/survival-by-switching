import glob
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
OUT_BASE = ROOT_DIR / 'outputs' / 'extinctionScaling'
PLOT_DIR = ROOT_DIR / 'plots' / 'extinction_scaling'

ALPHA_C = 1.5565


def read_scaling_dir(param_dir):
    """Return (alpha, Delta, tau, n_died, n_total) arrays from one alpha dir."""
    alpha = None
    rows = []
    with open(Path(param_dir) / 'scaling.tsv') as handle:
        for line in handle:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                if s.startswith('# alpha'):
                    alpha = float(s.split('\t')[1])
                continue
            if s.startswith('y_ext'):
                continue
            y_ext, Delta, n_died, n_total, tau = s.split('\t')
            rows.append((float(Delta), float(tau), float(n_died), float(n_total)))
    rows.sort()
    arr = np.array(rows)
    return alpha, arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def regime(alpha):
    if alpha < ALPHA_C - 1e-3:
        return 'sub-critical'
    if alpha > ALPHA_C + 1e-3:
        return 'super-critical'
    return 'critical'


# Three representative rates (sub / critical / super) for this teaching figure;
# the full alpha family lives in plot_critical_point.py.
REPRESENTATIVE = (1.24521, 1.55652, 1.86782)


def main():
    dirs = sorted(p for p in glob.glob(str(OUT_BASE / 'alpha_*')) if Path(p).is_dir())
    if not dirs:
        raise SystemExit(f"No scaling directories under {OUT_BASE}")
    runs = [read_scaling_dir(d) for d in dirs]
    runs = [r for r in runs if any(abs(r[0] - a) < 1e-3 for a in REPRESENTATIVE)]
    runs.sort(key=lambda r: r[0])

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

    cmap = plt.get_cmap('plasma')
    n = len(runs)
    colors = [cmap(0.12 + 0.66 * i / max(n - 1, 1)) for i in range(n)]

    fig = plt.figure(figsize=(7.2, 5.6))
    fig.patch.set_facecolor('none')
    gs = GridSpec(1, 1, left=0.15, right=0.96, top=0.96, bottom=0.14)
    ax = fig.add_subplot(gs[0, 0])

    for (alpha, Delta, tau, n_died, n_total), color in zip(runs, colors):
        # only trust barriers with enough recorded extinctions
        ok = n_died >= 50
        D, T = Delta[ok], tau[ok]
        reg = regime(alpha)
        ax.plot(D, T, 'o', color=color, markersize=7, zorder=3,
                label=rf'$\alpha={alpha:.3f}$  ({reg})')

        # Overlay the asymptotic (large-barrier) law. The robust contrast is
        # polynomial growth below alpha_c versus exponential growth above it;
        # exactly at alpha_c the persistence time is the marginal separatrix
        # (mean-based tau is heavy-tailed there, so we do not force a power law).
        tail = D >= 6.0
        Dt, Tt = D[tail], T[tail]
        Dfit = np.linspace(Dt.min(), Dt.max(), 100)
        if reg == 'sub-critical':
            a, b = np.polyfit(Dt, Tt, 1)              # tau ~ Delta (ballistic)
            ax.plot(Dfit, a * Dfit + b, '--', color=color, lw=1.6, zorder=2)
            txt = r'$\tau\propto\Delta$  (polynomial)'
        elif reg == 'super-critical':
            c, d = np.polyfit(Dt, np.log(Tt), 1)      # tau ~ exp(c Delta) (Kramers)
            ax.plot(Dfit, np.exp(c * Dfit + d), '--', color=color, lw=1.6, zorder=2)
            txt = rf'$\tau\propto e^{{{c:.2f}\,\Delta}}$  (exponential)'
        else:
            txt = 'marginal'
        ax.annotate(txt, xy=(D[-1], T[-1]), xytext=(6, 0),
                    textcoords='offset points', fontsize=15, color=color,
                    va='center')

    ax.set_yscale('log')
    ax.set_xlabel(r'Barrier height  $\Delta=\log(y_0/y_{\rm ext})\;\sim\;\log N$', labelpad=4)
    ax.set_ylabel(r'Mean extinction time  $\tau$', labelpad=4)
    ax.set_xlim(left=2.3)
    ax.patch.set_facecolor('none')
    ax.legend(frameon=False, fontsize=15, loc='upper left')

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    stem = 'extinction_time_scaling'
    for ext in ('svg', 'pdf'):
        outpath = PLOT_DIR / f'{stem}.{ext}'
        plt.savefig(outpath, dpi=300, format=ext,
                    bbox_inches='tight', facecolor='none', transparent=True)
        print(f'Saved: {outpath}')
    plt.close(fig)


if __name__ == '__main__':
    main()
