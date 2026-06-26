import glob
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# ── Locations ───────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
OUT_BASE = ROOT_DIR / 'outputs' / 'survivalProb'
PLOT_DIR = ROOT_DIR / 'plots' / 'survival_prob'


def read_param_dir(param_dir):
    """Pool the death times of every sim_*.tsv in one parameter directory and
    read the run metadata from the header. Returns (alpha, t_max, deaths, n)."""
    deaths = []
    alpha = None
    t_max = None
    for fname in sorted(glob.glob(str(Path(param_dir) / 'sim_*.tsv'))):
        with open(fname) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    if line.startswith('# alpha'):
                        alpha = float(line.split('\t')[1])
                    elif line.startswith('# t_max'):
                        t_max = float(line.split('\t')[1])
                    continue
                deaths.append(float(line))
    if alpha is None:
        # fall back to parsing the directory name, e.g. alpha_1.55652_tmax_600_...
        m = re.search(r'alpha_([0-9.]+)', Path(param_dir).name)
        alpha = float(m.group(1)) if m else np.nan
    deaths = np.asarray(deaths, dtype=float)
    if t_max is None:
        finite = deaths[deaths >= 0.0]
        t_max = float(finite.max()) if finite.size else 1.0
    return alpha, t_max, deaths


def survival_curve(deaths, t_grid):
    """S(t) = P(predator still alive at time t).

    A trajectory with death_time < 0 survived the whole horizon; otherwise it is
    alive at time t iff t < death_time."""
    n = deaths.size
    alive_forever = deaths < 0.0
    # alive at t  <=>  survived forever, or died strictly later than t
    counts = np.array([
        np.count_nonzero(alive_forever | (deaths > t)) for t in t_grid
    ], dtype=float)
    return counts / n


def main():
    param_dirs = sorted(p for p in glob.glob(str(OUT_BASE / '*')) if Path(p).is_dir())
    if not param_dirs:
        raise SystemExit(f"No parameter directories found under {OUT_BASE}")

    runs = [read_param_dir(p) for p in param_dirs]
    runs.sort(key=lambda r: r[0])  # order curves by alpha

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

    # Colour the curves on a perceptual ramp from slow (extinction) to fast
    # (survival) switching.
    cmap = plt.get_cmap('plasma')
    n = len(runs)
    colors = [cmap(0.12 + 0.66 * i / max(n - 1, 1)) for i in range(n)]

    fig = plt.figure(figsize=(7, 5.2))
    fig.patch.set_facecolor('none')
    gs = GridSpec(1, 1, left=0.16, right=0.96, top=0.95, bottom=0.16)
    ax = fig.add_subplot(gs[0, 0])

    t_max = max(r[1] for r in runs)
    t_grid = np.linspace(0.0, t_max, 400)

    for (alpha, _tm, deaths), color in zip(runs, colors):
        S = survival_curve(deaths, t_grid)
        ax.plot(t_grid, S, '-', color=color, linewidth=2.4,
                label=rf'$\alpha = {alpha:.3f}$')

    ax.set_xlabel(r'Time $t$', labelpad=4)
    ax.set_ylabel(r'Survival probability $S(t)$', labelpad=4)
    # ax.set_xlim(0.0, t_max)
    # ax.set_ylim(-0.02, 1.02)
    # ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(False)
    ax.patch.set_facecolor('none')
    ax.legend(frameon=False, fontsize=18, loc='lower left')

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    stem = f'survival_prob_n{len(runs)}_tmax{t_max:g}'
    for ext in ('svg', 'pdf'):
        outpath = PLOT_DIR / f'{stem}.{ext}'
        plt.savefig(outpath, dpi=300, format=ext,
                    bbox_inches='tight', facecolor='none', transparent=True)
        print(f'Saved: {outpath}')
    plt.close(fig)


if __name__ == '__main__':
    main()
