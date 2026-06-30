r"""y* and the invasion exponent lambda versus the switching rate alpha.

Two environments A and B are alternated.  In *deterministic* (periodic)
switching the system dwells exactly 1/alpha in each environment (period
2/alpha, half in A then half in B).  In *stochastic* switching the environment
is a symmetric telegraph process, toggling A<->B at rate alpha (mean dwell
1/alpha).

For each alpha we measure two things, for both switching protocols:

  * lambda  — the invasion / Floquet exponent.  With the predator rare (y -> 0)
              the prey relaxes to its resident attractor x(t) and the predator
              log-density grows as d(log y)/dt = x - gamma.  Hence
                  lambda = < x(t) - gamma(t) >
              averaged over the resident attractor.  lambda > 0 => the predator
              can invade and persist; lambda < 0 => extinction.

  * y*      — the long-time-averaged predator density of the full coexisting
              system, <y>.  It lifts off zero exactly where lambda crosses zero.

Both protocols share the same fast-switching limit (the time-averaged field,
lambda_fast = xbar_fast - gammabar > 0) and the same slow-switching limit
(long stays at each carrying capacity, lambda_slow = (K_A+K_B)/2 - gammabar < 0).
The survival-by-switching paradox: fast switching survives, slow switching dies.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from numba import njit


# ── Environments ────────────────────────────────────────────────────────────

def env_params():
    return {
        'A': {'r': 1.0, 'K': 1.0, 'gamma': 1.1},
        'B': {'r': 4.0, 'K': 2.0, 'gamma': 2.1},
    }


# ── Numba kernels ───────────────────────────────────────────────────────────

@njit(cache=True, inline='always')
def _logistic_rk4(x, r, K, dt):
    k1 = r * x * (1.0 - x / K)
    xa = x + 0.5 * dt * k1; k2 = r * xa * (1.0 - xa / K)
    xb = x + 0.5 * dt * k2; k3 = r * xb * (1.0 - xb / K)
    xc = x + dt * k3;       k4 = r * xc * (1.0 - xc / K)
    return x + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


@njit(cache=True, inline='always')
def _full_rk4(x, logy, r, K, g, dt):
    # derivative helper inlined (state = x, log y)
    y = np.exp(logy)
    k1x = r * x * (1.0 - x / K) - x * y;             k1l = x - g
    xa = x + 0.5 * dt * k1x; la = logy + 0.5 * dt * k1l; ya = np.exp(la)
    k2x = r * xa * (1.0 - xa / K) - xa * ya;          k2l = xa - g
    xb = x + 0.5 * dt * k2x; lb = logy + 0.5 * dt * k2l; yb = np.exp(lb)
    k3x = r * xb * (1.0 - xb / K) - xb * yb;          k3l = xb - g
    xc = x + dt * k3x;       lc = logy + dt * k3l;     yc = np.exp(lc)
    k4x = r * xc * (1.0 - xc / K) - xc * yc;          k4l = xc - g
    x2 = x + dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
    l2 = logy + dt / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l)
    if x2 < 1e-14:
        x2 = 1e-14
    if l2 < -700.0:
        l2 = -700.0
    elif l2 > 700.0:
        l2 = 700.0
    return x2, l2


# ---- deterministic (periodic) switching -------------------------------------

@njit(cache=True)
def lambda_periodic(alpha, rA, KA, gA, rB, KB, gB,
                    n_trans, n_meas, dt_max):
    P = 2.0 / alpha
    steps = int(np.ceil(P / dt_max))
    if steps % 2 == 1:
        steps += 1
    if steps < 2:
        steps = 2
    dt = P / steps
    x = 0.5 * (KA + KB)
    for _ in range(n_trans):
        for s in range(steps):
            if 2 * s < steps:
                x = _logistic_rk4(x, rA, KA, dt)
            else:
                x = _logistic_rk4(x, rB, KB, dt)
    acc = 0.0
    for _ in range(n_meas):
        for s in range(steps):
            if 2 * s < steps:
                acc += (x - gA) * dt
                x = _logistic_rk4(x, rA, KA, dt)
            else:
                acc += (x - gB) * dt
                x = _logistic_rk4(x, rB, KB, dt)
    return acc / (n_meas * P)


@njit(cache=True)
def ystar_periodic(alpha, rA, KA, gA, rB, KB, gB,
                   x0, y0, n_trans, n_meas, dt_max):
    P = 2.0 / alpha
    steps = int(np.ceil(P / dt_max))
    if steps % 2 == 1:
        steps += 1
    if steps < 2:
        steps = 2
    dt = P / steps
    x = x0
    logy = np.log(y0)
    for _ in range(n_trans):
        for s in range(steps):
            if 2 * s < steps:
                x, logy = _full_rk4(x, logy, rA, KA, gA, dt)
            else:
                x, logy = _full_rk4(x, logy, rB, KB, gB, dt)
    acc = 0.0
    for _ in range(n_meas):
        for s in range(steps):
            acc += np.exp(logy) * dt
            if 2 * s < steps:
                x, logy = _full_rk4(x, logy, rA, KA, gA, dt)
            else:
                x, logy = _full_rk4(x, logy, rB, KB, gB, dt)
    return acc / (n_meas * P)


# ---- stochastic (telegraph) switching ---------------------------------------

@njit(cache=True)
def lambda_stochastic(alpha, rA, KA, gA, rB, KB, gB,
                      T_trans, T_meas, dt, n_runs, seed):
    np.random.seed(seed)
    pswitch = alpha * dt
    nt = int(T_trans / dt)
    nm = int(T_meas / dt)
    total = 0.0
    for _ in range(n_runs):
        x = 0.5 * (KA + KB)
        env = 0 if np.random.random() < 0.5 else 1
        for _ in range(nt):
            if env == 0:
                x = _logistic_rk4(x, rA, KA, dt)
            else:
                x = _logistic_rk4(x, rB, KB, dt)
            if np.random.random() < pswitch:
                env = 1 - env
        acc = 0.0
        for _ in range(nm):
            if env == 0:
                acc += (x - gA) * dt
                x = _logistic_rk4(x, rA, KA, dt)
            else:
                acc += (x - gB) * dt
                x = _logistic_rk4(x, rB, KB, dt)
            if np.random.random() < pswitch:
                env = 1 - env
        total += acc / (nm * dt)
    return total / n_runs


@njit(cache=True)
def ystar_stochastic(alpha, rA, KA, gA, rB, KB, gB,
                     x0, y0, T_trans, T_meas, dt, n_runs, seed):
    np.random.seed(seed)
    pswitch = alpha * dt
    nt = int(T_trans / dt)
    nm = int(T_meas / dt)
    total = 0.0
    for _ in range(n_runs):
        x = x0
        logy = np.log(y0)
        env = 0 if np.random.random() < 0.5 else 1
        for _ in range(nt):
            if env == 0:
                x, logy = _full_rk4(x, logy, rA, KA, gA, dt)
            else:
                x, logy = _full_rk4(x, logy, rB, KB, gB, dt)
            if np.random.random() < pswitch:
                env = 1 - env
        acc = 0.0
        for _ in range(nm):
            acc += np.exp(logy) * dt
            if env == 0:
                x, logy = _full_rk4(x, logy, rA, KA, gA, dt)
            else:
                x, logy = _full_rk4(x, logy, rB, KB, gB, dt)
            if np.random.random() < pswitch:
                env = 1 - env
        total += acc / (nm * dt)
    return total / n_runs


# ── Theory: invasion exponent lambda(alpha) ──────────────────────────────────
# Closed-form logistic flow on the y = 0 boundary (overleaf eq. for X_E).

@njit(cache=True, inline='always')
def _Xflow(t, xin, r, K):
    e = np.exp(r * t)
    return K * xin * e / (K + xin * (e - 1.0))


# -- deterministic (periodic) Floquet exponent --------------------------------
# Period T = 2/alpha; each environment is held for T/2 = 1/alpha.  The boundary
# orbit turning points (x_a, x_b) solve the periodicity fixed point, and
#   lambda = (K1+K2-g1-g2)/2 - (alpha/2)(K1/r1 - K2/r2) ln(x_a/x_b).

@njit(cache=True)
def _det_orbit(alpha, rA, KA, rB, KB, iters):
    half = 1.0 / alpha
    xb = 0.5 * (KA + KB)
    for _ in range(iters):
        xa = _Xflow(half, xb, rA, KA)
        xb = _Xflow(half, xa, rB, KB)
    xa = _Xflow(half, xb, rA, KA)
    return xa, xb


@njit(cache=True)
def lambda_det_theory(alpha, rA, KA, gA, rB, KB, gB):
    xa, xb = _det_orbit(alpha, rA, KA, rB, KB, 600)
    return (0.5 * (KA + KB - gA - gB)
            - 0.5 * alpha * (KA / rA - KB / rB) * np.log(xa / xb))


# -- stochastic (Poisson) invasion exponent -----------------------------------
# Same form, with ln(x_a/x_b) replaced by its stationary average A(alpha) over
# the random boundary map (exponential dwell times).  No closed form, so A is
# estimated by iterating the closed-form flow with random dwells.

@njit(cache=True)
def _A_stoch(alpha, rA, KA, rB, KB, n_burn, n_samp, seed):
    np.random.seed(seed)
    xb = 0.5 * (KA + KB)
    for _ in range(n_burn):
        t1 = -np.log(np.random.random()) / alpha
        xa = _Xflow(t1, xb, rA, KA)
        t2 = -np.log(np.random.random()) / alpha
        xb = _Xflow(t2, xa, rB, KB)
    s = 0.0
    for _ in range(n_samp):
        t1 = -np.log(np.random.random()) / alpha
        xa = _Xflow(t1, xb, rA, KA)
        s += np.log(xa / xb)
        t2 = -np.log(np.random.random()) / alpha
        xb = _Xflow(t2, xa, rB, KB)
    return s / n_samp


@njit(cache=True)
def lambda_stoch_theory(alpha, rA, KA, gA, rB, KB, gB, n_burn, n_samp, seed):
    A = _A_stoch(alpha, rA, KA, rB, KB, n_burn, n_samp, seed)
    return (0.5 * (KA + KB - gA - gB)
            - 0.5 * alpha * (KA / rA - KB / rB) * A)


def _bisect_log(f, lo, hi, n=70):
    """Find the (single) sign change of an increasing f on [lo, hi] in log-x."""
    for _ in range(n):
        mid = np.sqrt(lo * hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return np.sqrt(lo * hi)


def critical_rates(params):
    A, B = params['A'], params['B']
    p = (A['r'], A['K'], A['gamma'], B['r'], B['K'], B['gamma'])
    ac_det = _bisect_log(lambda a: lambda_det_theory(a, *p), 0.05, 10.0)
    ac_sto = _bisect_log(
        lambda a: lambda_stoch_theory(a, *p, 5000, 1_000_000, 271828),
        0.05, 10.0)
    return ac_det, ac_sto


# ── Sweep ───────────────────────────────────────────────────────────────────

def compute_curves(params, alphas_dense, alphas_sim,
                   n_runs_lam, T_meas_lam, n_runs_y, T_meas_y):
    A, B = params['A'], params['B']
    p = (A['r'], A['K'], A['gamma'], B['r'], B['K'], B['gamma'])

    # -- dense theory lines (panel B) -----------------------------------------
    print("Theory curves (dense) ...")
    lam_det_th = np.array([lambda_det_theory(a, *p) for a in alphas_dense])
    # one fixed seed across alpha keeps the boundary-map estimate smooth in alpha
    lam_sto_th = np.array([lambda_stoch_theory(a, *p, 3000, 300_000, 99)
                           for a in alphas_dense])

    # -- dense deterministic y* via ODE limit cycle (panel A line) ------------
    y_det = np.empty_like(alphas_dense)
    for i, a in enumerate(alphas_dense):
        P = 2.0 / a
        n_tr = max(200, int(np.ceil(400.0 / P)))
        n_me = max(60,  int(np.ceil(60.0 / P)))
        y_det[i] = ystar_periodic(a, *p, 1.0, 0.05, n_tr, n_me, 0.01)
    y_det = np.maximum(y_det, 0.0)

    # -- dense stochastic y* via simulation (panel A line) --------------------
    print(f"Stochastic y* simulations (dense, n_runs={n_runs_y}) ...")
    y_sto = np.empty_like(alphas_dense)
    for i, a in enumerate(alphas_dense):
        y_sto[i] = ystar_stochastic(a, *p, 1.0, 0.05, 300.0, T_meas_y,
                                    0.0025, n_runs_y, 5000 + i)
    y_sto = np.maximum(y_sto, 0.0)

    # -- sparse simulation dots for lambda (panel B) -------------------------
    print(f"Lambda simulations (sparse, n_runs={n_runs_lam}) ...")
    lam_det_sim = np.empty_like(alphas_sim)
    lam_sto_sim = np.empty_like(alphas_sim)
    for i, a in enumerate(alphas_sim):
        P = 2.0 / a
        n_tr = max(200, int(np.ceil(400.0 / P)))
        n_me = max(60,  int(np.ceil(60.0 / P)))
        lam_det_sim[i] = lambda_periodic(a, *p, n_tr, n_me, 0.01)
        lam_sto_sim[i] = lambda_stochastic(a, *p, 300.0, T_meas_lam, 0.0025,
                                           n_runs_lam, 1000 + i)

    return dict(lam_det_th=lam_det_th, lam_sto_th=lam_sto_th,
                y_det=y_det, y_sto=y_sto,
                lam_det_sim=lam_det_sim, lam_sto_sim=lam_sto_sim)


# ── Plot ────────────────────────────────────────────────────────────────────

def plot(params=None, outdir=None, n_dense=100, n_sim=22,
         alpha_min=0.05, alpha_max=100.0,
         n_runs_lam=32, T_meas_lam=2000.0, n_runs_y=48, T_meas_y=2000.0):
    if params is None:
        params = env_params()

    plt.rcParams.update({
        'font.size':        24,
        'axes.titlesize':   24,
        'axes.labelsize':   26,
        'xtick.labelsize':  21,
        'ytick.labelsize':  21,
        'axes.linewidth':   1.2,
        'xtick.direction':  'in',
        'ytick.direction':  'in',
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'legend.fontsize':  17,
        'pdf.fonttype':     42,
        'svg.fonttype':     'none',
    })

    A, B = params['A'], params['B']
    gbar  = 0.5 * (A['gamma'] + B['gamma'])
    xfast = (A['r'] + B['r']) / (A['r'] / A['K'] + B['r'] / B['K'])
    lam_fast = xfast - gbar
    lam_slow = 0.5 * (A['K'] + B['K']) - gbar
    y_fast = 0.5 * (A['r'] * (1 - gbar / A['K']) + B['r'] * (1 - gbar / B['K']))

    alphas_dense = np.logspace(np.log10(alpha_min), np.log10(alpha_max), n_dense)
    alphas_sim   = np.logspace(np.log10(alpha_min), np.log10(alpha_max), n_sim)
    d = compute_curves(params, alphas_dense, alphas_sim,
                       n_runs_lam, T_meas_lam, n_runs_y, T_meas_y)

    ac_det, ac_sto = critical_rates(params)
    print(f"Critical rates (theory):  alpha_c_det = {ac_det:.3f}, "
          f"alpha_c_stoch = {ac_sto:.3f}")

    C_DET = '#1f3b73'   # deterministic — deep blue
    C_STO = '#c0392b'   # stochastic    — brick red
    C_REF = '#9aa0a6'   # asymptote guides

    # ── Two wide panels, stacked and sharing the x-axis ───────────────────────
    fig = plt.figure(figsize=(10, 9))
    gs = GridSpec(2, 1, figure=fig, hspace=0.10,
                  left=0.11, right=0.97, top=0.97, bottom=0.10)
    axY = fig.add_subplot(gs[0, 0])
    axL = fig.add_subplot(gs[1, 0], sharex=axY)

    def crit_lines(ax):
        """Vertical critical-rate guides (labelled only on the lower panel)."""
        ax.axvline(ac_det, color=C_DET, lw=1.0, ls=':', zorder=0)
        ax.axvline(ac_sto, color=C_STO, lw=1.0, ls=':', zorder=0)

    # consistent whitespace above the fast-limit line in both panels: leave the
    # same fraction GAP of the axis height between that line and the top spine.
    GAP = 0.15
    def top_lim(fast, bottom):
        return (fast - GAP * bottom) / (1.0 - GAP)

    # ── y* panel — both curves from simulation ────────────────────────────────
    axY.axhline(y_fast, color=C_REF, lw=1.2, ls='--', zorder=0)
    axY.text(alpha_max, y_fast, r' fast limit', color=C_REF,
             va='bottom', ha='right', fontsize=15)
    crit_lines(axY)
    axY.plot(alphas_dense, d['y_det'], '-', color=C_DET, lw=2.6,
             label='Deterministic')
    axY.plot(alphas_dense, d['y_sto'], '-', color=C_STO, lw=2.2,
             label='Stochastic')
    axY.set_xscale('log')
    axY.set_ylabel(r'$y^\star$')
    y_bot = -0.004
    axY.set_ylim(y_bot, top_lim(y_fast, y_bot))
    axY.tick_params(labelbottom=False)          # x-axis shared with panel below
    axY.legend(frameon=False, loc='upper left')

    # ── lambda panel — theory lines + simulation dots ─────────────────────────
    axL.axhline(0.0, color='#1a1a1a', lw=1.0, zorder=0)
    for val in (lam_fast, lam_slow):
        axL.axhline(val, color=C_REF, lw=1.2, ls='--', zorder=0)
    axL.text(alpha_max, lam_fast, r' fast limit', color=C_REF,
             va='bottom', ha='right', fontsize=15)
    axL.text(alpha_min, lam_slow, r' slow limit', color=C_REF,
             va='bottom', ha='left', fontsize=15)
    crit_lines(axL)
    # critical-rate labels: only here, upright, just right of each line and at
    # the same height as the "slow limit" text.
    axL.text(ac_det, lam_slow, r' $\alpha_x^\mathrm{det}$', color=C_DET,
             va='bottom', ha='left', fontsize=15)
    axL.text(ac_sto, lam_slow, r' $\alpha_x^\mathrm{stoch}$', color=C_STO,
             va='bottom', ha='left', fontsize=15)
    axL.plot(alphas_dense, d['lam_det_th'], '-', color=C_DET, lw=2.6,
             label='Deterministic')
    axL.plot(alphas_dense, d['lam_sto_th'], '-', color=C_STO, lw=2.2,
             label='Stochastic')
    axL.plot(alphas_sim, d['lam_det_sim'], 'o', color=C_DET, ms=6,
             markerfacecolor='white', markeredgewidth=1.6)
    axL.plot(alphas_sim, d['lam_sto_sim'], 's', color=C_STO, ms=6,
             markerfacecolor='white', markeredgewidth=1.6)
    axL.set_xscale('log')
    axL.set_xlabel(r'Switching rate  $\alpha$')
    axL.set_ylabel(r'$\lambda$')
    lam_bot = lam_slow - 0.012
    axL.set_ylim(lam_bot, top_lim(lam_fast, lam_bot))
    axL.legend(frameon=False, loc='upper left')

    # ── Save ──────────────────────────────────────────────────────────────────
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tag = (f"rA{A['r']:g}_KA{A['K']:g}_gA{A['gamma']:g}"
           f"_rB{B['r']:g}_KB{B['K']:g}_gB{B['gamma']:g}")
    base = f"ystar_lambda_vs_alpha_{tag}"
    for ext, kw in (('svg', dict(facecolor='none', transparent=True)),
                    ('pdf', dict(facecolor='none', transparent=True)),
                    ('png', dict(facecolor='white'))):
        fig.savefig(outdir / f"{base}.{ext}", dpi=200, bbox_inches='tight', **kw)
    plt.close(fig)
    print(f"Saved: {outdir / base}.svg / .pdf / .png")


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    plot(env_params(), outdir=here / 'plots' / 'ystar_lambda_vs_alpha')
