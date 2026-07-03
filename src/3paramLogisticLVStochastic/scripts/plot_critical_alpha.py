r"""Critical switching rate alpha_c versus the swept predator death rate gamma_2.

For each gamma_2 we locate the critical switching rate alpha_c at which the
predator's invasion exponent vanishes, lambda(alpha_c) = 0: below alpha_c the
predator dies out (switching too slow), above it the predator survives. We do
this for both switching protocols and by two routes each, mirroring exactly the
methods in ystar_lambda_vs_alpha.py:

                       lambda(alpha) obtained from ...
  --------------------------------------------------------------------------
  deterministic  theory : closed-form Floquet exponent; the boundary orbit
                          (x_a, x_b) is found by fixed-point iteration of the
                          exact logistic flow X_E, then
                          lambda = (K1+K2-g1-g2)/2 - (alpha/2)(K1/r1-K2/r2)ln(xa/xb).
  deterministic  sim.   : direct RK4 of the periodic-switching resident orbit,
                          time-averaging x(t) - gamma(t).
  stochastic     theory : same closed form with ln(xa/xb) replaced by its
                          stationary average A(alpha) = E[ln(xa/xb)] over
                          exponential dwells, iterating the exact flow X_E.
  stochastic     sim.   : direct RK4 of the telegraph-switching resident orbit,
                          time-averaging x(t) - gamma(t).

In every one of the four cases alpha_c is the root of the corresponding
lambda(alpha) = 0, found by the *same* log-domain bisection. Theory is drawn as
dense lines, simulation as sparse markers; the markers landing on the lines is
the validation.

Environments:  A (slow): r=1, K=1, gamma=1.1 ;  B (fast): r=4, K=2, gamma swept.
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
        'B': {'r': 4.0, 'K': 2.0},          # gamma_B is swept
    }


# ── Numba kernels (identical to ystar_lambda_vs_alpha.py) ────────────────────

@njit(cache=True, inline='always')
def _logistic_rk4(x, r, K, dt):
    """One RK4 step of the resident logistic flow dx/dt = r x (1 - x/K)."""
    k1 = r * x * (1.0 - x / K)
    xa = x + 0.5 * dt * k1; k2 = r * xa * (1.0 - xa / K)
    xb = x + 0.5 * dt * k2; k3 = r * xb * (1.0 - xb / K)
    xc = x + dt * k3;       k4 = r * xc * (1.0 - xc / K)
    return x + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

@njit(cache=True, inline='always')
def _Xflow(t, xin, r, K):
    # numerically stable logistic flow: avoid exp(r*t) overflow
    z = np.exp(-r * t)
    return K / (1.0 + (K / xin - 1.0) * z)


# -- deterministic (periodic) SIMULATION: lambda by direct RK4 ----------------
# Period P = 2/alpha, each environment held P/2 = 1/alpha. We relax onto the
# periodic resident orbit (n_trans periods) then time-average x(t) - gamma(t).

@njit(cache=True)
def lambda_periodic_sim(alpha, rA, KA, gA, rB, KB, gB, n_trans, n_meas, dt_max):
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


# -- stochastic (telegraph) SIMULATION: lambda by direct RK4 ------------------
# Symmetric two-state telegraph: toggle A<->B with probability alpha*dt per step
# (mean dwell 1/alpha). Average x(t) - gamma(t) over the resident orbit.

@njit(cache=True)
def lambda_stochastic_sim(alpha, rA, KA, gA, rB, KB, gB,
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


# -- deterministic THEORY: closed-form Floquet exponent -----------------------

@njit(cache=True)
def _det_orbit(alpha, rA, KA, rB, KB, iters):
    """Turning points (x_a, x_b) of the equal-duty boundary orbit, by
    fixed-point iteration of the two-stage closed-form map."""
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


# -- stochastic THEORY: boundary-map average A(alpha) -------------------------
# Same closed form, with ln(xa/xb) replaced by its stationary mean over the
# random boundary map (exponential dwells). One fixed seed across alpha keeps
# the estimate smooth in alpha, so the bisection below sees a clean sign change.

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


# ── Shared alpha_c root finder ───────────────────────────────────────────────
# lambda(alpha) is negative at slow switching (extinction) and positive at fast
# switching (survival), crossing zero once. We bracket that single crossing and
# bisect in log(alpha) -- the same routine for all four lambda(alpha) above, so
# deterministic and stochastic alpha_c are computed by an identical algorithm.

def alpha_c_root(f, lo=0.05, hi=10.0, n_bisect=40):
    """Root of f(alpha)=0 (f increasing through one zero) by log bisection.
    Returns nan if no sign change can be bracketed (e.g. gamma_2 outside the
    switching-rescue window, where the predator dies at every alpha)."""
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0.0:
        # Widen the fast side first (alpha_c can be large near the upper edge),
        # then the slow side, before giving up.
        bracketed = False
        for hi_try in (hi * 5.0, hi * 25.0, hi * 125.0):
            fhi_try = f(hi_try)
            if flo * fhi_try < 0.0:
                hi, fhi, bracketed = hi_try, fhi_try, True
                break
        if not bracketed:
            for lo_try in (lo / 5.0, lo / 25.0, lo / 125.0):
                flo_try = f(lo_try)
                if flo_try * fhi < 0.0:
                    lo, flo, bracketed = lo_try, flo_try, True
                    break
        if not bracketed:
            return np.nan
    for _ in range(n_bisect):
        mid = np.sqrt(lo * hi)                 # geometric (log) midpoint
        fmid = f(mid)
        if flo * fmid <= 0.0:
            hi = mid
        else:
            lo, flo = mid, fmid
    return np.sqrt(lo * hi)


def alpha_c_sim_mean(lam_of_alpha_seed, n_reps, seed0, **root_kw):
    """Simulation alpha_c averaged over independent seeds (mean, std). Each rep
    fixes a seed so lambda(alpha) is deterministic and the bisection is clean;
    different reps give independent estimates."""
    vals = []
    for rep in range(n_reps):
        ac = alpha_c_root(lambda a: lam_of_alpha_seed(a, seed0 + rep), **root_kw)
        if np.isfinite(ac):
            vals.append(ac)
    if not vals:
        return np.nan, np.nan
    return float(np.mean(vals)), float(np.std(vals))


# ── Sweep over gamma_2 ───────────────────────────────────────────────────────

def compute(params, gB_curve, gB_points,
            # deterministic
            det_dt=0.01, det_nbisect=40,
            # stochastic theory (boundary map)
            sto_burn=3000, sto_samp=150_000, sto_theory_seed=99, sto_th_nbisect=40,
            # stochastic simulation (telegraph ODE)
            sto_dt=0.0025, sto_Ttrans=300.0, sto_Tmeas=1200.0,
            sto_nruns=16, sto_nreps=3, sto_sim_seed0=5000, sto_sim_nbisect=24):

    A, B = params['A'], params['B']
    rA, KA, gA = A['r'], A['K'], A['gamma']
    rB, KB = B['r'], B['K']

    def det_counts(alpha):
        """Absolute-time-matched period counts (mirrors ystar file)."""
        P = 2.0 / alpha
        n_tr = max(200, int(np.ceil(400.0 / P)))
        n_me = max(60,  int(np.ceil(60.0 / P)))
        return n_tr, n_me

    # ── Warm up the JIT once, outside the timed loops ─────────────────────────
    lambda_det_theory(1.0, rA, KA, gA, rB, KB, 2.1)
    lambda_stoch_theory(1.0, rA, KA, gA, rB, KB, 2.1, 10, 10, 0)
    lambda_periodic_sim(1.0, rA, KA, gA, rB, KB, 2.1, 2, 2, 0.05)
    lambda_stochastic_sim(1.0, rA, KA, gA, rB, KB, 2.1, 1.0, 1.0, 0.01, 1, 0)

    # ── Theory lines (dense) ──────────────────────────────────────────────────
    print("Deterministic theory (closed form) ...")
    ac_det_th = np.array([
        alpha_c_root(lambda a: lambda_det_theory(a, rA, KA, gA, rB, KB, gB),
                     n_bisect=det_nbisect)
        for gB in gB_curve])

    print("Stochastic theory (boundary-map A) ...")
    ac_sto_th = np.array([
        alpha_c_root(lambda a: lambda_stoch_theory(a, rA, KA, gA, rB, KB, gB,
                                                   sto_burn, sto_samp, sto_theory_seed),
                     n_bisect=sto_th_nbisect)
        for gB in gB_curve])

    # ── Simulation points (sparse) ────────────────────────────────────────────
    print("Deterministic simulation (periodic RK4) ...")
    ac_det_sim = []
    for gB in gB_points:
        def f_det(a, _seed, gB=gB):
            n_tr, n_me = det_counts(a)
            return lambda_periodic_sim(a, rA, KA, gA, rB, KB, gB, n_tr, n_me, det_dt)
        m, _ = alpha_c_sim_mean(f_det, 1, 0, n_bisect=det_nbisect)   # deterministic
        ac_det_sim.append(m)
        print(f"  gamma_2={gB:.3f}: alpha_c(det sim) = {m}")
    ac_det_sim = np.array(ac_det_sim)

    print(f"Stochastic simulation (telegraph RK4, n_runs={sto_nruns}, "
          f"n_reps={sto_nreps}) ...")
    ac_sto_sim, ac_sto_sim_std = [], []
    for gB in gB_points:
        def f_sto(a, seed, gB=gB):
            return lambda_stochastic_sim(a, rA, KA, gA, rB, KB, gB,
                                         sto_Ttrans, sto_Tmeas, sto_dt, sto_nruns, seed)
        m, s = alpha_c_sim_mean(f_sto, sto_nreps, sto_sim_seed0,
                                n_bisect=sto_sim_nbisect)
        ac_sto_sim.append(m)
        ac_sto_sim_std.append(s)
        print(f"  gamma_2={gB:.3f}: alpha_c(stoch sim) = {m} +/- {s}")
    ac_sto_sim = np.array(ac_sto_sim)
    ac_sto_sim_std = np.array(ac_sto_sim_std)

    return dict(gB_curve=gB_curve, gB_points=gB_points,
                ac_det_th=ac_det_th, ac_sto_th=ac_sto_th,
                ac_det_sim=ac_det_sim,
                ac_sto_sim=ac_sto_sim, ac_sto_sim_std=ac_sto_sim_std)


# ── Plot ─────────────────────────────────────────────────────────────────────

def plot(params=None, outdir=None,
         gB_lo=2.00, gB_hi=2.24, n_dense=60, n_points=11,
         show_stoch_errorbars=False, **compute_kw):
    if params is None:
        params = env_params()

    plt.rcParams.update({
        'font.size':        22,
        'axes.titlesize':   22,
        'axes.labelsize':   26,
        'xtick.labelsize':  20,
        'ytick.labelsize':  20,
        'axes.linewidth':   1.2,
        'xtick.direction':  'in',
        'ytick.direction':  'in',
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'legend.fontsize':  16,
        'pdf.fonttype':     42,
        'svg.fonttype':     'none',
    })

    gB_curve  = np.linspace(gB_lo, gB_hi, n_dense)
    # keep the simulation points inset from the window edges, where alpha_c
    # diverges and the sign change is hard to bracket from a simulation.
    gB_points = np.linspace(gB_lo + 0.02, gB_hi - 0.03, n_points)

    d = compute(params, gB_curve, gB_points, **compute_kw)

    # sanity print at the standard operating point
    A, B = params['A'], params['B']
    p = (A['r'], A['K'], A['gamma'], B['r'], B['K'])
    ac_det_21 = alpha_c_root(lambda a: lambda_det_theory(a, *p, 2.1))
    ac_sto_21 = alpha_c_root(lambda a: lambda_stoch_theory(a, *p, 2.1, 5000, 10**6, 271828))
    print(f"\nAt gamma_2 = 2.1:  alpha_c(det) = {ac_det_21:.3f}, "
          f"alpha_c(stoch) = {ac_sto_21:.3f}")

    C_DET = '#1f3b73'   # deterministic — deep blue
    C_STO = '#c0392b'   # stochastic    — brick red

    fig = plt.figure(figsize=(7.4, 5.6))
    fig.patch.set_facecolor('none')
    gs = GridSpec(1, 1, left=0.17, right=0.96, top=0.95, bottom=0.16)
    ax = fig.add_subplot(gs[0, 0])
    ax.patch.set_facecolor('none')

    # theory lines
    ax.plot(d['gB_curve'], d['ac_det_th'], '-', color=C_DET, lw=2.6, zorder=3,
            label='Deterministic — theory')
    ax.plot(d['gB_curve'], d['ac_sto_th'], '-', color=C_STO, lw=2.6, zorder=3,
            label='Stochastic — theory')

    # simulation markers (open, matching the attached file's style)
    ax.plot(d['gB_points'], d['ac_det_sim'], 'o', color=C_DET, ms=8,
            markerfacecolor='white', markeredgewidth=1.8, zorder=5,
            label='Deterministic — simulation')
    if show_stoch_errorbars:
        ax.errorbar(d['gB_points'], d['ac_sto_sim'], yerr=d['ac_sto_sim_std'],
                    fmt='s', color=C_STO, ms=8, markerfacecolor='white',
                    markeredgewidth=1.8, capsize=4, elinewidth=1.4, zorder=5,
                    label='Stochastic — simulation')
    else:
        ax.plot(d['gB_points'], d['ac_sto_sim'], 's', color=C_STO, ms=8,
                markerfacecolor='white', markeredgewidth=1.8, zorder=5,
                label='Stochastic — simulation')

    ax.set_xlabel(r'$\gamma_2$', labelpad=4)
    ax.set_ylabel(r'Critical rate  $\alpha_c$', labelpad=4)
    ax.set_yscale('log')
    ax.grid(False)
    ax.legend(frameon=False, loc='upper left')

    # ── Save (filename carries the fixed parameters + the swept range) ────────
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tag = (f"rA{A['r']:g}_KA{A['K']:g}_gA{A['gamma']:g}"
           f"_rB{B['r']:g}_KB{B['K']:g}"
           f"_gBsweep{gB_lo:g}-{gB_hi:g}")
    outpath = outdir / f'critical_alpha_{tag}.svg'
    fig.savefig(outpath, dpi=300, bbox_inches='tight',
                facecolor='none', transparent=True)
    fig.savefig(outpath.with_suffix('.pdf'), dpi=300, bbox_inches='tight',
                facecolor='none', transparent=True)
    plt.close(fig)
    print(f"Saved: {outpath} (+ .pdf)")


if __name__ == '__main__':
    plot(env_params(), outdir='src/3paramLogisticLVStochastic/plots/critical_alpha')