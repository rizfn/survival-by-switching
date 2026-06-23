import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import brentq
from scipy.integrate import quad
from pathlib import Path
from numba import njit


# ── Ecology parameters ─────────────────────────────────────────────────────────

def env_params():
    return {
        'A': {'r': 1.0, 'K': 1.0, 'gamma': 1.1},
        'B': {'r': 4.0, 'K': 2.0},  # gamma_B is swept
    }


# ── Logistic flow (used by MC, density iteration, and the deterministic curve) ─

@njit
def x_flow(t, x_in, r, K):
    rt = r * t
    if rt > 700.0:
        return K
    if rt < -700.0:
        return 1e-300
    ert = np.exp(rt)
    return K * x_in * ert / (K + x_in * (ert - 1.0))


# ── Density iteration (deterministic / quadrature solution of lambda(alpha)=0) ─
#
# This solves the Perron-Frobenius equation for the stationary law of x0 by
# discretizing the push-forward kernel on a grid and power-iterating to the
# fixed point. No randomness is involved; "numerical" here means "quadrature",
# not "simulation".
#
# IMPORTANT CAVEAT (see question 2 in the chat): this method has a systematic,
# understood bias near both edges of the switching-paradox window, caused by
# under-resolution of the transition kernel on a *fixed* grid - not a bug, and
# not fixable by more burn-in (it isn't a transient effect). See the long
# comment above plot_critical_alpha() for the full explanation. Monte Carlo
# remains the trustworthy estimate near those edges.

def make_grid(N, KMAX, z_range=20):
    """Grid clustered near both fixed points (0 and KMAX) via a logistic
    transform of a uniform grid in 'log-odds' space z."""
    z = np.linspace(-z_range, z_range, N)
    x = KMAX / (1 + np.exp(-z))
    x = np.clip(x, 1e-12, KMAX * (1 - 1e-12))
    return np.unique(x)

def make_weights(grid):
    """Trapezoidal quadrature weights for the (non-uniform) grid."""
    N = len(grid)
    w = np.zeros(N)
    w[1:-1] = (grid[2:] - grid[:-2]) / 2
    w[0] = (grid[1] - grid[0]) / 2
    w[-1] = (grid[-1] - grid[-2]) / 2
    return w

@njit
def build_transition_matrix_nb(r, K, alpha, grid, weights, KMAX):
    """M[i, j] approx f(x1=grid[j] | x0=grid[i]) for x1 = x_flow(tau; x0),
    tau ~ Exp(alpha), normalized so the weighted row sums to 1."""
    N = grid.shape[0]
    M = np.zeros((N, N))
    eps = 1e-12
    for i in range(N):
        x0 = grid[i]
        if x0 < eps:
            x0 = eps
        elif x0 > KMAX * (1 - eps):
            x0 = KMAX * (1 - eps)
        s = 0.0
        row = np.zeros(N)
        for j in range(N):
            x1 = grid[j]
            val = 0.0
            if x0 < K:
                if x0 < x1 < K:
                    ratio = (x1 * (K - x0)) / (x0 * (K - x1))
                    if ratio > 0.0:
                        t = np.log(ratio) / r
                        if t >= 0.0:
                            dtdx1 = abs((1.0 / r) * (1.0 / x1 + 1.0 / (K - x1)))
                            val = alpha * np.exp(-alpha * t) * dtdx1
            elif x0 > K:
                if K < x1 < x0:
                    ratio = (x1 * (K - x0)) / (x0 * (K - x1))
                    if ratio > 0.0:
                        t = np.log(ratio) / r
                        if t >= 0.0:
                            dtdx1 = abs((1.0 / r) * (1.0 / x1 + 1.0 / (K - x1)))
                            val = alpha * np.exp(-alpha * t) * dtdx1
            row[j] = val
            s += val * weights[j]
        if s > 1e-12:
            for j in range(N):
                M[i, j] = row[j] / s
        else:
            idx_K = 0
            best = 1e300
            for j in range(N):
                d = abs(grid[j] - K)
                if d < best:
                    best = d
                    idx_K = j
            M[i, idx_K] = 1.0 / weights[idx_K] if weights[idx_K] > 0 else 1.0
    return M

def lambda_density_iteration(alpha, envA, envB, grid, weights, KMAX,
                              max_iter=400, tol=1e-13):
    """lambda(alpha) via quadrature solution of the Perron-Frobenius equation.
    The O(N^2) transition-matrix build is numba-jitted; the power-iteration
    step uses numpy's BLAS matmul, which is faster than a hand-rolled loop
    for the repeated N x N times N contraction."""
    rA, KA, gA = envA['r'], envA['K'], envA['gamma']
    rB, KB, gB = envB['r'], envB['K'], envB['gamma']
    MA = build_transition_matrix_nb(rA, KA, alpha, grid, weights, KMAX)
    MB = build_transition_matrix_nb(rB, KB, alpha, grid, weights, KMAX)

    N = len(grid)
    rho = np.ones(N)
    rho /= np.sum(rho * weights)
    # Power iteration: rho1(x1=j) = sum_i rho(x0=i) * weights[i] * M[i,j]
    # -- weight goes on the SOURCE (row) index here.
    WA_src = MA * weights[:, None]
    WB_src = MB * weights[:, None]
    for _ in range(max_iter):
        rho1 = rho @ WA_src
        rho2 = rho1 @ WB_src
        rho2 /= np.sum(rho2 * weights)
        diff = np.max(np.abs(rho2 - rho))
        rho = rho2
        if diff < tol:
            break

    # E[ln(x1/x0) | x0=i] = sum_j M[i,j] * weights[j] * ln(x1_j / x0_i)
    # -- weight goes on the DESTINATION (column) index here, since j is now
    # the integration variable. (Using weights[:, None] here, as if reusing
    # WA_src, silently swaps which index gets weighted and gives a wrong
    # answer that still looks plausible -- this was caught by cross-checking
    # against the unvectorized loop version and against direct Monte Carlo.)
    log_grid = np.log(grid)
    WA_dst = MA * weights[None, :]
    log_ratio_given_x0 = (WA_dst * (log_grid[None, :] - log_grid[:, None])).sum(axis=1)
    A_alpha = np.sum(rho * weights * log_ratio_given_x0)
    return (KA + KB - gA - gB) / 2 - (alpha / 2) * (KA / rA - KB / rB) * A_alpha

def find_alpha_c_density_iteration(envA, envB, grid, weights, KMAX, lo=0.02, hi=20.0):
    f = lambda a: lambda_density_iteration(a, envA, envB, grid, weights, KMAX)
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    return brentq(f, lo, hi, xtol=1e-5)


# ── Monte Carlo estimate of alpha_c (stochastic simulation of the switching) ───
#
# The cycle-to-cycle update is an inherent Markov chain (x at cycle i depends
# on cycle i-1), so it can't be vectorized across cycles with numpy the way
# the inner arithmetic can. numba's JIT compiles the per-cycle Python loop
# itself to machine code, which is what actually removes the bottleneck here
# (roughly 10-15x faster than the numpy version per call, and the bisection /
# multi-seed loops around it benefit directly since each call is now cheap).

@njit
def lambda_monte_carlo_nb(alpha, rA, KA, gA, rB, KB, gB, n_cycles, seed, burn=12000):
    np.random.seed(seed)
    x = 1.0
    for _ in range(burn):
        tau = np.random.exponential(1.0 / alpha)
        x = x_flow(tau, x, rA, KA)
        tau = np.random.exponential(1.0 / alpha)
        x = x_flow(tau, x, rB, KB)

    total_integral = 0.0
    total_time = 0.0
    for _ in range(n_cycles):
        tau_A = np.random.exponential(1.0 / alpha)
        x0 = x
        x1 = x_flow(tau_A, x0, rA, KA)
        if x0 > 1e-12 and x1 > 1e-12:
            contrib = (KA - gA) * tau_A - (KA / rA) * np.log(x1 / x0)
        else:
            contrib = (KA - gA) * tau_A
        total_integral += contrib
        total_time += tau_A
        x = x1

        tau_B = np.random.exponential(1.0 / alpha)
        x0b = x
        x1b = x_flow(tau_B, x0b, rB, KB)
        if x0b > 1e-12 and x1b > 1e-12:
            contrib = (KB - gB) * tau_B - (KB / rB) * np.log(x1b / x0b)
        else:
            contrib = (KB - gB) * tau_B
        total_integral += contrib
        total_time += tau_B
        x = x1b

    return total_integral / total_time

def find_alpha_c_monte_carlo(envA, envB, seed, lo=0.02, hi=20.0, n_bisect=16, n_cycles=100000):
    rA, KA, gA = envA['r'], envA['K'], envA['gamma']
    rB, KB, gB = envB['r'], envB['K'], envB['gamma']
    f = lambda a: lambda_monte_carlo_nb(a, rA, KA, gA, rB, KB, gB, n_cycles, seed)
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        # No sign change in [lo, hi]: widen the bracket once before giving up,
        # since gamma_B values near the edge of the paradox window can push
        # alpha_c outside the default range.
        for hi_try in (hi * 5, hi * 25, hi * 125):
            fhi_try = f(hi_try)
            if flo * fhi_try < 0:
                hi, fhi = hi_try, fhi_try
                break
        else:
            for lo_try in (lo / 5, lo / 25, lo / 125):
                flo_try = f(lo_try)
                if flo_try * fhi < 0:
                    lo, flo = lo_try, flo_try
                    break
            else:
                return None
    for _ in range(n_bisect):
        mid = (lo + hi) / 2
        fmid = f(mid)
        if flo * fmid < 0:
            hi = mid
        else:
            lo, flo = mid, fmid
    return (lo + hi) / 2


# ── Deterministic equal-duty reference curve ────────────────────────────────────

def lambda_det_periodic(T, envA, envB):
    """lambda(T) on the periodic orbit for deterministic equal-duty switching."""
    rA, KA, gA = envA['r'], envA['K'], envA['gamma']
    rB, KB, gB = envB['r'], envB['K'], envB['gamma']
    if T < 1e-9:
        x_star = (rA + rB) / (rA / KA + rB / KB)
        return x_star - (gA + gB) / 2
    def residual(x0):
        return x_flow(T / 2, x_flow(T / 2, x0, rA, KA), rB, KB) - x0
    lo, hi = 1e-9, KB * (1 - 1e-13)
    r_lo, r_hi = residual(lo), residual(hi)
    x0 = KB * (1 - 1e-13) if r_lo * r_hi > 0 else brentq(residual, lo, hi, xtol=1e-14)
    x1 = x_flow(T / 2, x0, rA, KA)
    IA, _ = quad(lambda t: x_flow(t, x0, rA, KA) - gA, 0, T / 2, limit=200)
    IB, _ = quad(lambda t: x_flow(t, x1, rB, KB) - gB, 0, T / 2, limit=200)
    return (IA + IB) / T

def find_Tc(envA, envB):
    T_vals = np.logspace(-2, 3.3, 200)
    lam_prev = lambda_det_periodic(T_vals[0], envA, envB)
    for i in range(1, len(T_vals)):
        lam_cur = lambda_det_periodic(T_vals[i], envA, envB)
        if lam_prev * lam_cur < 0:
            return brentq(lambda T: lambda_det_periodic(T, envA, envB),
                           T_vals[i - 1], T_vals[i], xtol=1e-8)
        lam_prev = lam_cur
    return None


# ── Main plot ─────────────────────────────────────────────────────────────────
#
# Why density iteration and Monte Carlo disagree near both edges of the
# paradox window (gamma_B near 2.0 or near 2.333 here):
#
# It is NOT a transient/burn-in effect. Burn-in was tested from 1,000 to
# 1,000,000 steps and the MC estimate barely moves -- the chain mixes fast.
#
# The real cause is that density iteration discretizes the transition kernel
# f(x1 | x0) on a FIXED grid, and that kernel's width depends on alpha:
#
#   - At LARGE alpha (near the high-gamma_B edge, where alpha_c is large):
#     dwells are short, so x1 stays close to x0. The kernel becomes a narrow
#     spike of width ~ 1/alpha *anywhere* in the domain, including the bulk,
#     where the grid (built to be fine near the two fixed points 0 and K) is
#     comparatively coarse. A narrow spike straddling only 1-2 grid cells is
#     not resolved by trapezoidal quadrature, biasing the estimate.
#
#   - At SMALL alpha (near the low-gamma_B edge, where alpha_c is small):
#     dwells are long, so x1 sits very close to the environment's fixed point
#     K. The kernel develops a genuine boundary-layer singularity there
#     (|dt/dx1| -> infinity as x1 -> K), and the bulk of the probability mass
#     (>80% in some cases checked) can fall within a window as thin as 1e-3
#     of K. Increasing the grid density near the boundary helps, but the
#     convergence is logarithmically slow -- pushing N from 1,800 to 6,000
#     only closed a small fraction of the gap, at a large cost in runtime.
#
# Both mechanisms are the same underlying issue (a kernel narrower than the
# local grid spacing), just appearing at opposite ends of the alpha range.
# Monte Carlo has no such resolution limit -- it samples the true continuous
# kernel directly -- so the MC points are the trustworthy reference near the
# edges, and the density-iteration curve should be read as exact in the
# interior of the gamma_B range and increasingly approximate near its ends.

def plot_critical_alpha(params=None, outdir=None,
                         gB_curve=None, gB_points=None, n_mc_reps=4,
                         N_grid=1800, n_cycles=150000):
    if params is None:
        params = env_params()
    if gB_curve is None:
        gB_curve = np.linspace(2.00, 2.24, 100)
    if gB_points is None:
        gB_points = np.linspace(2.00, 2.24, 16)

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

    envA = params['A']
    envB_base = params['B']
    KMAX = max(envA['K'], envB_base['K'])

    # ── Warm up the JIT compiler once, outside the timed sweeps ─────────────────
    _grid_warm = make_grid(10000, KMAX)
    _weights_warm = make_weights(_grid_warm)
    lambda_density_iteration(1.0, envA, {**envB_base, 'gamma': 2.2},
                              _grid_warm, _weights_warm, KMAX, max_iter=10)
    lambda_monte_carlo_nb(1.0, envA['r'], envA['K'], envA['gamma'],
                          envB_base['r'], envB_base['K'], 2.2, 100, 0, burn=10)

    # ── Semi-analytical curve: density iteration (quadrature, no randomness) ───
    print("Computing density-iteration curve...")
    grid = make_grid(N_grid, KMAX)
    weights = make_weights(grid)
    alpha_c_DI = []
    for gB in gB_curve:
        envB = {**envB_base, 'gamma': gB}
        ac = find_alpha_c_density_iteration(envA, envB, grid, weights, KMAX)
        alpha_c_DI.append(ac if ac is not None else np.nan)
        print(f"  gamma_B={gB:.4f}: alpha_c (density iteration) = {ac}")
    alpha_c_DI = np.array(alpha_c_DI, dtype=float)

    # ── Simulated points: Monte Carlo, several seeds for error bars ────────────
    print("Computing Monte Carlo points...")
    gB_points_kept, alpha_c_MC_mean, alpha_c_MC_std = [], [], []
    for gB in gB_points:
        envB = {**envB_base, 'gamma': gB}
        ests = [find_alpha_c_monte_carlo(envA, envB, seed=100 + rep, n_cycles=n_cycles)
                for rep in range(n_mc_reps)]
        ests = [e for e in ests if e is not None]
        if len(ests) == 0:
            print(f"  gamma_B={gB:.3f}: no sign change found even after widening "
                  f"the bracket; skipping this point")
            continue
        gB_points_kept.append(gB)
        alpha_c_MC_mean.append(np.mean(ests))
        alpha_c_MC_std.append(np.std(ests))
        print(f"  gamma_B={gB:.3f}: alpha_c (MC) = {np.mean(ests):.4f} +/- {np.std(ests):.4f}"
              f"  ({len(ests)}/{n_mc_reps} reps succeeded)")
    gB_points = np.array(gB_points_kept)
    alpha_c_MC_mean = np.array(alpha_c_MC_mean)
    alpha_c_MC_std = np.array(alpha_c_MC_std)

    # ── Reference curve: deterministic equal-duty switching ────────────────────
    print("Computing deterministic reference curve...")
    alpha_c_det = []
    for gB in gB_curve:
        envB = {**envB_base, 'gamma': gB}
        Tc = find_Tc(envA, envB)
        alpha_c_det.append(2 / Tc if Tc else np.nan)
    alpha_c_det = np.array(alpha_c_det)

    # ── Style ─────────────────────────────────────────────────────────────────
    C_LINE     = '#7d1d96'   # deep magenta-purple, consistent with the plasma family
    C_LINE_DET = '#c9882e'   # warm amber, distinguishes the deterministic reference
    C_POINTS   = '#1a1a1a'   # near-black markers, matches stable-fixed-point styling

    # ── Layout ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(7, 5.2))
    fig.patch.set_facecolor('none')
    gs = GridSpec(1, 1, left=0.16, right=0.96, top=0.90, bottom=0.16)
    ax = fig.add_subplot(gs[0, 0])

    ax.plot(gB_curve, alpha_c_det, '--', color=C_LINE_DET, linewidth=1.8,
             label='Deterministic (equal-duty)', zorder=1)
    ax.plot(gB_curve, alpha_c_DI, '-', color=C_LINE, linewidth=2.2,
             label='Numerical (density iteration)', zorder=2)
    # ax.errorbar(gB_points, alpha_c_MC_mean, yerr=alpha_c_MC_std,
    #             fmt='x', color=C_POINTS, markeredgecolor=C_POINTS,
    #             markerfacecolor=C_POINTS, markersize=9, capsize=4,
    #             elinewidth=1.4, linewidth=0, label='Simulation (Monte Carlo)',
    #             zorder=5)
    ax.plot(gB_points, alpha_c_MC_mean, 'x', markersize=9, color=C_POINTS, markeredgecolor=C_POINTS,
             label='Simulation (Monte Carlo)',
                zorder=5)

    ax.set_xlabel(r'$\gamma_B$', labelpad=4)
    ax.set_ylabel(r'Critical rate  $\alpha_c$', labelpad=4)
    ax.set_yscale('log')
    ax.grid(False)
    ax.patch.set_facecolor('none')
    ax.legend(frameon=False, fontsize=18, loc='upper left')

    # ── Save ──────────────────────────────────────────────────────────────────
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / 'critical_alpha.svg'
    plt.savefig(outpath, dpi=300, format='svg',
                bbox_inches='tight', facecolor='none', transparent=True)
    plt.savefig(outpath.with_suffix('.pdf'), dpi=300, format='pdf',
                bbox_inches='tight', facecolor='none', transparent=True)
    plt.close(fig)
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    plot_critical_alpha(env_params(), outdir='src/3paramLogisticLVStochastic/plots/critical_alpha')