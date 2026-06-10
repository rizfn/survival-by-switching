"""
Phase-space limit cycle plots + R(x) shading.

Two models:
  Base LV:   dx/dt = rx(1-x/K) - xy,        dy/dt = xy - gamma*y
  Shielding: dx/dt = rx(1-x/K) - xye^{-ax}, dy/dt = xye^{-ax} - y

Standard case parameters (from paper):
  rA=1, KA=1, gA=1.01,  rB=4, KB=2, gB=2.2
  (gA ~ KA so predator barely excluded in A; gB > KB so excluded in B)

Shielding case parameters (as specified):
  r=1, gamma=1 in both environments
  KA=1, KB=2, alphaA=0, alphaB=0.35

R(x) formula (paper eq. 41):
  R(x) = h2(x)/f2(x) - h1(x)/f1(x)
  where hE(x) = d/dy GE(x,0)  [linearised predator growth near y=0]
        fE(x) = FE(x,0)       [prey dynamics on y=0 axis]

  Base LV:   hE(x) = x - gE,      fE(x) = rE*x*(1 - x/KE)
  Shielding: hE(x) = x*e^{-aE*x} - 1, fE(x) = rE*x*(1 - x/KE)

The integration domain for R is [x*_A, x*_B] where x*_E is the
prey equilibrium on the y=0 axis in environment E:
  Base LV:   x*_E = gE   (from dy/dt=0 at small y: x=gamma)
  Shielding: x*_E = KE   (prey goes to K when predator absent)

Wait — more carefully:
  The paper assumes a switching cycle where prey moves between x*_A and x*_B
  driven by the two environments. For R(x)>0 to be meaningful, we need the
  *prey nullcline fixed points on y=0*:
    Base LV: prey fixed point on y=0 is x=K (logistic), but the predator
             nullcline (x=gamma) determines where y can grow. The cycle
             bounces between gamma_A and gamma_B on y=0 approximately.
    Shielding: prey fixed point on y=0 is x=K; predator nullcline is
               x*e^{-alpha*x}=1.

  For the shading, we shade x in [min(x*_A,x*_B), max(x*_A,x*_B)] where
  x*_E is the x-coordinate at which the predator nullcline crosses y=0,
  i.e. hE(x*_E)=0:
    Base LV:   x*_E = gE
    Shielding: x*_E solves x*e^{-a*x}=1  (=1 for both since gamma=1)
               alphaA=0 -> x*_A = 1
               alphaB=0.35 -> solve numerically
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import brentq

# ── RK4 integrator ────────────────────────────────────────────────────────────

def make_rhs(rA, KA, gA, alphaA, rB, KB, gB, alphaB):
    def rhs(x, y, env):
        if env == 'A':
            r, K, g, a = rA, KA, gA, alphaA
        else:
            r, K, g, a = rB, KB, gB, alphaB
        interaction = x * y * np.exp(-a * x)
        dx = r * x * (1.0 - x / K) - interaction
        dy = interaction - g * y
        return dx, dy
    return rhs

def rk4_step(x, y, t, T, p, rhs_fn, dt):
    env = lambda tt: 'A' if (tt % T) < p * T else 'B'
    e = env(t)
    k1x, k1y = rhs_fn(x,            y,            e)
    k2x, k2y = rhs_fn(x+.5*dt*k1x,  y+.5*dt*k1y,  env(t+.5*dt))
    k3x, k3y = rhs_fn(x+.5*dt*k2x,  y+.5*dt*k2y,  env(t+.5*dt))
    k4x, k4y = rhs_fn(x+dt*k3x,     y+dt*k3y,     env(t+dt))
    xn = max(x + (dt/6)*(k1x+2*k2x+2*k3x+k4x), 1e-14)
    yn = max(y + (dt/6)*(k1y+2*k2y+2*k3y+k4y), 1e-14)
    return xn, yn

def simulate_limit_cycle(T, rhs_fn, p=0.5,
                          x0=0.5, y0=0.05,
                          dt=5e-4, transient=2000, measure_periods=8):
    x, y, t = x0, y0, 0.0
    # transient
    n_trans = int(transient / dt)
    for _ in range(n_trans):
        x, y = rk4_step(x, y, t, T, p, rhs_fn, dt)
        t += dt
    # record one limit cycle
    xs, ys = [x], [y]
    n_meas = int(measure_periods * T / dt)
    for _ in range(n_meas):
        x, y = rk4_step(x, y, t, T, p, rhs_fn, dt)
        t += dt
        xs.append(x)
        ys.append(y)
    return np.array(xs), np.array(ys)


# ── R(x) helper ───────────────────────────────────────────────────────────────

def R_func(x, rA, KA, gA, alphaA, rB, KB, gB, alphaB):
    """
    R(x) = h2(x)/f2(x) - h1(x)/f1(x)
    hE(x) = x*exp(-aE*x) - gE   [= x - gE for base LV with a=0]
    fE(x) = rE*x*(1 - x/KE)
    """
    h1 = x * np.exp(-alphaA * x) - gA
    h2 = x * np.exp(-alphaB * x) - gB
    f1 = rA * x * (1.0 - x / KA)
    f2 = rB * x * (1.0 - x / KB)
    with np.errstate(divide='ignore', invalid='ignore'):
        t1 = np.where(np.abs(f1) > 1e-12, h1 / f1, np.nan)
        t2 = np.where(np.abs(f2) > 1e-12, h2 / f2, np.nan)
    return t2 - t1


# ── Main plot function ────────────────────────────────────────────────────────

def make_plot(rA, KA, gA, alphaA, rB, KB, gB, alphaB,
              T_values, title, outpath,
              x_range=(0.0, 3.0), y_range=(0.0, 0.8),
              dt=5e-4, transient=2000, measure_periods=8):

    rhs_fn = make_rhs(rA, KA, gA, alphaA, rB, KB, gB, alphaB)

    # R(x) shading domain: everywhere.

    fig, ax = plt.subplots(figsize=(7, 5.5))

    # ── R(x) shading ──────────────────────────────────────────────────────────
    x_shade = np.linspace(max(x_range[0] * 0.98, 1e-4), x_range[1] * 1.02, 2000)
    R_vals  = R_func(x_shade, rA, KA, gA, alphaA, rB, KB, gB, alphaB)
    y_top   = y_range[1]

    ax.fill_between(x_shade, 0, y_top,
                    where=(np.isfinite(R_vals) & (R_vals > 0)),
                    alpha=0.20, color='#2ca02c', zorder=0,
                    label=r'$R(x)>0$')
    ax.fill_between(x_shade, 0, y_top,
                    where=(np.isfinite(R_vals) & (R_vals <= 0)),
                    alpha=0.10, color='#d62728', zorder=0,
                    label=r'$R(x)\leq 0$')

    # ── Limit cycles ──────────────────────────────────────────────────────────
    cmap   = plt.get_cmap('plasma')
    T_log  = np.log10(T_values)
    T_norm = (T_log - T_log.min()) / (T_log.max() - T_log.min() + 1e-12)

    for i, T in enumerate(T_values):
        print(f"  T = {T}")
        # Use the midpoint between the two carrying capacities as the initial condition for the limit cycle simulation.
        x0 = 0.5 * (KA + KB)
        xs, ys = simulate_limit_cycle(T, rhs_fn, p=0.5,
                                       x0=x0, y0=0.02,
                                       dt=dt, transient=transient,
                                       measure_periods=measure_periods)
        color = cmap(0.12 + 0.76 * T_norm[i])
        ax.plot(xs, ys, color=color, lw=1.6, alpha=0.88,
                label=f'$T={T}$', zorder=3)

    # Fixed points on y=0 axis
    for xfp, lbl, xoff in [(KA, r'$K_A$', -0.05), (KB, r'$K_B$', +0.05)]:
        ax.scatter([xfp], [0], marker='o', color='black', s=50, zorder=5,
                   clip_on=False)
        ax.text(xfp + xoff, -0.025, lbl, ha='center', va='top', fontsize=10)

    # ── Styling ───────────────────────────────────────────────────────────────
    ax.set_xlabel('Prey $x$', fontsize=13)
    ax.set_ylabel('Predator $y$', fontsize=13)
    ax.set_title(title, fontsize=12, pad=8)
    ax.set_xlim(x_range)
    ax.set_ylim(y_range[0] - 0.04, y_top)
    ax.tick_params(labelsize=11)
    ax.grid(True, linestyle='--', alpha=0.22)
    ax.set_axisbelow(True)

    green_patch = mpatches.Patch(color='#2ca02c', alpha=0.45, label=r'$R(x)>0$')
    red_patch   = mpatches.Patch(color='#d62728', alpha=0.30, label=r'$R(x)\leq 0$')
    handles, labels = ax.get_legend_handles_labels()
    cyc_h = [h for h, l in zip(handles, labels) if l.startswith('$T')]
    cyc_l = [l for l in labels if l.startswith('$T')]
    ax.legend([green_patch, red_patch] + cyc_h,
              [r'$R(x)>0$', r'$R(x)\leq 0$'] + cyc_l,
              fontsize=9, frameon=True, framealpha=0.85,
              loc='upper right')

    fig.tight_layout()
    fig.savefig(outpath, dpi=220, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {outpath}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import os
    os.makedirs('src/general/plots', exist_ok=True)

    T_values = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]

    # ── 1. Base LV (paper parameters) ────────────────────────────────────────
    # rA=1, KA=1, gA=1.01  ->  predator barely excluded (x*_A = gA = 1.01 > KA=1)
    # rB=4, KB=2, gB=2.2   ->  predator excluded       (x*_B = gB = 2.2  > KB=2)
    # Limit cycle bounces between x*_A ~ 1 and x*_B ~ 2.2 on the y=0 boundary
    print("Base LV:")
    make_plot(
        rA=1.0, KA=1.0, gA=1.01, alphaA=0.0,
        rB=4.0, KB=2.0, gB=2.2,  alphaB=0.0,
        T_values=T_values,
        title=r'Base LV:  $\alpha_A=\alpha_B=0$',
        outpath='src/general/plots/phase_space_base_LV.pdf',
        x_range=(0.0, 2.8),
        y_range=(0.0, 0.55),
        transient=3000,
    )

    # ── 2. Shielding model ────────────────────────────────────────────────────
    # r=1, gamma=1 in both environments; environments differ in K and alpha.
    # alphaA=0:    x*_A solves x*e^0 = 1  ->  x*_A = 1  (= KA, so marginally excluded)
    # alphaB=0.35: x*_B solves x*e^{-0.35x} = 1  (numerically ~ 1.58)
    #              KB=2 > x*_B, so predator is excluded in B (x*_B < KB)
    print("\nShielding:")
    make_plot(
        rA=1.0, KA=1.0, gA=1.0, alphaA=0.0,
        rB=1.0, KB=2.0, gB=1.0, alphaB=0.35,
        T_values=T_values,
        title=r'Shielding:  $\alpha_A=0,\ \alpha_B=0.35$',
        outpath='src/general/plots/phase_space_shielding.pdf',
        x_range=(0.0, 2.5),
        y_range=(0.0, 0.55),
        transient=3000,
    )