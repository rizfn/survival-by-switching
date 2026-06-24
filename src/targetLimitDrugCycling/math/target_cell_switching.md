# Switching-induced persistence in the target-cell-limited model

A second, structurally distinct instance of the switching-rescue mechanism, built on the
standard within-host (virus-dynamics) model of Nowak & May and Perelson. Unlike the
Lotka–Volterra case, the two extinction fixed points here are **not** set by a carrying
capacity but by a **source–sink balance** between target-cell production and turnover.
We work out only the two analytically clean regimes: slow switching (each drug in
isolation) and infinitely fast switching (the time-averaged drug).

---

## 1. The model

The standard target-cell-limited model has three variables: uninfected target cells $T$,
infected cells $I$, and free pathogen $V$. The environment $\sigma(t)\in\{1,2\}$ labels
which of two drug regimens is currently applied, and each regimen fixes a set of
parameters:

$$
\frac{dT}{dt} = \lambda_\sigma - d_\sigma\,T - \beta_\sigma\,T V,\qquad
\frac{dI}{dt} = \beta_\sigma\,T V - \delta_\sigma\,I,\qquad
\frac{dV}{dt} = p_\sigma\,I - c_\sigma\,V .
$$

The terms are the textbook ones: target cells are produced at constant rate
$\lambda_\sigma$ and die at rate $d_\sigma$; infection is the mass-action contact
$\beta_\sigma TV$; infected cells die at $\delta_\sigma$; pathogen is produced at
$p_\sigma$ and cleared at $c_\sigma$. The pathogen (the "loser") is the pair $(I,V)$.

Note immediately that $I=V=0$ is invariant: if the pathogen is absent it cannot
reappear. This is the within-host meaning of the first general-framework assumption
($G_i(x,0)=0$, no recovery from extinction).

### 1.1 Fast-pathogen reduction to two dimensions

Free pathogen turns over much faster than cells ($c_\sigma$ large), so we apply a
quasi-steady-state approximation to $V$. Setting $dV/dt \approx 0$ gives

$$
V \approx \frac{p_\sigma}{c_\sigma}\,I .
$$

Substitute this into the infection term, $\beta_\sigma T V \approx
\frac{\beta_\sigma p_\sigma}{c_\sigma}\,T I$, and define the single effective
**infectivity**

$$
a_\sigma \equiv \frac{\beta_\sigma\, p_\sigma}{c_\sigma}.
$$

The model collapses to a two-variable system, one per environment:

$$
\boxed{\;\frac{dT}{dt} = \lambda_\sigma - d_\sigma\,T - a_\sigma\,T I,\qquad
\frac{dI}{dt} = \big(a_\sigma\,T - \delta_\sigma\big)\,I\;}
\tag{1}
$$

This is the object we analyse. It has the consumer–resource *skeleton* of
Lotka–Volterra (a bilinear $TI$ coupling and a threshold predator), but the resource
term is the **linear source–sink** $\lambda_\sigma - d_\sigma T$, not a logistic
$rT(1-T/K)$. That single difference is what makes this a genuinely different model
rather than a relabelling.

Two combinations of parameters will do all the work:

$$
T^{0}_\sigma \equiv \frac{\lambda_\sigma}{d_\sigma}\quad\text{(pathogen-free setpoint)},
\qquad
T^{*}_\sigma \equiv \frac{\delta_\sigma}{a_\sigma}\quad\text{(invasion threshold)} .
$$

$T^0_\sigma$ is the target-cell level the host relaxes to with no infection; $T^*_\sigma$
is the target-cell level the pathogen needs to break even (its net per-capita growth is
$a_\sigma T-\delta_\sigma$, which is zero at $T=T^*_\sigma$).

---

## 2. Slow switching: each drug acts as a fixed environment

When the switching period is very long ($T_{\text{period}}\to\infty$), the system spends
long stretches in a single environment and equilibrates there. So the slow-switching
question is simply: *does each drug, held constant, clear the infection?* We answer it by
finding the fixed points of (1) for a fixed $\sigma$ and testing their stability. I drop
the subscript $\sigma$ within this subsection.

### 2.1 Fixed points

Set the right-hand sides of (1) to zero.

**Pathogen-free fixed point.** Put $I=0$. Then $dI/dt=0$ automatically, and
$dT/dt = \lambda - dT = 0$ gives

$$
\big(T,I\big) = \Big(\tfrac{\lambda}{d},\,0\Big) = \big(T^{0},\,0\big).
$$

**Endemic fixed point.** Put $I\neq 0$. Then $dI/dt=0$ forces $aT-\delta=0$, i.e.
$T=\delta/a=T^{*}$. Substituting into $dT/dt=0$:

$$
\lambda - d\,T^{*} - a\,T^{*} I = 0
\;\Longrightarrow\;
I = \frac{\lambda - d\,T^{*}}{a\,T^{*}}
= \frac{\lambda}{\delta} - \frac{d}{a}.
$$

This is biologically meaningful ($I>0$) only if $\lambda/\delta > d/a$, which we will see
is exactly the condition for the pathogen to survive.

### 2.2 Stability of the pathogen-free state

Linearise (1) about a general point. The Jacobian is

$$
J(T,I)=
\begin{pmatrix}
-d - aI & -aT\\[2pt]
aI & aT-\delta
\end{pmatrix}.
$$

Evaluate at the pathogen-free point $(T^{0},0)=(\lambda/d,\,0)$:

$$
J\big(T^{0},0\big)=
\begin{pmatrix}
-d & -a\lambda/d\\[2pt]
0 & a\lambda/d-\delta
\end{pmatrix}.
$$

Because this matrix is upper-triangular, its eigenvalues are the diagonal entries:

$$
\mu_1 = -d \;(<0),\qquad
\mu_2 = \frac{a\lambda}{d}-\delta .
$$

The first eigenvalue is always negative (target cells relax back to $T^0$). The
pathogen-free state is therefore stable — **the drug clears the infection** — precisely
when $\mu_2<0$:

$$
\frac{a\lambda}{d}-\delta < 0
\;\Longleftrightarrow\;
\boxed{\,R_0 \equiv \dfrac{a\,\lambda}{d\,\delta} < 1\,}.
\tag{2}
$$

This $R_0$ is the standard within-host basic reproduction number; restoring the original
parameters, $R_0 = \beta p \lambda/(c d \delta)$, the textbook expression. It has a clean
geometric reading:

$$
R_0 = \frac{a\lambda}{d\delta}
= \frac{\lambda/d}{\delta/a}
= \frac{T^{0}}{T^{*}} ,
$$

so **clearance ($R_0<1$) means the pathogen-free setpoint sits below the invasion
threshold**, $T^{0}<T^{*}$. (When $R_0<1$ the endemic $I=\lambda/\delta-d/a$ is negative,
i.e. not feasible, consistent with extinction being the only outcome.)

### 2.3 The slow-switching condition

For the pathogen to be a "loser in both drugs", we require (2) in each environment:

$$
\boxed{\,R_0^{(1)}=\dfrac{a_1\lambda_1}{d_1\delta_1}<1
\quad\text{and}\quad
R_0^{(2)}=\dfrac{a_2\lambda_2}{d_2\delta_2}<1\,}
\tag{3}
$$

equivalently $T^0_1<T^*_1$ and $T^0_2<T^*_2$.

---

## 3. Infinitely fast switching: the time-averaged drug

Now take the opposite limit, $T_{\text{period}}\to 0$, with a fraction $p$ of each period
spent in environment 1 and $(1-p)$ in environment 2. Over one (vanishingly short) period
the state $(T,I)$ barely changes, so integrating (1) over a period gives

$$
\Delta T \approx T_{\text{period}}\Big[p\,f_1(T,I) + (1-p)\,f_2(T,I)\Big],
$$

and likewise for $I$, where $f_\sigma$ are the right-hand sides. Dividing by
$T_{\text{period}}$, the dynamics follow the **time-averaged vector field** — this is
standard fast-switching averaging. Writing the $p$-weighted means

$$
\bar\lambda = p\lambda_1+(1-p)\lambda_2,\quad
\bar d = p d_1+(1-p)d_2,\quad
\bar a = p a_1+(1-p)a_2,\quad
\bar\delta = p\delta_1+(1-p)\delta_2,
$$

the averaged system is again a target-cell model:

$$
\frac{dT}{dt} = \bar\lambda - \bar d\,T - \bar a\,T I,\qquad
\frac{dI}{dt} = \big(\bar a\,T - \bar\delta\big) I .
\tag{4}
$$

We take the symmetric duty cycle $p=\tfrac12$ from here on (matching the main-text
50–50 analysis); a general $p$ just keeps the weights above.

### 3.1 Persistence under fast switching

System (4) is identical in form to (1), so we can reuse Section 2 with barred
parameters. The pathogen-free point is $(\bar\lambda/\bar d,\,0)$, and the pathogen now
**invades and persists** when its eigenvalue is *positive*:

$$
\frac{\bar a\,\bar\lambda}{\bar d}-\bar\delta > 0
\;\Longleftrightarrow\;
\boxed{\,\bar R_0 \equiv \dfrac{\bar a\,\bar\lambda}{\bar d\,\bar\delta} > 1\,}.
\tag{5}
$$

When (5) holds, the averaged system carries a feasible, stable endemic equilibrium

$$
T^\star = \frac{\bar\delta}{\bar a},\qquad
I^\star = \frac{\bar\lambda}{\bar\delta}-\frac{\bar d}{\bar a}
= \frac{\bar a\bar\lambda-\bar d\bar\delta}{\bar a\,\bar\delta} > 0 .
$$

(The 2-D target-cell model is a competitive system whose pathogen-free and endemic states
exchange stability transcritically at $\bar R_0=1$, so $\bar R_0>1$ implies trajectories
approach the endemic state — the pathogen survives.)

---

## 4. The combined criterion and why it can be met

Putting (3) and (5) together, **switching rescues a pathogen that each drug eliminates**
whenever

$$
\boxed{\;\frac{a_1\lambda_1}{d_1\delta_1}<1,\qquad
\frac{a_2\lambda_2}{d_2\delta_2}<1,\qquad
\frac{(a_1{+}a_2)(\lambda_1{+}\lambda_2)}{(d_1{+}d_2)(\delta_1{+}\delta_2)}>1 \;}
\tag{6}
$$

(the last expression is $\bar R_0$ at $p=\tfrac12$, since $\bar a=\tfrac12(a_1{+}a_2)$,
etc.).

**Why these are compatible.** Write $R_0 = u\,v$ with $u=\lambda/d$ (the setpoint $T^0$)
and $v=a/\delta$ (the inverse threshold $1/T^*$). Each drug has $R_0^{(i)}=u_i v_i<1$. The
fast-switching number, however, is **not** the average of the $R_0^{(i)}$: it is

$$
\bar R_0 = \frac{\bar\lambda}{\bar d}\cdot\frac{\bar a}{\bar\delta}
= \underbrace{\Big(\tfrac{\lambda_1+\lambda_2}{d_1+d_2}\Big)}_{\text{between }u_1,\,u_2}
\;\cdot\;
\underbrace{\Big(\tfrac{a_1+a_2}{\delta_1+\delta_2}\Big)}_{\text{between }v_1,\,v_2}.
$$

Each factor is a *ratio of sums*, which by the mediant inequality lies between the two
environments' values. If the two factors are **anti-correlated across drugs** — drug 1
gives a high setpoint $u$ but a low $v$ (hard for the pathogen), drug 2 the reverse — then
each product $u_i v_i$ is kept small (high $\times$ low), while the averaged factors are
both middling, and their product (middling $\times$ middling) can exceed $1$. This is the
ratio-of-means versus mean-of-ratios gap; it is the target-cell version of the same
Parrondo-type inequality that drives the Lotka–Volterra and SIS+Q examples.

This also makes precise the caveat for drug cycling: if the two drugs differ **only** in
their direct effect on the pathogen ($a$, i.e. $\beta,p,c$) and leave target-cell
turnover ($\lambda,d$) untouched, then $T^0_1=T^0_2$, the two factors stop being
independent, and no rescue is possible. **Rescue requires the drugs to leave the host
(target-cell) environment in genuinely different states**, not merely to kill the
pathogen by different amounts.

---

## 5. A concrete example

| parameter | $\lambda$ | $d$ | $a=\beta p/c$ | $\delta$ | $T^0=\lambda/d$ | $T^*=\delta/a$ | $R_0=T^0/T^*$ |
|---|---|---|---|---|---|---|---|
| **drug 1** | 4.0 | 1.0 | 1.0 | 5.0 | 4.0 | 5.0 | **0.80** |
| **drug 2** | 0.2 | 1.0 | 5.0 | 1.25 | 0.2 | 0.25 | **0.80** |
| **average** | 2.1 | 1.0 | 3.0 | 3.125 | 2.1 | 1.04 | **2.02** |

Both drugs clear the infection ($R_0^{(i)}=0.8<1$), yet fast 50–50 alternation gives
$\bar R_0\approx 2.0>1$ and a stable endemic load $I^\star\approx 0.34$.

**Intuition (transient overshoot).** Drug 1 *fills the tank*: it drives target cells
toward the high setpoint $T^0_1=4$, but it is itself non-permissive, with threshold
$T^*_1=5>T^0_1$. Drug 2 is highly permissive ($T^*_2=0.25$) but *starves the tank*,
pulling $T$ toward $T^0_2=0.2$. Because target cells relax on the timescale $1/d=1$ while
switching is fast, $T$ never falls to drug 2's meagre setpoint; it lingers near the level
stocked by drug 1, which is $\sim 16\times$ above drug 2's threshold. During every drug-2
phase the pathogen then grows at rate $a_2(T-T^*_2)\approx 5(4-0.25)\approx 19$. The
pathogen feeds on drug 1's target cells under drug 2's permissive conditions — with no
resistance, tolerance, or mutation.

Numerical integration of (1) confirms the criterion: under either drug alone the infected
population decays to extinction, while under fast switching it climbs to the predicted
endemic load.

---

## 6. Relation to the general criterion

System (1) satisfies the three assumptions of the general framework: (i) $I=0$ is
invariant; (ii) each environment has a stable pathogen-free fixed point at
$T^0_\sigma=\lambda_\sigma/d_\sigma$ (the eigenvalue $-d_\sigma<0$); and (iii) the pathogen
is excluded there, $h_\sigma(T^0_\sigma)=a_\sigma T^0_\sigma-\delta_\sigma<0$, i.e.
$R_0^{(\sigma)}<1$. With $f_\sigma(T)=\lambda_\sigma-d_\sigma T$ and
$h_\sigma(T)=a_\sigma T-\delta_\sigma$, the switching rescue function is

$$
R(T) = \frac{a_2 T-\delta_2}{\lambda_2-d_2 T} - \frac{a_1 T-\delta_1}{\lambda_1-d_1 T},
$$

and a periodic cycle can rescue the pathogen iff $R(T)>0$ somewhere on
$[\,T^0_1,\,T^0_2\,]$ — the same geometric statement as in the Lotka–Volterra case,
now with the resource self-dynamics $f_\sigma$ being source–sink rather than logistic.
