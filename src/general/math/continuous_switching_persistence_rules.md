# Switching-induced persistence in a *continuously* changing environment

**A step-by-step derivation, generalizing Noronha & Kaneko (2026).**

---

## 0. The result, stated up front (no sup/inf needed)

We have one "loser" species $y$ (the predator) that goes extinct in **every** frozen
environment. The environment is now a single scalar $s(t)$ that varies **continuously**
in time and drives all parameters. The question: can some periodic protocol $s(t)$ keep
$y$ alive?

The clean answer is:

> **A continuously-varying environment can rescue $y$ if and only if there exist two
> environments $s_a, s_b$ in the accessible continuum whose *ordinary two-environment
> rescue function* (the paper's $R$) is positive somewhere. The optimal protocol is
> bang-bang between exactly those two environments.**

So the continuum buys you **nothing new except the freedom to choose which two
environments to alternate.** Everything reduces to the paper's discrete criterion applied
to a freely chosen pair. The "envelope" function $R_{\text{cont}}$ that appears below is
just an efficient way to *search all pairs at once*; it is not a new mechanism, and we
prove it equals the best-pair criterion exactly.

Everything below builds this statement carefully, tracking every sign.

---

## 1. Model and assumptions

Two species $x$ (prey) and $y$ (predator), one scalar environment $s(t)$:

$$
\dot x = F\big(x,y;\,s(t)\big), \qquad
\dot y = G\big(x,y;\,s(t)\big).
$$

All parameters are smooth functions of $s$. For the concrete Lotka–Volterra model of the
paper,

$$
F(x,y;s) = r(s)\,x\Big(1-\tfrac{x}{K(s)}\Big) - xy, \qquad
G(x,y;s) = xy - \gamma(s)\,y .
$$

We work on and near the extinction axis $y=0$, so define the two **boundary fields**

$$
\boxed{\,f_s(x) := F(x,0;s) = r(s)\,x\Big(1-\tfrac{x}{K(s)}\Big)\,}
\qquad
\boxed{\,h_s(x) := \partial_y G(x,0;s) = x-\gamma(s)\,}.
$$

$f_s(x)$ is the prey's speed when the predator is absent; $h_s(x)$ is the predator's
*per-capita* growth rate when it is vanishingly rare.

We keep the paper's three assumptions, imposed **at every frozen value of $s$**:

- **(A1) The axis is invariant.** $G(x,0;s)=0$ for all $x,s$. The predator cannot appear
  from nothing. (For LV: $G(x,0;s)=x\cdot 0-\gamma\cdot 0=0$. ✓)
- **(A2) Each frozen environment has a boundary fixed point.** There is a $K_s$ with
  $f_s(K_s)=0$ and $f_s'(K_s)<0$ (stable along the axis). (For LV: $K_s=K(s)$. ✓)
- **(A3) The predator dies in every frozen environment.** The per-capita growth at the
  prey fixed point is negative: $h_s(K_s)<0$, i.e. $\gamma(s)>K_s$ for all $s$.

Assumption (A3) is the whole point: there is **no** value of $s$ at which $y$ survives.
We ask whether *moving* $s$ can rescue it.

---

## 2. "Extinction" is not a point — it's a moving orbit

Because the axis $y=0$ is invariant (A1), the extinct state is the prey running on the
boundary by itself:

$$
\dot x = f_{s(t)}(x).
$$

The prey chases the moving fixed point $K_{s(t)}$. If $s(t)$ is periodic with period $T$,
then after transients $x(t)$ settles onto a $T$-periodic trajectory $x^\star(t)$, with
$x^\star(t+T)=x^\star(t)$. This is a limit cycle living on the boundary. Our question
becomes: **is this boundary orbit stable or unstable in the transverse ($y$) direction?**
If a tiny bit of predator grows, $y$ persists; if it shrinks, $y$ dies.

---

## 3. Transverse stability = a Floquet exponent

Put a tiny predator perturbation on the boundary orbit and linearize the $y$-equation
about $y=0$. Taylor-expand $G$ in $y$ at fixed $x$:

$$
G(x,y;s) = \underbrace{G(x,0;s)}_{=0\ \text{by (A1)}}
          + y\,\underbrace{\partial_y G(x,0;s)}_{=\,h_s(x)}
          + \mathcal{O}(y^2).
$$

The constant term vanishes exactly (A1), which is why this linearization is clean: there
is no forcing, only multiplicative growth. Along the orbit,

$$
\dot y = h_{s(t)}\big(x^\star(t)\big)\,y + \mathcal{O}(y^2).
$$

This is a scalar linear ODE with a time-periodic coefficient. Separating variables,
$\dfrac{d}{dt}\ln y = h_{s(t)}(x^\star(t))$, and integrating over one full period:

$$
\ln\frac{y(T)}{y(0)}
= \int_0^T h_{s(t)}\big(x^\star(t)\big)\,dt \;=:\; \Delta\ln y .
$$

Define the **invasion (Floquet) exponent** $\Lambda := \Delta\ln y / T$. Then

$$
\boxed{\;y \text{ persists } \iff \Lambda>0 \iff \Delta\ln y=\int_0^T h_{s(t)}(x^\star(t))\,dt>0.\;}
$$

Nothing here required $s(t)$ to be a square wave. This is the exact criterion for *any*
continuous protocol.

---

## 4. Change variables: from time $t$ to prey density $x$

Doing the $t$-integral directly is hard because $x^\star(t)$ has no closed form for a
general protocol. The trick (which is what makes the whole theory tractable) is to
integrate over $x$ instead of $t$. On the boundary $\dot x = f_{s(t)}(x)$, so

$$
dt = \frac{dx}{f_{s(t)}(x)}.
$$

Define the **per-density weight**

$$
\boxed{\,w(x;s) := \frac{h(x;s)}{f(x;s)} = \frac{\text{predator per-capita growth at }x}
{\text{prey boundary speed at }x}\,}.
$$

Its meaning: $w\,dx = h\,dt$ is the little slice of predator log-growth the orbit picks up
while the prey passes through $[x,x+dx]$. Slow passage (small $f$) means a lot of time
there, hence a big slice — this is why $w$ has $f$ in the denominator.

---

## 5. The per-cycle log-gain and the rescue function

Assume for now the orbit is **single-humped**: over one period $x^\star$ rises from a
minimum $x_a$ to a maximum $x_b$ (the **up-sweep**) and falls back (the **down-sweep**).
We handle the general multi-hump case rigorously in §7; single-humped is all we need for
the main result, and the optimal protocol turns out to be single-humped anyway.

**Up-sweep** ($x:\,x_a\to x_b$, so $dx>0$). Here $x$ is increasing, so $f>0$, which by (A2)
means the current carrying capacity satisfies $K_s>x$. Call this environment $U$:

$$
\big[\Delta\ln y\big]_{\uparrow}
= \int_{\text{up}} h\,dt
= \int_{x_a}^{x_b} \frac{h_U(x)}{f_U(x)}\,dx
= \int_{x_a}^{x_b} w\big(x;s_\uparrow(x)\big)\,dx,
\qquad f_U>0.
$$

**Down-sweep** ($x:\,x_b\to x_a$, so $dx<0$). Here $x$ decreasing, $f<0$, meaning $K_s<x$.
Call it $D$. Flipping the limits to run $x_a\to x_b$ introduces a minus sign:

$$
\big[\Delta\ln y\big]_{\downarrow}
= \int_{x_b}^{x_a} \frac{h_D(x)}{f_D(x)}\,dx
= -\int_{x_a}^{x_b} w\big(x;s_\downarrow(x)\big)\,dx,
\qquad f_D<0.
$$

Adding the two branches:

$$
\boxed{\;\Delta\ln y
= \int_{x_a}^{x_b}\Big[\,w\big(x;s_\uparrow(x)\big) - w\big(x;s_\downarrow(x)\big)\,\Big]\,dx
= \int_{x_a}^{x_b} R_{\text{prot}}(x)\,dx.\;}
$$

The **realized rescue function** of a given protocol,

$$
R_{\text{prot}}(x)
= \frac{h\big(x;s_\uparrow(x)\big)}{f\big(x;s_\uparrow(x)\big)}
- \frac{h\big(x;s_\downarrow(x)\big)}{f\big(x;s_\downarrow(x)\big)},
$$

has *exactly* the structure of the paper's rescue function — a difference of two $h/f$'s.
The only new feature is that the two environments being compared are whichever ones the
protocol visits on the up- vs down-sweep at that particular $x$, and they may drift with
$x$.

---

## 6. Recovering the paper (discrete two-environment case)

If $s(t)$ is a two-value square wave alternating environments $1$ and $2$ with
$K_1<K_2$, then on the up-sweep the prey is being pulled toward the *higher* fixed point,
so $s_\uparrow\equiv 2$; on the down-sweep, toward the lower, so $s_\downarrow\equiv 1$.
Then

$$
R_{\text{prot}}(x)=\frac{h_2(x)}{f_2(x)}-\frac{h_1(x)}{f_1(x)}
= \underbrace{\frac{\partial_y G_2(x,0)}{F_2(x,0)}-\frac{\partial_y G_1(x,0)}{F_1(x,0)}}_{=\,R(x)\text{ of the paper (Eq. 13)}},
$$

and $\Delta\ln y=\int_{x_a}^{x_b}R(x)\,dx$. The paper's theorem — *rescue is possible iff
$R(x)>0$ on some subinterval of $[K_1,K_2]$* — is the special case of our formula in which
only two environments are accessible. Good: the continuous theory contains the discrete
one as a limit, not as a separate story.

---

## 7. General (multi-hump) protocols: the crossing-count bound

To make the necessity ("if no pair works, nothing works") airtight, drop the
single-hump assumption. For an arbitrary periodic boundary orbit, fix a level $x$ and
count how the orbit crosses it. On a periodic orbit the ascending and descending
crossings alternate and are equal in number; call it $m(x)$. Ascending crossings have
$f>0$ (so $K_s>x$), descending have $f<0$ (so $K_s<x$). Summing the $w\,dx$ slices over
all crossings,

$$
\Delta\ln y = \int \Bigg[\;\sum_{i=1}^{m(x)} w\big(x;s_i^\uparrow\big)
- \sum_{j=1}^{m(x)} w\big(x;s_j^\downarrow\big)\;\Bigg]\,dx .
$$

Now bound each term by the best it could possibly be. Every ascending environment obeys
$K_s>x$, so $w(x;s_i^\uparrow)\le \sup_{\,K_s>x} w(x;s)=:W_+(x)$. Every descending one obeys
$K_s<x$, so $w(x;s_j^\downarrow)\ge \inf_{\,K_s<x} w(x;s)=:W_-(x)$. Hence

$$
\sum_i w_i^\uparrow-\sum_j w_j^\downarrow \le m(x)\big[W_+(x)-W_-(x)\big],
$$

and therefore, for **every** protocol whatsoever,

$$
\boxed{\;\Delta\ln y \;\le\; \int m(x)\,\big[\,W_+(x)-W_-(x)\,\big]\,dx
= \int m(x)\,R_{\text{cont}}(x)\,dx,\;}
$$

where we have defined the **envelope rescue function**

$$
\boxed{\;R_{\text{cont}}(x) := \sup_{\,s:\,K_s>x}\frac{h(x;s)}{f(x;s)}
\;-\; \inf_{\,s:\,K_s<x}\frac{h(x;s)}{f(x;s)}.\;}
$$

Because $m(x)\ge 0$, the bound gives immediately:

> **If $R_{\text{cont}}(x)\le 0$ for all $x$, then $\Delta\ln y\le 0$ for every protocol —
> no continuous environment can rescue $y$.** (Necessity.)

Extra humps ($m>1$) only add more crossings; they can never beat the per-crossing
envelope. This is *why* wiggling the environment cleverly can't outperform a clean
two-environment alternation.

---

## 8. The clean condition: rescue ⟺ some pair works

The sup/inf looks heavy, but it collapses to a pair statement. Suppose
$R_{\text{cont}}(x_0)>0$ at some $x_0$. Let $s_\star$ (nearly) attain the sup and
$s_{\star\star}$ (nearly) attain the inf at $x_0$. These are two specific environments with
$K_{s_{\star\star}}<x_0<K_{s_\star}$, and by construction

$$
R_{\text{cont}}(x_0) = \underbrace{w(x_0;s_\star)-w(x_0;s_{\star\star})}_{\text{pairwise } R \text{ of the pair } (s_{\star\star},\,s_\star)} > 0 .
$$

So $R_{\text{cont}}(x_0)>0$ **is** the paper's pairwise $R$ of one particular pair, evaluated
at $x_0$. Conversely any pair's $R$ is $\le R_{\text{cont}}$ pointwise (the pair competes in
the same sup/inf). Therefore

$$
\boxed{\;\exists x:\ R_{\text{cont}}(x)>0
\;\;\Longleftrightarrow\;\;
\exists\ \text{a pair } (s_a,s_b)\subset\text{continuum with }R_{s_a s_b}(x)>0 \text{ somewhere}.\;}
$$

**Numerically this equivalence is exact:** computing $R_{\text{cont}}$ two ways — as the
sup/inf envelope, and as a brute-force maximum over all pairs — agrees to machine zero
($\max|\Delta| = 0$; see verification script). So the honest, un-ugly statement of the
theorem is the pair version in §0. The envelope is merely the search.

**Sufficiency.** If some pair's $R$ is positive on an interval, the paper already proves a
rescuing bang-bang cycle exists for that pair (choose $x_a,x_b$ inside the positive
interval and the switching rate above the pair's critical value). Combined with §7's
necessity, we have a genuine *if and only if*.

---

## 9. Why the optimal rescuer is bang-bang between two environments

The bound $\Delta\ln y\le\int m(x)R_{\text{cont}}(x)\,dx$ is saturated by taking $m(x)=1$
(a single hump) and placing the up-sweep in the sup-environment and the down-sweep in the
inf-environment. If a single pair $(s_a,s_b)$ attains the sup and inf across the whole
positive interval — which happens whenever the extremes of $h/f$ are achieved at fixed
environments, as in the standard LV parameters — then the maximizer is literally a
two-state square wave. Continuously sweeping through intermediate environments can only
match, never beat it, and adding humps strictly wastes orbit length. Hence:

> The best a continuum can do is pick the single most favorable pair and alternate it
> bang-bang. The continuum's only advantage over a fixed two-environment system is the
> **freedom to choose that pair**.

This is the resolution of the "ugliness": the sup/inf is real, but it is optimized by a
two-point (bang-bang) control, so the final mechanism is as simple as the discrete one.

---

## 10. Worked check on the Lotka–Volterra numbers

Standard parameters: env 1 $(K_1,r_1,\gamma_1)=(1,1,1.1)$, env 2 $(K_2,r_2,\gamma_2)=(2,4,2.1)$,
$K_1<K_2$. Both die in isolation: $\gamma_1=1.1>K_1=1$ and $\gamma_2=2.1>K_2=2$. ✓ (A3).

For $x\in(1,2)$ the only valid pair is $U=$ env 2 (needs $K>x$: $K_2=2>x$) and $D=$ env 1
(needs $K<x$: $K_1=1<x$). At $x=1.5$:

$$
w_U = \frac{h_2}{f_2}=\frac{1.5-2.1}{4\cdot 1.5\,(1-1.5/2)}=\frac{-0.6}{1.5}=-0.400,
\qquad f_U=+1.5>0\ \checkmark
$$
$$
w_D = \frac{h_1}{f_1}=\frac{1.5-1.1}{1\cdot 1.5\,(1-1.5/1)}=\frac{0.4}{-0.75}=-0.5333,
\qquad f_D=-0.75<0\ \checkmark
$$
$$
R(1.5)=w_U-w_D=-0.400-(-0.5333)=\mathbf{+0.1333}>0 .
$$

Positive, so rescue is possible — consistent with the paper's finite critical rate
$\alpha_c\approx0.716$. The maximum over $x\in(1,2)$ is $\max_x R = +0.13395$, attained near
$x\approx1.53$. Both signs came out exactly as the branch orientation of §5 predicts.

---

## 11. Numerical confirmation

Three independent checks (script `continuous_switching.py`, verification in `verify.py`):

1. **Sign / hand-calc.** $R(1.5)=+0.13333$ analytically = code. ✓
2. **Envelope = best pair.** $\max_x|R_{\text{cont}}-R_{\text{best pair}}|=0$ to machine
   precision. ✓ (§8 equivalence.)
3. **Interior-pair counterexample** (proving the continuum is *strictly* more powerful than
   endpoints alone). Construct a bent path
   $s\!\in\![0,0.5]$: env 1 $\to$ good env $(K,r,\gamma)=(2,4,2.1)$;
   $s\!\in\![0.5,1]$: good env $\to$ bad env $(2,0.5,2.5)$. Every frozen environment dies
   ($\gamma(s)>K(s)$ throughout). The **endpoint** pair $s\in\{0,1\}$ gives
   $\max_x R=-4.70<0$ (cannot rescue), yet $\max_x R_{\text{cont}}=+0.134>0$: the rescuing
   pair is the **interior** environment $s=0.5$ against $s=0$. Full 2-D simulation confirms
   the endpoint protocol drives $y\to10^{-15}$ (extinction) while the interior protocol
   holds $y\approx0.033$ (persistence).

Figures (publication style, parameters in filenames):
`R_of_x_std_K1_1_K2_2_r1_1_r2_4_g1_1.1_g2_2.1`,
`interior_pair_rescue_r_1_4_0.5_K_1_2_2_g_1.1_2.1_2.5`,
`lambda_vs_period_squareness_k3_k6_k12_std_params`.

---

## 12. The one genuine subtlety: dwell distribution, not just rate

The identity $\Delta\ln y=\int_{x_a}^{x_b}R\,dx$ looks rate-independent because the time
variable was integrated out. It is **not**, and here is exactly why: the turning points
$x_a,x_b$ *themselves depend on the protocol*. The weight $w=h/f$ has integrable
singularities at $x=K_s$ (where $f=0$). As switching slows, the orbit lingers near a
frozen fixed point $K_i$ for a long real time with $h(K_i)=K_i-\gamma_i<0$; equivalently
$x_a\to K_1$ and $x_b\to K_2$ push the integral into those singular endpoints where the
negative contributions dominate. This is why **slow switching is always lethal** and why a
finite critical rate exists.

Consequence for protocol *shape*: $R_{\text{cont}}(x)>0$ guarantees a rescuer *exists*, but
a given protocol realizes it only if it approaches **bang-bang** — long dwells at the two
chosen environments, fast transitions between. A smooth **sinusoid** does the opposite
(lingers at its extremes, rushes through the middle), behaves like slow switching, and
never rescued the standard parameters at any period tested. Quantitatively, smoothing a
square wave as $s(t)\propto\tanh(k\sin\omega t)$ and increasing the squareness $k$ drives
the critical rate down toward the paper's bang-bang value:

$$
k=3:\ \alpha_c\approx1.44 \quad\to\quad k=6:\ 0.99 \quad\to\quad k=12:\ 0.85
\quad\longrightarrow\quad \alpha_c^{\text{paper}}=0.716 .
$$

So the paper's discrete result is precisely the $k\to\infty$ limit of the continuous
theory, and smoother protocols must cycle faster (or fail).

---

## Summary

- The transverse stability of the extinct boundary orbit is governed exactly by
  $\Delta\ln y=\int_0^T h_{s(t)}(x^\star(t))\,dt$ (Floquet, §3).
- Changing variables $t\to x$ turns this into $\int_{x_a}^{x_b} R_{\text{prot}}(x)\,dx$
  with $R_{\text{prot}}=w_\uparrow-w_\downarrow$, $w=h/f$ (§4–5) — the paper's rescue
  function, now with freely chosen environments on each branch.
- A crossing-count bound gives, for *every* protocol,
  $\Delta\ln y\le\int m(x)\,R_{\text{cont}}(x)\,dx$, so rescue is possible **iff**
  $R_{\text{cont}}(x)>0$ somewhere (§7), which is **iff some pair of environments satisfies
  the paper's ordinary criterion** (§8, verified exact).
- The optimal rescuer is bang-bang between that single best pair (§9); the continuum's
  only gift is the freedom to choose the pair.
- Rate still matters through the turning points $x_a,x_b$: slow/smooth loses, and the
  discrete $\alpha_c$ is the bang-bang limit (§12).

**Answer to the original question:** yes, a continuously-varying environment admits a clean
switching-rescue condition — it is the paper's condition, optimized over all pairs drawn
from the continuum, with a rigorous envelope $R_{\text{cont}}$ and no genuinely new
mechanism beyond pair selection.