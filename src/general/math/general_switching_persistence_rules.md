# Rigorous Mathematical Theory of Switching Rescue

## 1. Goal and Mathematical Setup

We consider a system tracking prey $x(t)$ and predator $y(t)$ switching between two environments $\sigma(t) \in \{1, 2\}$:

$$
\frac{dx}{dt} = F_{\sigma(t)}(x,y), \qquad \frac{dy}{dt} = G_{\sigma(t)}(x,y)
$$

**Assumption 1 (Invariant Extinction Boundary):** 
In both environments, if $y=0$, it remains equal to $0$. Hence:
$$ G_1(x,0) = 0 \quad \text{and} \quad G_2(x,0) = 0 \quad \text{for all } x \ge 0. $$
This allows us to easily study the purely-prey dynamics:
$$ f_i(x) := F_i(x, 0) $$

**Assumption 2 (Stable Single-Environment Equilibria):**
Assume each isolated environment $i$ has a stable positive prey equilibrium $x = K_i$.
$$ f_i(K_i) = 0 \quad \text{and} \quad f_i'(K_i) < 0 $$
Without loss of generality, assume $K_1 < K_2$.

**Assumption 3 (Fixed-Environment Extinction):**
The predator $y$ goes extinct in each isolated environment. The initial linear growth of $y$ along the boundary is given by $h_i(x) := \partial_y G_i(x,0)$. So the extinction premise is exactly:
$$ h_1(K_1) < 0 \quad \text{and} \quad h_2(K_2) < 0 $$

We are looking for necessary and sufficient conditions *on these functions directly* such that $y$ can survive when $\sigma(t)$ alternates between 1 and 2.

---

## 2. Deriving the Periodic Switching Exponent

Let us translate the dynamics of $y$ directly into spatial integrals of $F$ and $G$. 

Suppose we construct a periodic switching cycle between bounds $x_a$ and $x_b$, such that $K_1 \le x_a < x_b \le K_2$. 
- The cycle switches to environment $2$ at $x_a$. Since $x < K_2$, the prey grows ($f_2(x) > 0$), driving the system up to $x_b$.
- At $x_b$, the cycle switches to environment $1$. Since $x > K_1$, the prey shrinks ($f_1(x) < 0$), driving the system back down to $x_a$.

While near $y=0$, the growth of $y$ is exactly $\frac{dy}{dt} = h_{\sigma(t)}(x) y$. 
We compute the net logarithmic gain $\Delta \ln y$ over one full cycle:

1. **Upward phase (Env 2):** 
   Since $dt = \frac{dx}{f_2(x)}$, the growth is $\int \limits_{x_a}^{x_b} h_2(x) dt = \int \limits_{x_a}^{x_b} \frac{h_2(x)}{f_2(x)} dx$.
2. **Downward phase (Env 1):** 
   Since $dt = \frac{dx}{f_1(x)}$, the growth is $\int \limits_{x_b}^{x_a} h_1(x) dt$. 
   Reversing the integration limits multiplies by $-1$: $\int \limits_{x_b}^{x_a} \frac{h_1(x)}{f_1(x)} dx = - \int \limits_{x_a}^{x_b} \frac{h_1(x)}{f_1(x)} dx$.

**Total Gain per Cycle:**
Summing the phases gives the central rigorous measure of switching rescue:

$$ \Delta \ln y = \int_{x_a}^{x_b} \underbrace{\left[ \frac{\partial_y G_2(x,0)}{F_2(x,0)} - \frac{\partial_y G_1(x,0)}{F_1(x,0)} \right]}_{:=\, R(x)} dx $$

Here, $R(x)$ is the **Switching Rescue Function**, fully determined by fixed functions $F$ and $G$. 

---

## 3. The Necessary and Sufficient Condition

Using the derivation above, we can prove highly specific conditions for coexistence.

### Theorem 1 (Existence of Rescue)
*Given a system satisfying Assumptions 1–3, there exists a periodic switching cycle that rescues species $y$ if and only if there exists some interval inside $(K_1, K_2)$ where the integral of $R(x)$ is strictly positive.*

**Proof:**
- **Sufficiency:** If such an interval $[x_a, x_b]$ exists, we enforce a local switching law triggered at exactly $x_a$ and $x_b$. By the derivation in Section 2, the per-cycle gain of $y$ near zero is uniformly positive ($\Delta \ln y > 0$). Floquet theory ensures $y$ initially grows radially outward, implying survival (permanence) assuming global boundedness.
- **Necessity:** If $R(x) \le 0$ for all $x \in (K_1, K_2)$, then for *any* interval $x_a < x_b$, the integral $\int_{x_a}^{x_b} R(x) dx \le 0$. Since all periodic trajectories are inherently forced to oscillate bounds between $K_1$ and $K_2$, the net-gain per cycle generated is uniformly $\le 0$. Hence, $y$ inevitably tracks toward extinction and no cycle can possibly rescue it. $\blacksquare$

### Corollary 1: A Strict Null Condition
If the rescue function is uniformly negative or zero everywhere:
$$ \frac{\partial_y G_2(x,0)}{F_2(x,0)} \le \frac{\partial_y G_1(x,0)}{F_1(x,0)} \quad \text{for all } x \in (K_1, K_2) $$
then **no period switching between these environments can save $y$**.

---

## 4. Connecting the Limit to Fast Switching

A very fast switching cycle means $x_b - x_a \to 0$. In this limit, the system hovers indefinitely close to an averaged equilibrium $x^*$. 

Does the Rescue Function $R(x)$ successfully predict the fast switching conditions?

### Theorem 2 (Fast Switching Threshold)
*Fast switching with duty cycle $p$ yields an averaged fixed point $x^*$. The predator $y$ survives at this limit if and only if $R(x^*) > 0$.*

**Proof derived rigorously from F and G:**
In the high-frequency limit with duty cycle $p \in (0,1)$, the prey sits at an averaged fixed point $x^*$ where spatially continuous movement ceases:
$$ p f_1(x^*) + (1-p) f_2(x^*) = 0 \implies \frac{p}{1-p} = -\frac{f_2(x^*)}{f_1(x^*)} $$
At this fixed point, the average linear growth rate for $y$ over infinite time simplifies to:
$$ \bar{H}(x^*) = p h_1(x^*) + (1-p) h_2(x^*) $$

To analyze when $\bar{H}(x^*) > 0$, notice $1-p > 0$ enabling safe factor division:
$$ \frac{p}{1-p}h_1(x^*) + h_2(x^*) > 0 $$

Substitute the temporal ratio $\frac{p}{1-p}$ previously derived:
$$ -\frac{f_2(x^*)}{f_1(x^*)} h_1(x^*) + h_2(x^*) > 0 $$

Because $x^*$ requires both vectors to push against one another towards local equilibrium, $x^*$ structurally falls tightly between bounded coordinates $K_1$ and $K_2$. Therefore, necessarily $f_2(x^*) > 0$ and $f_1(x^*) < 0$. 
Because structurally $f_2(x^*) > 0$, we can safely divide the entire resultant inequality by $f_2(x^*)$:
$$ -\frac{h_1(x^*)}{f_1(x^*)} + \frac{h_2(x^*)}{f_2(x^*)} > 0 $$

Rearranging yields:
$$ \frac{h_2(x^*)}{f_2(x^*)} - \frac{h_1(x^*)}{f_1(x^*)} > 0 \implies R(x^*) > 0 $$

**Conclusion of Proof:** 
The analytical requirement that the time-averaged growth rate $\bar{H}(x^*)$ is positive is fundamentally an **algebraic restatement** of evaluating the spatial switching Rescue Function $R(x)$ at $x^*$.

---

## 5. Summary Mathematical Checklist

By transforming evaluations over time $t$ logically to operations over distance $x$, the entire generalized concept condenses securely underneath a single mathematical entity: $R(x)$.

**Calculate the Switching Rescue Function exactly:**
$$ R(x) = \frac{\partial_y G_2(x,0)}{F_2(x,0)} - \frac{\partial_y G_1(x,0)}{F_1(x,0)} $$

The mathematical constraints dictating rescue behavior:
- **Rule 1 (Dead System):** If $R(x) \le 0$ everywhere between the equilibria, switching fundamentally cannot induce persistence.
- **Rule 2 (Rescuable System):** If $R(x) > 0$ somewhere between the equilibria, any positional bounds $[x_a, x_b]$ entirely inside this positive region structurally guarantees mathematically derived survival for species $y$. 
   - A wide interval integrates positive growth for moderate-speed dynamic cycles. 
   - Limiting the coordinate bounds directly onto a focal point $x^*$ equates mathematically identically to high-frequency fast-switching survival.