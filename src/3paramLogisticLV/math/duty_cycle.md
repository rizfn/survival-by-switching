# Duty Cycle Analysis: Unequal Time Switching

If the system spends a fraction $p$ of its time in Environment A and $(1-p)$ of its time in Environment B, the dynamics and fixed points scale smoothly. This unequal time weighting is known as the **duty cycle**.

## 1. The Averaged Vector Field Equations

In the fast-switching limit, the system follows an averaged vector field weighted by the duty cycle $p \in [0, 1]$:

$$
\begin{align}
\frac{dx}{dt} &= p \left[ r_A x \left(1 - \frac{x}{K_A}\right) - x y \right] + (1-p) \left[ r_B x \left(1 - \frac{x}{K_B}\right) - x y \right] \\
\frac{dy}{dt} &= p (x y - \gamma_A y) + (1-p) (x y - \gamma_B y)
\end{align}
$$

Grouping the terms, we define the parameter averages as functions of $p$:
*   $\bar{r}(p) = p r_A + (1-p) r_B$
*   $\bar{\gamma}(p) = p \gamma_A + (1-p) \gamma_B$
*   $\bar{c}(p) = p \frac{r_A}{K_A} + (1-p) \frac{r_B}{K_B}$  (the averaged prey self-competition term)

The fast-switching system is then exactly a new, effective Lotka-Volterra interaction:
$$
\begin{align}
\frac{dx}{dt} &= \bar{r}(p) x - \bar{c}(p) x^2 - x y \\
\frac{dy}{dt} &= y \big( x - \bar{\gamma}(p) \big)
\end{align}
$$

## 2. Shifted Fixed Points and Persistence

Setting the derivatives to zero yields the new fixed points parameterized by $p$.

**Extinction Fixed Point ($y=0$):**
$$ x = \bar{x}_{fast}(p) = \frac{\bar{r}(p)}{\bar{c}(p)} = \frac{p r_A + (1-p) r_B}{p \frac{r_A}{K_A} + (1-p) \frac{r_B}{K_B}} $$

**Coexistence Fixed Point ($y>0$):**
The predator equation forces the prey perfectly to the averaged death rate:
$$ x^*(p) = \bar{\gamma}(p) $$

Substituting this back into the prey equation to find the steady-state predator density:
$$ y^*(p) = \bar{r}(p) - \bar{\gamma}(p) \bar{c}(p) $$

### The Transcritical Persistence Condition
The condition for survival is $y^*(p) > 0$. This implies:
$$ \bar{r}(p) - \bar{\gamma}(p) \bar{c}(p) > 0 \implies \frac{\bar{r}(p)}{\bar{c}(p)} > \bar{\gamma}(p) $$
$$ \bar{x}_{fast}(p) > \bar{\gamma}(p) $$

Therefore, the persistence condition dynamically shifts with $p$. The Floquet exponent at the $y=0$ manifold simply adapts to:
$$ \lambda_{fast}(p) = \bar{x}_{fast}(p) - \bar{\gamma}(p) $$

For the predator to survive under fast switching, there must exist some range of $p$ strictly inside $(0, 1)$ where $\lambda_{fast}(p) > 0$.

## 3. Optimal Duty Cycle for Persistence

Is there an "optimal" duty cycle? It depends on your metric, but a natural choice is minimizing extinction risk by **maximizing the steady-state predator density $y^*(p)$**.

Let's look at the shape of $y^*(p)$ structurally:
$$ y^*(p) = \bar{r}(p) - \bar{\gamma}(p) \bar{c}(p) $$

Notice that $\bar{r}(p)$, $\bar{\gamma}(p)$, and $\bar{c}(p)$ are all strictly **linear functions** of $p$. 
Therefore, the product $\bar{\gamma}(p) \bar{c}(p)$ creates a quadratic term in $p$. 

This means that $y^*(p)$ is a **parabola**. 
Because we assume the predator goes extinct in each isolated environment:
*   At $p = 0$ (all Env B): $y^*(0) < 0$
*   At $p = 1$ (all Env A): $y^*(1) < 0$

If there is any window of survival where $y^*(p) > 0$, the parabola must open downwards (it is concave). Because the boundary points are both negative, a strictly positive arc necessitates an absolute **peak**. 

Taking the derivative to find the optimum:
$$ \frac{d}{dp} y^*(p) = \frac{d\bar{r}}{dp} - \left( \frac{d\bar{\gamma}}{dp}\bar{c}(p) + \bar{\gamma}(p)\frac{d\bar{c}}{dp} \right) = 0 $$

To evaluate this comprehensively, let's redefine the averages explicitly in terms of the differences between the environments.
Let $\Delta r = r_A - r_B$, $\Delta \gamma = \gamma_A - \gamma_B$, and $\Delta c = \frac{r_A}{K_A} - \frac{r_B}{K_B}$. We also denote the Environment B baselines as $r_B$, $\gamma_B$, and $c_B = \frac{r_B}{K_B}$.

Now we can write the linear averages as:
1. $\bar{r}(p) = p \Delta r + r_B \implies \frac{d\bar{r}}{dp} = \Delta r$
2. $\bar{\gamma}(p) = p \Delta \gamma + \gamma_B \implies \frac{d\bar{\gamma}}{dp} = \Delta \gamma$
3. $\bar{c}(p) = p \Delta c + c_B \implies \frac{d\bar{c}}{dp} = \Delta c$

Substitute these directly back into our derivative equation:
$$ \Delta r - \Big[ \Delta \gamma (p \Delta c + c_B) + (p \Delta \gamma + \gamma_B) \Delta c \Big] = 0 $$

Next, distribute the terms inside the parentheses:
$$ \Delta r - \Big[ p \Delta \gamma \Delta c + c_B \Delta \gamma + p \Delta \gamma \Delta c + \gamma_B \Delta c \Big] = 0 $$

Combine the identical $p \Delta \gamma \Delta c$ terms:
$$ \Delta r - \Big[ 2p \Delta \gamma \Delta c + c_B \Delta \gamma + \gamma_B \Delta c \Big] = 0 $$

Now, distribute the negative sign and isolate the term containing $p$:
$$ \Delta r - c_B \Delta \gamma - \gamma_B \Delta c = 2p \Delta \gamma \Delta c $$

Finally, dividing by $2 \Delta \gamma \Delta c$ yields the exact closed-form solution for the optimum duty cycle $p_{opt}^{y^*}$:
$$ p_{opt}^{y^*} = \frac{\Delta r - c_B \Delta \gamma - \gamma_B \Delta c}{2 \Delta \gamma \Delta c} $$

Does max $y^*$ minimize the critical switching frequency?

If we rely on hand-waving, generating the largest possible steady-state population $y^*$ implies we are maximizing the "safety margin" against extinction, so this should allow for the maximum possible switching period $T_{crit}$ (or minimum frequency $f_{crit}$).

However, mathematically, this is **not strictly true**. There are three distinct definitions of "optimal," and they yield slightly different duty cycles.

### A. Maximizing $y^*$ vs. Maximizing $\lambda_{fast}$
The fast-switching invasion exponent $\lambda_{fast}(p)$ determines whether the predator can grow from zero. It is related to $y^*(p)$ by:
$$ \lambda_{fast}(p) = \bar{x}_{fast}(p) - \bar{\gamma}(p) = \frac{y^*(p)}{\bar{c}(p)} $$
Taking the derivative of $\lambda_{fast}(p)$ with respect to $p$ gives:
$$ \frac{d\lambda_{fast}}{dp} = \frac{\frac{dy^*}{dp} \bar{c}(p) - y^*(p) \frac{d\bar{c}}{dp}}{\bar{c}(p)^2} $$
At $p_{opt}^{y^*}$, we know $\frac{dy^*}{dp} = 0$. Therefore, the slope of the invasion exponent at this point is:
$$ \left. \frac{d\lambda_{fast}}{dp} \right|_{p=p_{opt}^{y^*}} = - \frac{y^* \Delta c}{\bar{c}(p)^2} $$
Unless $\Delta c = 0$ (meaning the prey self-competition $r_i/K_i$ is identical in both environments), this derivative is non-zero. 
**Conclusion 1:** The duty cycle that maximizes the initial invasion speed ($\lambda_{fast}$) is distinct from the duty cycle that maximizes the final resting population ($y^*$).

### B. Minimizing Critical Frequency ($f_{crit}$)
The finite critical switching frequency is found by setting the exact periodic invasion condition to zero. From general theory, $f_{crit}$ is bounded by evaluating the integral of the rescue function $R(x)$ over the exact limit cycle bounds $x_a(p)$ and $x_b(p)$.

The limit cycle bounds $x_a, x_b$ span asymmetrically depending on the nonlinear prey dynamics. Minimizing $f_{crit}(p)$ (maximizing $T_{crit}$) requires stretching this interval to exactly cover the widest possible region where $R(x) > 0$. 

While $\lambda_{fast}(p)$ governs the *infinitesimal* limit ($x_a \to x_b \to \bar{x}_{fast}$), $f_{crit}$ incorporates the *global* nonlinearity of the prey vector field $f(x)$ away from the fixed point. The $p$ that centers the fast-switching state over the optimal position for infinitesimal growth will generally not perfectly balance against the finite-time distortions tracked by the spatial integrals $\int \frac{h(x)}{f(x)} dx$.

**Final Conclusion:** 
While $p_{opt}^{y^*}$, $p_{opt}^{\lambda}$, and $p_{opt}^{T_{crit}}$ are generally very close to one another—because they all rely on leveraging the positive portions of the rescue function $R(x)$—they are formally mathematically distinct quantities. The hand-waving intuition provides a solid biological approximation, but an exact mathematical equivalence between max $y^*$ and min $f_{crit}$ is false unless strict symmetric constraints on the environments are imposed
When we shift the duty cycle to $p_{opt}$, we maximize the fast-switching growth average. 

The finite switching critical period ($T_{crit}$) roughly estimates the maximum time you can spend without the positive phases failing to compensate for the negative phases. 
* By setting $p$ to $p_{opt}$, we maximize the "slack" in the system's average growth. 
* A maximized $y^*$ and $\lambda_{fast}$ typically means the system is furthest from the extinction boundary, granting maximum tolerance for slower switching. 
* Therefore, the optimal duty cycle $p_{opt}$ that maximizes $y^*$ generally correlates with maximizing the allowable switching period (lowering the critical switching frequency).