## Floquet Analysis for Linear Interaction Models

Following the [demonstration of switching-driven persistence in the exponential efficacy-loss model](./floquet_analysis.md), a natural question arises: does this paradoxical survival rely on the specific non-linear exponential form ($e^{-\alpha x}$), or can simpler, linear modifications to the Lotka-Volterra dynamics produce the same phenomenon?

To answer this, we investigate two alternative models with linear parameters: alternating the predator's intrinsic death rate ($\gamma_E$), and alternating the predator's attack rate ($\delta_E$).

### Model 1: The Death-Rate Switching Model

We begin by varying the intrinsic death rate of the predator, $\gamma_E$, while keeping the standard Lotka-Volterra linear interaction. The system is defined as:

$$
\begin{align*}
    \frac{dx}{dt} &= x\left(1-\frac{x}{K_E}\right) - x y \\
    \frac{dy}{dt} &= x y - \gamma_E y = y(x - \gamma_E)
\end{align*}
$$

With equal-time switching between environments $E \in \{A, B\}$:
*   **Environment A:** $K_A = 1, \gamma_A = 1$
*   **Environment B:** $K_B = 2, \gamma_B > 0$

#### The Floquet (Invasion) Exponent

When the predator density is negligible ($y \approx 0$), the prey dynamics decouple to $dx/dt = x(1 - x/K_E)$. The Floquet exponent for predator invasion is the time-averaged relative growth rate:
$$
\lambda = \lim_{T_{total} \to \infty} \frac{1}{T_{total}} \int_0^{T_{total}} \left[ x_0(t) - \gamma(t) \right] dt
$$

#### Slow Switching Limit ($T \to \infty$)

Under slow switching, the prey spends the vast majority of its time at the environmental carrying capacities ($K_A$ and $K_B$). The invasion exponent is the simple average of the local growth rates:
$$
\lambda_{slow} = \frac{1}{2} (g_A + g_B)
$$
Where:
*   $g_A = K_A - \gamma_A = 1 - 1 = 0$
*   $g_B = K_B - \gamma_B = 2 - \gamma_B$

Because $g_A = 0$, the predator is linearly neutrally stable in Environment A isolated, but strictly goes extinct due to non-linear effects (any $y > 0$ decreases $x < 1$, pushing $dy/dt < 0$).
For the predator to fail to invade under the slow switching limit, we require $\lambda_{slow} < 0$:
$$
\frac{1}{2} (0 + 2 - \gamma_B) < 0 \implies 1 - \frac{\gamma_B}{2} < 0 \implies \gamma_B > 2
$$

#### Fast Switching Limit ($T \to 0$)

Under extreme high-frequency switching, the prey experiences the *average* logistic vector field:
$$
\frac{dx}{dt} \approx \frac{1}{2}x(1-x) + \frac{1}{2}x\left(1-\frac{x}{2}\right) = x - \frac{3}{4}x^2
$$
The effective carrying capacity for the fast-switching prey is $\bar{x} = \frac{4}{3} \approx 1.333$.
The invasion exponent evaluates this averaged prey density against the alternating death rates:
$$
\lambda_{fast} = \bar{x} - \frac{\gamma_A + \gamma_B}{2} = \frac{4}{3} - \frac{1 + \gamma_B}{2} = \frac{5}{6} - \frac{\gamma_B}{2}
$$
For the predator to successfully invade in the fast-switching environment, we need $\lambda_{fast} > 0$:
$$
\frac{5}{6} - \frac{\gamma_B}{2} > 0 \implies \gamma_B < \frac{5}{3} \approx 1.667
$$

#### Conclusion on Model 1

To achieve switching-driven persistence, the predator must die out in the slow limit ($\gamma_B > 2$) but survive in the fast limit ($\gamma_B < 1.667$). These two conditions are mutually exclusive. It is mathematically impossible to find a parameter $\gamma_B$ that permits a survival window.

---

### Model 2: The Attack-Rate Switching Model

Instead of varying the death rate, we fix $\gamma = 1$ in all environments and allow the interaction strength (attack rate) $\delta_E$ to fluctuate.

$$
\begin{align*}
    \frac{dx}{dt} &= x\left(1-\frac{x}{K_E}\right) - \delta_E x y \\
    \frac{dy}{dt} &= \delta_E x y - \gamma y = y(\delta_E x - \gamma)
\end{align*}
$$

With equal-time switching between environments $E \in \{A, B\}$:
*   **Environment A:** $K_A = 1, \delta_A = 1$
*   **Environment B:** $K_B = 2, \delta_B > 0$

#### Slow Switching Limit ($T \to \infty$)

The local growth rates are:
*   $g_A = \delta_A K_A - \gamma = 1(1) - 1 = 0$
*   $g_B = \delta_B K_B - \gamma = 2\delta_B - 1$

For the predator to fail to invade under the slow switching limit, we require $\lambda_{slow} < 0$:
$$
\lambda_{slow} = \frac{1}{2} (g_A + g_B) = \frac{1}{2} (2\delta_B - 1) = \delta_B - 0.5 < 0 \implies \delta_B < 0.5
$$

#### Fast Switching Limit ($T \to 0$)

The averaged prey density is unaffected by the predator on the $y \approx 0$ invasion manifold, remaining at $\bar{x} = 4/3$. 
The invasion exponent evaluates time spent in both environments with this constant prey density:
$$
\lambda_{fast} = \frac{1}{2}(\delta_A \bar{x} - \gamma) + \frac{1}{2}(\delta_B \bar{x} - \gamma) = \bar{x} \frac{\delta_A + \delta_B}{2} - \gamma
$$
$$
\lambda_{fast} = \frac{4}{3} \left( \frac{1 + \delta_B}{2} \right) - 1 = \frac{2}{3}(1 + \delta_B) - 1 = \frac{2}{3}\delta_B - \frac{1}{3}
$$
For the predator to successfully invade in the fast-switching environment, we need $\lambda_{fast} > 0$:
$$
\frac{2}{3}\delta_B - \frac{1}{3} > 0 \implies 2\delta_B - 1 > 0 \implies \delta_B > 0.5
$$

#### Conclusion on Model 2

Once again, the required conditions for paradox persistence point in entirely opposite directions. The slow-limit extinction bound ($\delta_B < 0.5$) and the fast-limit survival bound ($\delta_B > 0.5$) collide exactly at $\delta_B = 0.5$, leaving zero overlap for a survival window. 

---

### The Averaged Prey Density and the Necessity of Non-Linearity

Both linear implementations—altering either the baseline death rate or the interaction strength—completely failed to produce switching-driven persistence. We can prove that *any* linear parameter addition to this alternating system is guaranteed to fail by analyzing the temporal averaging of the decoupled prey density, $x_0(t)$.

For equal-time switching, we can apply the variable substitution $u(t) = 1/x_0(t)$. The decoupled prey equation $dx/dt = x(1 - x/K_E(t))$ transforms exactly into a linear differential equation:
$$
\frac{du}{dt} = -u + \frac{1}{K_E(t)}
$$

Because this differential equation is strictly linear in $u$, the time-averaged value of $u(t)$ over one full period, $\langle u \rangle$, evaluates exclusively to the simple average of the fluctuating forcing term, **entirely independent of the switching frequency**:
$$
\langle u \rangle = \frac{1}{2} \left(\frac{1}{K_A} + \frac{1}{K_B}\right) = \frac{1}{2} \left(1 + \frac{1}{2}\right) = \frac{3}{4}
$$

However, the predator's linear invasion fitness depends on the time-averaged actual prey density, $\langle x_0(t) \rangle$, which mathematically equals $\langle 1/u(t) \rangle$.

Because the function $f(u) = 1/u$ is strictly convex for $u > 0$ (its second derivative is positive), Jensen's Inequality dictates that the average of the function is strictly greater than the function of the average: $\langle 1/u \rangle \ge 1/\langle u \rangle$.

*   **When switching is slow ($T \to \infty$):** $u(t)$ traces out large temporal extremes, maximizing its variance. This maximizes the convexity advantage, resulting in the highest possible average effective prey density: $\langle x_0(t) \rangle_{slow} = \frac{1 + 2}{2} = 1.5$.
*   **When switching is fast ($T \to 0$):** $u(t)$ barely reacts to the changing environments and collapses toward a constant state where its temporal variance is zero. Without variance to exploit the convex curve, the average collapses identically to the function of the average: $\langle x_0(t) \rangle_{fast} \to 1/\langle u \rangle = \frac{4}{3} \approx 1.333$.

The profound consequence is that **as environmental switching happens faster, the effective average prey density strictly drops**. 

If the predator's per-capita growth rate is a purely linear function of the prey density (e.g., $dy/dt = y(c_1(E) x - c_2(E))$), the predator will always face a mathematically *worse* environment during high-frequency switching than during slow-frequency switching. If the predator cannot survive the slow limit, it is guaranteed by Jensen's inequality to die out even harder in the fast limit.

Therefore, switching-driven persistence strictly necessitates a **non-linear mechanism**—such as the $e^{-\alpha_E x}$ term in the exponential efficacy-loss model. The non-linearity must allow the beneficial effects of rapid temporal averaging to dynamically override and rescue the predator from the intrinsic foundational penalty on time-averaging prey density.
