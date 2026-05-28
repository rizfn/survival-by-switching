## Floquet Theory and Invasion Analysis for Changing Environments

To understand how the predator ($y$) can survive in a fluctuating environment even when it faces extinction in each static environment, we analyze its ability to "invade" from close to zero density ($y \approx 0$) using Floquet theory.

### The Decoupled Prey Dynamics

The system is defined as:
$$
\begin{align*}
    \frac{dx}{dt} &= x\left(1-\frac{x}{K_E}\right) - x y e^{-\alpha_E x} \\
    \frac{dy}{dt} &= x y e^{-\alpha_E x} - y
\end{align*}
$$
With environments $E \in \{A, B\}$:
*   **Environment A:** $K_A = 1, \alpha_A = 0$
*   **Environment B:** $K_B = 2, \alpha_B > 0$

When the predator density is negligible ($y \approx 0$), the prey dynamics decouple and track the standard logistic growth:
$$
\frac{dx}{dt} = x \left(1 - \frac{x}{K(t)}\right)
$$
Let $x_0(t)$ be the stable periodic solution to this equation as it alternates between environments with period $T$.

### The Floquet (Invasion) Exponent

The stability of the extinction state ($y = 0$) is determined by the long-term relative growth rate of a small predator population. This is the transverse Lyapunov exponent (or Floquet exponent):
$$
\lambda = \lim_{T_{total} \to \infty} \frac{1}{T_{total}} \int_0^{T_{total}} \left[ x_0(t) e^{-\alpha(t) x_0(t)} - 1 \right] dt
$$
If $\lambda > 0$, the predator can invade and survive. If $\lambda \le 0$, the predator is driven to extinction.

### Slow Switching Limit ($T \to \infty$)

If the environments switch very slowly, the prey tracks the stable equilibria of each environment: $x_0(t) \approx K_A$ in Env A, and $x_0(t) \approx K_B$ in Env B. Assuming equal time in both environments ($50\%$ duty cycle), the slow limit invasion exponent is the simple average of the static equilibrium growth rates:
$$
\lambda_{slow} = \frac{1}{2} (g_A + g_B)
$$
Where:
*   $g_A = K_A e^{-\alpha_A K_A} - 1 = 1 \cdot e^0 - 1 = 0$
*   $g_B = K_B e^{-\alpha_B K_B} - 1 = 2 e^{-2 \alpha_B} - 1$

*Note on Environment A:* While $g_A = 0$ means the predator is linearly neutrally stable, non-linear terms ensure extinction. Specifically, any $y > 0$ strictly suppresses $x$ below $K_A=1$, yielding a strictly negative growth rate ($dy/dt = y(x-1) < 0$ since $x < 1$). Thus, the predator goes extinct in Env A isolated.

For the predator to fail to invade across the long-term slow switching average, we need $\lambda_{slow} < 0$. Since $g_A = 0$, this requires $g_B < 0$:
$$
2 e^{-2 \alpha_B} - 1 < 0 \implies \alpha_B > \frac{\ln(2)}{2} \approx 0.3465
$$
In other words, a low value of $\alpha_B$ leads to the creation of a coexistence fixed point in environment B. If $\alpha_B > 0.3465$, the predator dies out under the slow switching limit.

### Fast Switching Limit ($T \to 0$)

If the environments switch extremely rapidly, the prey population doesn't have time to react to the fluctuations. Instead, it equilibrates to the *average* logistic vector field:
$$
\frac{dx}{dt} \approx \frac{1}{2}x\left(1 - \frac{x}{K_A}\right) + \frac{1}{2}x\left(1 - \frac{x}{K_B}\right) = x - \frac{3}{4}x^2
$$
Solving $\frac{dx}{dt} = 0$, the effective carrying capacity for the fast-switching prey is:
$$
\bar{x} = \frac{4}{3} \approx 1.333
$$
The invasion exponent in the high-frequency limit utilizes this averaged prey density against the alternating $\alpha$ values:
$$
\lambda_{fast} = \frac{1}{2} \left[ \bar{x} e^{-\alpha_A \bar{x}} - 1 \right] + \frac{1}{2} \left[ \bar{x} e^{-\alpha_B \bar{x}} - 1 \right]
$$
Substituting $\bar{x} = \frac{4}{3}$ and $\alpha_A = 0$:
$$
\lambda_{fast} = \frac{1}{2} \left( \frac{4}{3} - 1 \right) + \frac{1}{2} \left( \frac{4}{3} e^{-\frac{4}{3}\alpha_B} - 1 \right) = \frac{1}{6} + \frac{2}{3} e^{-\frac{4}{3}\alpha_B} - \frac{1}{2} = \frac{2}{3} e^{-\frac{4}{3}\alpha_B} - \frac{1}{3}
$$
For the predator to successfully invade in the fast-switching environment, we need $\lambda_{fast} > 0$:
$$
\frac{2}{3} e^{-\frac{4}{3}\alpha_B} - \frac{1}{3} > 0 \implies e^{-\frac{4}{3}\alpha_B} > \frac{1}{2} \implies \alpha_B < \frac{3}{4} \ln(2) \approx 0.5198
$$

### The Switching-exclusive survival window

Combining the two limits, we find an analytical result. If we choose the parameter $\alpha_B$ such that:
$$
\frac{\ln(2)}{2} < \alpha_B < \frac{3}{4} \ln(2) 
$$
Or numerically:
$$
0.3465 < \alpha_B < 0.5198
$$
Then the predator:
1. **Dies out** in Environment A isolated ($g_A = 0$ causes non-linear extinction).
2. **Dies out** in Environment B isolated ($g_B < 0$).
3. **Dies out** if switching slowly ($\lambda_{slow} < 0$).
4. **Persists and invades** if switching rapidly ($\lambda_{fast} > 0$).

Therefore, any choice of $\alpha_B$ in this permissible range (such as $\alpha_B = 0.35$) resolves the paradox by demonstrating that the nonlinearity of the invasion exponent under the mathematical effects of averaging provides a strict survival window.

### Critical Switching Frequency ($T_c$)

Because $\lambda_{slow} < 0$ and $\lambda_{fast} > 0$ for a valid $\alpha_B$ (e.g., $\alpha_B = 0.35$), it follows by continuity that there is a critical switching period $T_c$ (and corresponding critical frequency $f_c = 1/T_c$). 

*   If $T > T_c$ (frequency is too low), the population spends too much time dwelling in the isolated environments, and the predator is driven to extinction.
*   If $T < T_c$ (frequency is sufficiently fast), the temporal averaging effect dominates, and the predator successfully invades and coexists.

To find this critical switching period $T_c = 2\tau$ (where $\tau$ is the fixed dwell time in each environment), we set the exact Floquet exponent to zero:
$$
\lambda(\tau) = \frac{1}{2\tau} \left( \int_0^\tau \left[ x_A(t) - 1 \right] dt + \int_0^\tau \left[ x_B(t) e^{-\alpha_B x_B(t)} - 1 \right] dt \right) = 0
$$

We can elegantly simplify the integral for Environment A ($K_A=1, \alpha_A=0$). Using the decoupled prey equation $\frac{dx}{dt} = x(1-x) \implies dt = \frac{dx}{x(1-x)}$:
$$
\begin{align*}
\int_0^\tau (x-1) dt &= \int_{x_A(0)}^{x_A(\tau)} \frac{x-1}{x(1-x)} dx \\
&= \int_{x_A(0)}^{x_A(\tau)} -\frac{1}{x} dx \\
&= -\ln \left( \frac{x_A(\tau)}{x_A(0)} \right)
\end{align*}
$$

The integral for Environment B ($K_B=2$), however, incorporates the nonlinear efficacy loss $e^{-\alpha_B x}$. Substituting $dt = \frac{dx}{x(1-x/2)}$ gives:
$$
\int_0^\tau \left[ x(t) e^{-\alpha_B x(t)} - 1 \right] dt = \int_{x_B(0)}^{x_B(\tau)} \frac{x e^{-\alpha_B x} - 1}{x\left(1 - \frac{x}{2}\right)} dx
$$

We can use partial fractions to decompose the integrand:
$$
\begin{align*}
\frac{x e^{-\alpha_B x} - 1}{x\left(1 - \frac{x}{2}\right)} &= 2 \frac{x e^{-\alpha_B x} - 1}{x(2-x)} \\
&= \frac{2e^{-\alpha_B x}}{2-x} - \frac{2}{x(2-x)} \\
&= \frac{2e^{-\alpha_B x}}{2-x} - \left( \frac{1}{x} + \frac{1}{2-x} \right)
\end{align*}
$$

The first term can be integrated using the Exponential Integral function, defined as $\text{Ei}(z) = \int_{-\infty}^z \frac{e^t}{t} dt$. Substituting $u = \alpha_B(2-x)$ gives:
$$
\int \frac{2e^{-\alpha_B x}}{2-x} dx = -2e^{-2\alpha_B} \text{Ei}(\alpha_B(2-x))
$$

Therefore, completing the integral evaluates to:
$$
\begin{align*}
\int \left( \frac{2x e^{-\alpha_B x} - 2}{x(2-x)} \right) dx &= -2e^{-2\alpha_B} \text{Ei}(\alpha_B(2-x)) - \ln(x) + \ln(2-x) \\
&= -2e^{-2\alpha_B} \text{Ei}(\alpha_B(2-x)) + \ln\left(\frac{2-x}{x}\right)
\end{align*}
$$

Setting the exact Floquet exponent $\lambda(\tau) = 0$ requires the sum of the Environment A and Environment B integrals to be zero:
$$
-\ln \left( \frac{x_A(\tau)}{x_A(0)} \right) + \left[ -2e^{-2\alpha_B} \text{Ei}(\alpha_B(2-x)) + \ln\left(\frac{2-x}{x}\right) \right]_{x_B(0)}^{x_B(\tau)} = 0
$$

Thus, for a fixed $\alpha_B$, finding the critical time period $T_c = 2\tau$ requires resolving the stable periodic boundary values for the prey $(x_A(0), x_A(\tau))$ and utilizing a numerical root-finding algorithm to determine the precise $\tau$ that satisfies this equation.
