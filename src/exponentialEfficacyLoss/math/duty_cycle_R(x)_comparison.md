
# Analytical Threshold for $y^* \leq 0$ Under Fast Switching

## Setup

The shielding model with $\gamma_A = \gamma_B = \gamma$, $\alpha_A = 0$:

$$\frac{dx}{dt} = r_E x\left(1 - \frac{x}{K_E}\right) - xye^{-\alpha_E x}, \qquad \frac{dy}{dt} = xye^{-\alpha_E x} - \gamma y$$

With environments $E \in \{A, B\}$ and $\alpha_A = 0$.

The shielding model averaged over a duty cycle $p$ (fraction of time in env A):

$$\frac{dx}{dt} = p\!\left[x\!\left(1-\frac{x}{K_A}\right) - xy\right] + (1-p)\!\left[x\!\left(1-\frac{x}{K_B}\right) - xye^{-\alpha_B x}\right]$$

$$\frac{dy}{dt} = p\!\left[xy - \gamma y\right] + (1-p)\!\left[xye^{-\alpha_B x} - \gamma y\right]$$

with $\alpha_A = 0$, $\gamma_A = \gamma_B = \gamma$, $r_A = r_B = 1$.

## Solving for the Coexistence Fixed Point

**From $dy/dt = 0$** (dividing by $y > 0$), define $h_A(x) = x - \gamma$ and $h_B(x) = xe^{-\alpha_B x} - \gamma$:

$$p\,h_A(x^*) + (1-p)\,h_B(x^*) = 0 \tag{1}$$

Solving for $p$ explicitly:

$$p = \frac{\gamma - x^* e^{-\alpha_B x^*}}{x^*(1 - e^{-\alpha_B x^*})} \tag{2}$$

For $p \in (0,1)$ to be valid, $h_A$ and $h_B$ must have opposite signs at $x^*$, which requires $x^* \in (K_A, K_B)$ (where $h_A > 0$ and $h_B < 0$).

**From $dx/dt = 0$** (dividing by $x > 0$):

$$y^* = \frac{p(1 - x^*/K_A) + (1-p)(1 - x^*/K_B)}{p + (1-p)e^{-\alpha_B x^*}} \tag{3}$$

## Simplifying $y^*$

Writing $s = e^{-\alpha_B x^*}$ and using $1 - p = (x^* - \gamma)/[x^*(1-s)]$ from equation (2), the denominator of $y^*$ simplifies cleanly:

$$p + (1-p)s = \frac{\gamma - x^*s}{x^*(1-s)} + \frac{(x^* - \gamma)s}{x^*(1-s)} = \frac{\gamma(1-s)}{x^*(1-s)} = \frac{\gamma}{x^*}$$

Substituting into equation (3):

$$y^* = \frac{x^*}{\gamma} \cdot \frac{(\gamma - x^*s)(1 - x^*/K_A) + (x^* - \gamma)(1 - x^*/K_B)}{x^*(1-s)}$$

$$\boxed{y^* = \frac{(\gamma - x^*s)\!\left(1 - \dfrac{x^*}{K_A}\right) + (x^* - \gamma)\!\left(1 - \dfrac{x^*}{K_B}\right)}{\gamma(1 - s)}} \tag{4}$$

The denominator $\gamma(1-s) > 0$ always (since $s = e^{-\alpha_B x^*} < 1$ for $\alpha_B, x^* > 0$). So the sign of $y^*$ is entirely determined by the numerator.

## Condition $y^* = 0$

Setting the numerator of (4) to zero:

$$(\gamma - x^*s)\!\left(1 - \frac{x^*}{K_A}\right) + (x^* - \gamma)\!\left(1 - \frac{x^*}{K_B}\right) = 0$$

Now substitute $K_A = \gamma = 1$, so $1 - x^*/K_A = 1 - x^*$:

$$(1 - x^*s)(1 - x^*) + (x^* - 1)(1 - x^*/K_B) = 0$$

Factor out $(x^* - 1)$, noting $(1 - x^*) = -(x^* - 1)$:

$$(x^* - 1)\!\left[(1 - x^*/K_B) - (1 - x^*s)\right] = 0$$

This gives two cases:

1. $x^* = 1 = K_A$ — the degenerate boundary case where the fixed point has collapsed to $K_A$.
2. $1 - x^*/K_B = 1 - x^*s$, i.e.:

$$\frac{x^*}{K_B} = x^* s = x^* e^{-\alpha_B x^*} \implies e^{-\alpha_B x^*} = \frac{1}{K_B}$$

$$\implies \alpha_B x^* = \ln K_B \implies x^* = \frac{\ln K_B}{\alpha_B} \tag{5}$$

## The Critical $\alpha_B^*$

Equation (5) gives the $x^*$ at which the coexistence fixed point has $y^* = 0$. As $\alpha_B$ increases, $x^* = \ln K_B / \alpha_B$ decreases. The fixed point ceases to be feasible — i.e. $y^* \leq 0$ for **all** $p$ — when $x^*$ exits the interval $(K_A, K_B)$ through the left boundary $K_A$:

$$x^* = K_A \implies \frac{\ln K_B}{\alpha_B^*} = K_A \implies \boxed{\alpha_B^* = \frac{\ln K_B}{K_A}}$$

With $K_A = 1,\ K_B = 2,\ \gamma = 1$:

$$\alpha_B^* = \ln 2 \approx 0.693$$

For $\alpha_B < \alpha_B^*$: the curve $x^* = \ln K_B / \alpha_B > K_A$ lies inside $(K_A, K_B)$, meaning there exists a $p \in (0,1)$ giving $y^* > 0$.

For $\alpha_B \geq \alpha_B^*$: $x^* \leq K_A$, the fixed point is outside the feasible region, and $y^* \leq 0$ for all $p$.

## Agreement with the $R(x)$ Criterion

The switching rescue function is:

$$R(x) = \frac{h_B(x)}{f_B(x)} - \frac{h_A(x)}{f_A(x)}$$

where $f_E(x) = r_E x(1 - x/K_E)$. On the interval $(K_A, K_B)$: $f_B > 0$ and $f_A < 0$, so $R(x) > 0$ requires $h_A$ and $h_B$ to have opposite signs — exactly the same condition needed for fast switching to rescue the predator.

The zero of $R(x)$ exiting the interval through $K_A$ is therefore the threshold condition.

Setting $R(K_A) = 0$. Since $\gamma = K_A = 1$, we have $h_A(K_A) = K_A - \gamma = 0$ and $f_A(K_A) = 0$ simultaneously, so we evaluate the limit via L'Hôpital:

$$\frac{h_A(x)}{f_A(x)} = \frac{x - \gamma}{x(1 - x/K_A)} \xrightarrow{x \to K_A} \frac{1}{-1} = -1$$

Thus:

$$R(K_A) = 0 \implies \frac{h_B(K_A)}{f_B(K_A)} = -1$$

$$\frac{K_A e^{-\alpha_B K_A} - \gamma}{K_A\!\left(1 - K_A/K_B\right)} = -1$$

Substituting $K_A = 1,\ K_B = 2,\ \gamma = 1$:

$$\frac{e^{-\alpha_B} - 1}{1 \cdot \tfrac{1}{2}} = -1 \implies e^{-\alpha_B} - 1 = -\tfrac{1}{2} \implies e^{-\alpha_B} = \tfrac{1}{2}$$

$$\boxed{\alpha_B^* = \ln 2 \approx 0.693}$$


## Summary

| Quantity | Value |
|---|---|
| $\alpha_B^*$ (threshold) | $\ln(K_B / K_A) = \ln 2 \approx 0.693$ |
| $x^*$ at threshold | $K_A = 1$ |
| $y^*$ at threshold | $0$ |
| $p^*$ at threshold | $1$ (degenerate) |

For any $\alpha_B < \ln 2$, there exists a duty cycle $p \in (0,1)$ such that the averaged system has a coexistence fixed point with $y^* > 0$. At $\alpha_B = \ln 2$, this fixed point merges with the boundary at $x^* = K_A$, $y^* = 0$, and for $\alpha_B > \ln 2$ no such fixed point exists.
