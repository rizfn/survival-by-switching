# Fast-Switching Threshold for the Base LV Model

## Setup

The base Lotka-Volterra model with environment-dependent parameters $r_E, K_E, \gamma_E$:

$$\frac{dx}{dt} = r_E x\!\left(1 - \frac{x}{K_E}\right) - xy, \qquad \frac{dy}{dt} = xy - \gamma_E y$$

We fix $r_1 = K_1 = \gamma_1 = 1$ and $r_2 = 4$, $K_2 = 2$, and treat $\gamma_2$ as the free parameter.

## Coexistence Fixed Point of the Averaged System

With duty cycle $p$ (fraction of time in env 1), the averaged system is:

$$\frac{dx}{dt} = p\!\left[r_1 x\!\left(1-\frac{x}{K_1}\right) - xy\right] + (1-p)\!\left[r_2 x\!\left(1-\frac{x}{K_2}\right) - xy\right] = 0$$

$$\frac{dy}{dt} = p\!\left[xy - \gamma_1 y\right] + (1-p)\!\left[xy - \gamma_2 y\right] = 0$$

**From $dy/dt = 0$** (dividing by $y > 0$), with $h_E(x) = x - \gamma_E$:

$$p(x^* - \gamma_1) + (1-p)(x^* - \gamma_2) = 0$$

Since this is linear in $x^*$, it solves exactly:

$$\boxed{x^* = p\gamma_1 + (1-p)\gamma_2} \tag{1}$$

The coexistence fixed point sits at the weighted average of the two predator nullclines.

**From $dx/dt = 0$** (dividing by $x > 0$):

$$p\!\left[r_1\!\left(1 - \frac{x^*}{K_1}\right) - y^*\right] + (1-p)\!\left[r_2\!\left(1 - \frac{x^*}{K_2}\right) - y^*\right] = 0$$

Solving for $y^*$:

$$\boxed{y^* = p\, r_1\!\left(1 - \frac{x^*}{K_1}\right) + (1-p)\, r_2\!\left(1 - \frac{x^*}{K_2}\right)} \tag{2}$$

## Simplifying $y^*$

Substitute $K_1 = \gamma_1 = r_1 = 1$, $r_2 = 4$, $K_2 = 2$, and $x^* = p + (1-p)\gamma_2$ into equation (2).

First, note that since $\gamma_1 = K_1 = 1$:

$$1 - \frac{x^*}{K_1} = 1 - x^* = 1 - p - (1-p)\gamma_2 = (1-p)(1-\gamma_2)$$

And:

$$1 - \frac{x^*}{K_2} = 1 - \frac{p + (1-p)\gamma_2}{2} = \frac{2 - p - (1-p)\gamma_2}{2}$$

Substituting:

$$y^* = p(1-p)(1-\gamma_2) + (1-p) \cdot 4 \cdot \frac{2 - p - (1-p)\gamma_2}{2}$$

$$= (1-p)\!\left[p(1-\gamma_2) + 2\!\left(2 - p - (1-p)\gamma_2\right)\right]$$

$$= (1-p)\!\left[p - p\gamma_2 + 4 - 2p - 2(1-p)\gamma_2\right]$$

$$= (1-p)\!\left[4 - p - \gamma_2\!\left(p + 2(1-p)\right)\right]$$

$$\boxed{y^* = (1-p)\!\left[4 - p - \gamma_2(2-p)\right]} \tag{3}$$

Since $(1-p) > 0$ for $p \in [0,1)$, the sign of $y^*$ is determined entirely by the bracket.

## Condition $y^* = 0$

Setting the bracket to zero:

$$4 - p - \gamma_2(2-p) = 0 \implies \gamma_2 = \frac{4-p}{2-p} \tag{4}$$

So for a given $\gamma_2$, the coexistence fixed point has $y^* = 0$ at duty cycle $p^0 = \frac{2\gamma_2 - 4}{\gamma_2 - 1}$ (inverting equation 4). The fixed point has $y^* > 0$ iff $\gamma_2 < (4-p)/(2-p)$.

## The Critical $\gamma_2^*$

For $y^* > 0$ to be achievable by **some** duty cycle, we need:

$$\gamma_2 < \sup_{p \,\in\, (p_{\min},\, 1)} \frac{4-p}{2-p}$$

The function $(4-p)/(2-p)$ is strictly increasing in $p$ (its derivative is $(6-2p)/(2-p)^2 > 0$), so the supremum is approached as $p \to 1^-$:

$$\sup_p \frac{4-p}{2-p} = \lim_{p\to 1} \frac{4-p}{2-p} = \frac{3}{1} = 3$$

Therefore:

- For $\gamma_2 < 3$: there exists a $p \in (0,1)$ such that $y^* > 0$ — fast switching can sustain the predator.
- For $\gamma_2 \geq 3$: $(4-p)/(2-p) \leq 3 \leq \gamma_2$ for all valid $p$ — no duty cycle gives $y^* > 0$.

$$\boxed{\gamma_2^* = 3}$$

## Independent Check via $R(x)$

The switching rescue function for the base LV model is:

$$R(x) = \frac{h_2(x)}{f_2(x)} - \frac{h_1(x)}{f_1(x)}$$

with $h_E(x) = x - \gamma_E$ and $f_E(x) = r_E x(1 - x/K_E)$.

**$R(x)$ zero analytically.** On $(K_1, K_2) = (1,2)$: $f_1 < 0$ and $f_2 > 0$, and $h_1(x)/f_1(x)$ has a removable singularity at $K_1 = \gamma_1 = 1$ resolved by L'Hôpital:

$$\frac{h_1(x)}{f_1(x)} = \frac{x-1}{x(1-x)} \xrightarrow{x\to 1} \frac{1}{-1} = -1 \implies -\frac{h_1}{f_1} \to +1$$

For the interior, with $r_1=K_1=1$ and $r_2=4$, $K_2=2$:

$$\frac{h_1}{f_1} = \frac{x-1}{x(1-x)} = \frac{-(1-x)}{x(1-x)} = -\frac{1}{x}$$

$$\frac{h_2}{f_2} = \frac{x - \gamma_2}{4x(1-x/2)} = \frac{x-\gamma_2}{2x(2-x)}$$

Setting $R(x) = 0$:

$$\frac{x-\gamma_2}{2x(2-x)} = -\frac{1}{x} \implies \frac{x-\gamma_2}{2(2-x)} = -1 \implies x - \gamma_2 = -2(2-x) = 2x-4$$

$$\boxed{x^*_{R=0} = 4 - \gamma_2} \tag{5}$$

**The zero of $R(x)$ exits $(K_1, K_2)$ through $K_1 = 1$ when:**

$$4 - \gamma_2 = 1 \implies \gamma_2^* = 3$$

This is exactly the threshold found from the fixed-point calculation.

## Summary

| $\gamma_2$ | $R(x)=0$ at $x$ | $R>0$ region | $\max_p\, y^*$ |
|---|---|---|---|
| $2.0$ | $x = 2.0 = K_2$ | empty (boundary) | $0$ |
| $2.5$ | $x = 1.5$ | $(1,\ 1.5)$ | $> 0$ |
| $2.8$ | $x = 1.2$ | $(1,\ 1.2)$ | $> 0$ |
| $3.0$ | $x = 1.0 = K_1$ | empty (boundary) | $0$ |
| $3.5$ | $x = 0.5 < K_1$ | empty | $< 0$ |

The two approaches agree exactly:

$$\gamma_2^* = 3 = 4 - K_1 = r_2 K_2 - (r_2 - r_1) K_1$$

which in general parameters is the solution to $R(K_1) = 0$, i.e. the zero of $R$ touching $K_1$ from the right.

## General Formula

For arbitrary $r_1, K_1, \gamma_1, r_2, K_2$, the same calculation gives:

$$x^* = p\gamma_1 + (1-p)\gamma_2, \quad y^* = (1-p)\!\left[\gamma_2\!\left(\frac{r_1}{K_1} - \frac{r_2}{K_2}\right)\cdot\frac{x^*}{\gamma_2-\gamma_1} + \cdots\right]$$

and the threshold from $R(K_1) = 0$ is:

$$\frac{h_2(K_1)}{f_2(K_1)} = \frac{h_1(K_1)}{f_1(K_1)} \implies \frac{K_1 - \gamma_2}{r_2 K_1(1-K_1/K_2)} = -\frac{1}{r_1} \implies \gamma_2^* = K_1 + \frac{r_2 K_1(1-K_1/K_2)}{r_1}$$

Substituting $r_1=K_1=\gamma_1=1$, $r_2=4$, $K_2=2$: $\gamma_2^* = 1 + 4\cdot 1\cdot(1-1/2)/1 = 1 + 2 = 3$. ✓