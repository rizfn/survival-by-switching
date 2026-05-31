# Fixed Point Analysis of the Averaged Vector Field (Fast Switching Limit)

Persistence in the fast-switching limit mathematically implies that the **time-averaged vector field** gives birth to a stable, attracting interior fixed point where strictly $y > 0$.

Let's rigorously analyze the fixed points of this averaged field to see how they directly connect back to your $\lambda_{fast} > 0$ Floquet condition.

## 1. The Averaged Vector Field Equations

In the infinite-frequency limit (interval $\to 0$), the populations react to the algebraic average of the two distinct Lotka-Volterra vector fields.

For environments A and B, the averaged system is:

$$
\begin{align}
\frac{dx}{dt} &= \frac{1}{2} \left[ r_A x \left(1 - \frac{x}{K_A}\right) - x y \right] + \frac{1}{2} \left[ r_B x \left(1 - \frac{x}{K_B}\right) - x y \right] \\
\frac{dy}{dt} &= \frac{1}{2} (x y - \gamma_A y) + \frac{1}{2} (x y - \gamma_B y)
\end{align}
$$

We can simplify these by grouping the predator and prey terms:

$$
\begin{align}
\frac{dx}{dt} &= \frac{x}{2} \left[ r_A \left(1 - \frac{x}{K_A}\right) + r_B \left(1 - \frac{x}{K_B}\right) \right] - x y \\
\frac{dy}{dt} &= y \left( x - \frac{\gamma_A + \gamma_B}{2} \right)
\end{align}
$$

## 2. Solving for the Fixed Points 

To find the fixed points (equilibriums), we set $\frac{dx}{dt} = 0$ and $\frac{dy}{dt} = 0$.

From the predator equation (Eq. 2), we have two possible states:
1. $y = 0$ (Predator Extinction)
2. $x = \frac{\gamma_A + \gamma_B}{2}$ (Predator Survival/Coexistence)

### State 1: The Extinction Manifold ($y = 0$)
If $y = 0$, we substitute this into Eq. 1 and cancel out an $x$ (ignoring the trivial origin $0,0$):
$$
0 = \frac{1}{2} \left[ r_A \left(1 - \frac{x}{K_A}\right) + r_B \left(1 - \frac{x}{K_B}\right) \right]
$$
This is precisely the equation we solved in the Floquet analysis for the fast-switching effective prey density.
$$
\bar{x}_{fast} = \frac{r_A + r_B}{\frac{r_A}{K_A} + \frac{r_B}{K_B}}
$$
This gives us the **Prey-Only Fixed Point**: $(x = \bar{x}_{fast}, y = 0)$.

### State 2: Coexistence / Persistence Fixed Point ($y > 0$)
For a true coexistence fixed point to exist, the density of the prey is forced by the predator equation to equal the average death rate:
$$
x^* = \bar{\gamma} = \frac{\gamma_A + \gamma_B}{2}
$$

To find the resting predator density $y^*$, substitute $x = \bar{\gamma}$ back into Eq. 1 and set it to 0:
$$
0 = \frac{\bar{\gamma}}{2} \left[ r_A \left(1 - \frac{\bar{\gamma}}{K_A}\right) + r_B \left(1 - \frac{\bar{\gamma}}{K_B}\right) \right] - \bar{\gamma} y^*
$$

Assuming $\bar{\gamma} > 0$, we divide out $\bar{\gamma}$ and isolate $y^*$:
$$
y^* = \frac{1}{2} \left[ r_A \left(1 - \frac{\bar{\gamma}}{K_A}\right) + r_B \left(1 - \frac{\bar{\gamma}}{K_B}\right) \right]
$$

## 3. The Condition for $y^* > 0$

For the system to persist, this fixed point must be physically meaningful; we need strictly $y^* > 0$.

When does this happen? Let's check the inequality:
$$
y^* > 0 \implies \frac{1}{2} \left[ (r_A + r_B) - \bar{\gamma} \left( \frac{r_A}{K_A} + \frac{r_B}{K_B} \right) \right] > 0
$$

Which simplifies to:
$$
r_A + r_B > \bar{\gamma} \left( \frac{r_A}{K_A} + \frac{r_B}{K_B} \right)
$$

Substitute back the definition of $\bar{x}_{fast}$:
$$
\bar{x}_{fast} = \frac{r_A + r_B}{\frac{r_A}{K_A} + \frac{r_B}{K_B}} > \bar{\gamma}
$$

This perfectly echoes back to the Floquet exponent! In Floquet theory, $\lambda_{fast} = \bar{x}_{fast} - \bar{\gamma}$. 

Thus: **The interior fixed point lifts vertically off the $y=0$ extinction axis ($y^* > 0$) if and only if $\lambda_{fast} > 0$.** 

## 4. Applying Your Canvas Vector Parameters

Let's plug in the exact parameters from your Javascript simulation canvas:
* **Env 1:** $K_A = 1, r_A = 1, \gamma_A = 1$
* **Env 2:** $K_B = 1.5, r_B = 2.3, \gamma_B = 1.6$

The average death rate sets the stable prey coordinate:
$$
x^* = \frac{1 + 1.6}{2} = 1.3
$$

Now calculate the resulting predator density $y^*$ at this $x$:
$$
y^* = \frac{1}{2} \left[ 1 \left(1 - \frac{1.3}{1}\right) + 2.3 \left(1 - \frac{1.3}{1.5}\right) \right]
$$
$$
y^* = \frac{1}{2} \left[ (1 - 1.3) + 2.3 \left(1 - 0.8666...\right) \right]
$$
$$
y^* = \frac{1}{2} \left[ -0.3 + 2.3 \left(\frac{2}{15}\right) \right]
$$
$$
y^* = \frac{1}{2} \left[ -0.3 + \frac{4.6}{15} \right] = \frac{1}{2} \left[ -0.3 + 0.30666... \right]
$$
$$
y^* = \frac{1}{2} \left[ \frac{1}{150} \right] = \frac{1}{300} \approx 0.00333
$$

The averaged paradox yields an attracting sink exactly at **$(1.3, \frac{1}{300})$**.

Because $y^* > 0$, the averaged geometry inherently generates a stable coexistence point, keeping the particles permanently trapped parallel to (but strictly above) extinction!