# Floquet Analysis for the Growth-Rate Switching Model

We have established that linear modifications to the death rate ($\gamma$) or interaction strength ($\delta$) strictly fail to produce switching-driven persistence because fast switching inherently drops the time-averaged prey density.

However, could alternating the prey's intrinsic **growth rate ($r_E$)** between environments bypass this penalty, allowing for survival in a purely linear Lotka-Volterra interaction model? 

The answer is **yes**.

## 1. The Growth-Rate Switching System

Consider the entirely linear interaction model ($\alpha = 0$), where the environments distinctively alter both the carrying capacity ($K_E$) and the prey's intrinsic growth rate ($r_E$).

$$
\begin{align*}
    \frac{dx}{dt} &= r_E x\left(1-\frac{x}{K_E}\right) - x y \\
    \frac{dy}{dt} &= x y - \gamma_E y = y(x - \gamma_E)
\end{align*}
$$

With equal-time switching between environments $E \in \{A, B\}$:
*   **Environment A:** $K_A = 1, \gamma_A = 1, r_A = 1$
*   **Environment B:** $K_B = 2, \gamma_B > 0, r_B > 0$

## 2. Linearization at $y = 0$

On the extinction manifold ($y \approx 0$), the prey dynamics decouple into the switching logistic equation:
$$
\frac{dx_0}{dt} = r(t) x_0\left(1 - \frac{x_0}{K(t)}\right)
$$

The Floquet exponent for predator invasion remains the time-averaged relative growth rate:
$$
\lambda = \lim_{T \to \infty} \frac{1}{T} \int_0^{T} \left[ x_0(t) - \gamma(t) \right] dt
$$

## 3. Slow and Fast Limits

### Slow Switching Limit ($T \to \infty$)
Under slow switching, the prey strictly tracks the local carrying capacities $K_A$ and $K_B$. The growth rate $r_E$ dictates how quickly it reaches equilibrium, but at $T \to \infty$, the percentage of time spent *in transit* approaches zero. 
Thus, $x_0(t)$ essentially alternates perfectly between $x=1$ and $x=2$.

The invasion exponent is the average of the isolated growths:
$$
\lambda_{slow} = \frac{1}{2} (g_A + g_B) = \frac{1}{2} \left[ (K_A - \gamma_A) + (K_B - \gamma_B) \right]
$$
$$
\lambda_{slow} = \frac{1}{2}(0 + 2 - \gamma_B) = 1 - \frac{\gamma_B}{2}
$$

For the predator to be neutrally stable in isolated Env A (and therefore perish due to non-linear feedbacks), and to perish strictly in isolated Env B, we require $\gamma_B > K_B = 2$.
Thus, if **$\gamma_B > 2$**, both environments act as sinks ($\lambda_{slow} < 0$).

### Fast Switching Limit ($T \to 0$)
Under extreme rapid switching, the prey reacts to the *averaged* vector field rather than tracking the individual equilibriums. Because the factor $r_E$ multiplies the entire vector field, it directly weights the averaging process.

$$
\frac{dx}{dt} \approx \frac{1}{2} \left[ r_A x \left(1 - \frac{x}{K_A}\right) + r_B x \left(1 - \frac{x}{K_B}\right) \right]
$$

Solving for the effective rapid-switching equilibrium $\bar{x}$ where $dx/dt = 0$:
$$
r_A \left(1 - \bar{x}\right) + r_B \left(1 - \frac{\bar{x}}{2}\right) = 0
$$
$$
r_A + r_B = \bar{x} \left(r_A + \frac{r_B}{2}\right) \implies \bar{x} = \frac{r_A + r_B}{r_A + \frac{r_B}{2}} = 2 \frac{r_A + r_B}{2r_A + r_B}
$$

The Floquet exponent evaluates this effective density against the average death rate:
$$
\lambda_{fast} = \bar{x} - \frac{\gamma_A + \gamma_B}{2} = 2 \frac{r_A + r_B}{2r_A + r_B} - \frac{1 + \gamma_B}{2}
$$

## 4. The Switching-Exclusive Survival Window

Can we achieve $\lambda_{slow} < 0$ and $\lambda_{fast} > 0$ simultaneously by tuning our free parameters $r_B$ and $\gamma_B$?

Yes. If we set the growth rate in Environment B to be highly explosive, for instance **$r_B = 10$**:
$$
\bar{x} = 2 \frac{1 + 10}{2 + 10} = 2 \frac{11}{12} = \frac{11}{6} \approx 1.833
$$

Now calculate the survival exponent boundary:
$$
\lambda_{fast} > 0 \implies \frac{11}{6} - \frac{1 + \gamma_B}{2} > 0 \implies \frac{11}{3} > 1 + \gamma_B \implies \gamma_B < \frac{8}{3} \approx 2.667
$$

Combining our conditions, if we select any $\gamma_B$ such that:
$$ 2 < \gamma_B < 2.667 $$

(For example, **$\gamma_B = 2.2$**), we orchestrate the paradox exactly:
1. **Dies out** in Environment A isolated ($K_A - \gamma_A = 0 \rightarrow$ non-linear extinction).
2. **Dies out** in Environment B isolated ($K_B - \gamma_B = 2 - 2.2 = -0.2 < 0$).
3. **Dies out** under slow switching ($\lambda_{slow} = 1 - 1.1 = -0.1 < 0$).
4. **Persists and invades** under fast switching ($\lambda_{fast} = 1.833 - \frac{3.2}{2} = 0.233 > 0$).

## 5. Why Growth-rate Switching Bypasses the Averaging Penalty

In the earlier death-rate analysis, we noted via Jensen's Inequality that high-frequency switching strictly depressed the effective prey density ($\bar{x} = 1.333$) relative to the slow-frequency average ($1.5$), locking out any linear persistence window. 

Why did switching $r_E$ break that mathematical wall?

If we re-apply the substitution $u(t) = 1/x(t)$ to the new decoupled equation:
$$
\frac{du}{dt} = -r(t) u + \frac{r(t)}{K(t)}
$$

The vital difference is that **$r(t)$ now acts as a time-varying relaxation rate multiplying $u$**. The differential equation is linear, but its decay constant is oscillating. 

This means the temporal average $\langle u \rangle$ is *no longer simply the arithmetic average of the forcing terms*. By assigning a massive dynamical pull to Environment B ($r_B = 10$) compared to Environment A ($r_A = 1$), the high-frequency averaged vector field is wildly skewed toward the high-resource state ($K_B = 2$). 

This disproportionate, asymmetric mathematical "weighting" dynamically drags the fast-averaging prey density $\bar{x}$ upwards to $1.833$, decisively crushing the slow-averaging baseline of $1.5$. By ensuring there are vastly *more* prey available in the fast limit than the slow limit, linear survival parameters easily fit into the resulting surplus.

## 6. Can We Reduce This to Two Parameters?

The model above successfully orchestrates the paradox by fluctuating **three** parameters across the environments: the carrying capacity ($K_E$), the growth rate ($r_E$), and the death rate ($\gamma_E$). 

A natural question is: can we keep one of these constant and drive the paradox with only **two** fluctuating parameters? 

Mathematically, **no**. It is demonstrably impossible algebraically. To establish the paradox, we require three strict criteria:
1. Both isolated environments must be sinks.
2. The slow-switching average must be a sink ($\lambda_{slow} < 0$).
3. The fast-switching average must be a source ($\lambda_{fast} > 0$).

Let's test keeping each parameter constant:

### Scenario 1: Constant Growth Rate ($r_A = r_B$)
As proven in the previous section via Jensen's Inequality, if $r$ is constant, the fast-switching prey average is strictly penalized compared to the slow-switching average ($\bar{x}_{fast} < \langle x \rangle_{slow}$). If the predator dies in the slow limit, it will die even harder in the fast limit. Thus, **$r$ must switch** to artificially weight $\bar{x}_{fast}$ upwards.

### Scenario 2: Constant Carrying Capacity ($K_A = K_B = K$)
If the carrying capacity is fixed, the prey density never fluctuates; it rests permanently at $x(t) = K$. The environment offers no variation in resources. Because the prey density is completely constant, temporal averaging does nothing: $\lambda_{slow} = \lambda_{fast}$. If the predator dies at slow frequencies, it dies at fast frequencies. Thus, **$K$ must switch** to create rest-state variance in the first place.

### Scenario 3: Constant Death Rate ($\gamma_A = \gamma_B = \gamma$)
To satisfy criterion #1 (isolated sinks), the constant death rate must be high enough to kill the predator even in the most resource-rich environment. This requires $\gamma \ge \max(K_A, K_B)$.

However, the fast-averaging resource density, $\bar{x}$, is a weighted average between $K_A$ and $K_B$. By mathematical definition, any average is strictly less than the maximum of its components: $\bar{x} < \max(K_A, K_B)$. 

Therefore, the fast-switching Floquet exponent is inherently capped:
$$ \lambda_{fast} = \bar{x} - \gamma < \max(K_A, K_B) - \gamma \le 0 $$
The predator can never accumulate a high enough averaged prey density to overcome a unified death rate that suppresses the richest environment. Thus, **$\gamma$ must switch** to dodge this mathematical floor.

### Conclusion
In a completely linear interaction framework, **exactly three parameters must oscillate**. You need $K$ to create resource variance, $r$ to break the penalty of Jensen's inequality (dynamically weighting the fast-average upwards), and $\gamma$ to dodge the death rate cap so that the environments independently still act as sinks.

## 7. Derivation of the General Parametric Boundaries

We can securely derive the exact algebraic boundaries required for persistence using only the mathematical definitions of the six parameters ($K_A, K_B, r_A, r_B, \gamma_A, \gamma_B$).

Assume without loss of generality that Environment B is the more resource-rich environment, so $K_B > K_A$.

To achieve the paradox, we need:
1.  **Isolated Sinks:** $\gamma_A \ge K_A$ and $\gamma_B > K_B$. 
    *(Note: This strictly guarantees that the slow-switching average is also a sink, since $\frac{\gamma_A + \gamma_B}{2} > \frac{K_A + K_B}{2}$.)*
2.  **Fast-Switching Source:** $\lambda_{fast} > 0 \implies \bar{x}_{fast} > \frac{\gamma_A + \gamma_B}{2}$.

Combining these constraints bounds the averaged death rate, forming the baseline **survival window**:
$$ \frac{K_A + K_B}{2} < \frac{\gamma_A + \gamma_B}{2} < \bar{x}_{fast} $$

### The Universal Existence Constraint

For this survival window to exist at all, the upper bound must strictly exceed the lower bound:
$$ \frac{K_A + K_B}{2} < \bar{x}_{fast} $$

Substitute the harmonic form of $\bar{x}_{fast}$:
$$ \frac{K_A + K_B}{2} < \frac{r_A + r_B}{\frac{r_A}{K_A} + \frac{r_B}{K_B}} $$

Multiply the numerator and denominator of the RHS by $K_A K_B$:
$$ \frac{K_A + K_B}{2} < \frac{K_A K_B (r_A + r_B)}{r_A K_B + r_B K_A} $$

Cross-multiply to remove fractions (we can safely do this since all environmental parameters are strictly positive):
$$ (K_A + K_B)(r_A K_B + r_B K_A) < 2 K_A K_B (r_A + r_B) $$

Expand both sides algebraically:
$$ r_A K_A K_B + r_B K_A^2 + r_A K_B^2 + r_B K_A K_B < 2 r_A K_A K_B + 2 r_B K_A K_B $$

Subtract the RHS from the LHS to set the inequality to zero:
$$ r_B K_A^2 - r_B K_A K_B + r_A K_B^2 - r_A K_A K_B < 0 $$

Group the terms by rate parameters ($r_A$ and $r_B$) to factor them:
$$ r_B K_A (K_A - K_B) + r_A K_B (K_B - K_A) < 0 $$

Factor out a $-1$ from the first term to match $(K_B - K_A)$:
$$ -r_B K_A (K_B - K_A) + r_A K_B (K_B - K_A) < 0 $$

Pull out the common $(K_B - K_A)$ multiplier:
$$ (r_A K_B - r_B K_A)(K_B - K_A) < 0 $$

Because we defined Environment B as the resource-rich state, we know definitively that $(K_B - K_A) > 0$. Therefore, to satisfy the $< 0$ inequality, the other binomial must be strictly negative:
$$ r_A K_B - r_B K_A < 0 \implies r_A K_B < r_B K_A $$

Finally, dividing both sides by $K_A K_B$ produces the fundamental law:
$$ \frac{r_A}{K_A} < \frac{r_B}{K_B} $$

### Necessary and Sufficient Conditions

To guarantee switching-driven persistence in this purely linear model (where the predator dies in both isolated environments and under slow switching, but thrives under fast switching), the parameter set must satisfy a strictly necessary and sufficient system of inequalities.

Assuming without loss of generality that Environment B has the higher carrying capacity ($K_B > K_A$):

1. **Isolated Extinctions:** 
   $$ \gamma_A \ge K_A $$
   $$ \gamma_B \ge K_B $$
   *(Note: At least one of these must be a strict inequality, $>$, to ensure the slow-switching average is strictly negative).*

2. **Fast-Switching Survival (The Window):**
   $$ \frac{\gamma_A + \gamma_B}{2} < \frac{r_A + r_B}{\frac{r_A}{K_A} + \frac{r_B}{K_B}} $$

These are the only algebraic rules you need to check. If you satisfy 1 and 2, the paradox is mathematically guaranteed. 

However, structurally, a parameter set can *only* satisfy rule 2 without violating rule 1 if the underlying parameters also obey two structural laws:
* **The Asymmetric Growth Law:** $\frac{r_A}{K_A} < \frac{r_B}{K_B}$ (Otherwise the right side of rule 2 is mathematically smaller than $\frac{K_A + K_B}{2}$, violating rule 1).
* **The Gamma Oscillation Law:** $\gamma_A < \gamma_B$ (Otherwise a flat $\gamma$ must be $\ge K_B$, which strictly exceeds the right side of rule 2).

Let's verify the working setup from Section 4 against the necessary and sufficient conditions:
* **Parameters:** $K_A=1, K_B=2, r_A=1, r_B=10, \gamma_A=1, \gamma_B=2.2$
* Rule 1 (Extinctions): $1 \ge 1$ and $2.2 > 2$ 
* Rule 2 (Invasion): 
  $$ \frac{1 + 2.2}{2} < \frac{1 + 10}{1/1 + 10/2} $$
  $$ 1.6 < 1.833 $$