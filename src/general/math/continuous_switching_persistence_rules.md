# Rigorous Theory of Continuous Environmental Rescue

## 1. Goal and mathematical setup

$$
\Lambda[\sigma] \approx \frac{1}{T}\int_0^T \lambda_{\mathrm{fr}}(\sigma(t))\,dt.
$$
We want conditions under which:

1. for every constant environment $\sigma(t)\equiv\sigma_0$, the predator goes extinct;
2. for some nonconstant periodic or chaotic $\sigma(t)$, the predator persists.
$$
	\text{persistence} \iff \text{top exponent } > 0.
$$
### Assumption 1 (Invariant extinction boundary)

The boundary $y=0$ is invariant for every frozen environment:
$$
G(x,0,\sigma)=0 \quad \text{for all } x\ge 0,\ \sigma\in\Sigma.
$$
Hence the predator-free prey equation is
$$
\frac{dx}{dt}=f(x,\sigma(t)),\qquad f(x,\sigma):=F(x,0,\sigma).
$$

### Assumption 2 (Frozen prey equilibrium)

For each constant $\sigma\in\Sigma$, the prey-only equation
$$
\dot x=f(x,\sigma)
$$
has a unique globally attracting positive equilibrium $x=K(\sigma)$.

### Assumption 3 (Boundary predator growth rate)

Define the predator per-capita growth rate on the boundary by
$$
h(x,\sigma):=\partial_y G(x,0,\sigma).
$$
Then near $y=0$ the predator equation is exactly
$$
\dot y = h(x(t),\sigma(t))\,y.
$$

### Assumption 4 (Frozen extinction)

For every constant environment $\sigma\in\Sigma$,
$$
\lambda_{\mathrm{fr}}(\sigma):=h(K(\sigma),\sigma)<0.
$$
This means every frozen environment is a sink for the predator.

The problem is to decide when a time-dependent $\sigma(t)$ can reverse this and produce persistence.

## 2. Constant $\sigma$: extinction is necessary and sufficient

If $\sigma(t)\equiv\sigma_0$ is constant, then the prey converges to $K(\sigma_0)$ and the predator linearization is
$$
\dot y = h(K(\sigma_0),\sigma_0)\,y.
$$

Therefore:

- if $\lambda_{\mathrm{fr}}(\sigma_0)<0$, then $y(t)\to 0$ exponentially;
- if $\lambda_{\mathrm{fr}}(\sigma_0)>0$, then the predator invades;
- if $\lambda_{\mathrm{fr}}(\sigma_0)=0$, the frozen system is critical.

So the frozen-environment extinction condition is simply
$$
\lambda_{\mathrm{fr}}(\sigma)<0 \quad \text{for all } \sigma\in\Sigma.
$$

This is necessary and sufficient for extinction in every constant environment.

## 3. Periodic $\sigma(t)$: exact invasion criterion

Now assume $\sigma(t)$ is $T$-periodic and $C^1$.

Under the same dissipativity assumptions on the prey equation, the boundary system admits a positive $T$-periodic solution $x_*(t)$ solving
$$
\dot x = f(x,\sigma(t)),\qquad x_*(t+T)=x_*(t).
$$

Along this orbit, the predator satisfies
$$
\dot y = h(x_*(t),\sigma(t))\,y.
$$

Its Floquet exponent is
$$
\Lambda[\sigma] := \frac{1}{T}\int_0^T h(x_*(t),\sigma(t))\,dt.
$$

### Theorem 1 (Periodic persistence criterion)

For a periodic forcing $\sigma(t)$, the predator persists near the boundary if and only if
$$
\Lambda[\sigma] > 0.
$$
It goes extinct if and only if
$$
\Lambda[\sigma] < 0.
$$

Thus the exact necessary and sufficient condition for rescue by a continuous periodic environment is positivity of the boundary Floquet exponent.

### Corollary 1 (Extinction under uniform negativity)

If there exists $\eta>0$ such that
$$
h(x_*(t),\sigma(t))\le -\eta\qquad\text{for all }t,
$$
then $\Lambda[\sigma]\le -\eta<0$ and extinction is unavoidable.

### Corollary 2 (Sufficient condition for rescue)

If
$$
\int_0^T h(x_*(t),\sigma(t))\,dt>0,
$$
then the periodic forcing rescues the predator.

## 4. What makes rescue possible?

To have the phenomenon “dies for every constant $\sigma$ but survives for a varying $\sigma(t)$,” two statements must hold simultaneously:

1. **Frozen extinction:**
	$$
	\lambda_{\mathrm{fr}}(\sigma)<0\qquad\text{for all constant }\sigma\in\Sigma.
	$$

2. **Forced invasion:** there exists a nonconstant periodic $\sigma(t)$ such that
	$$
	\Lambda[\sigma]>0.
	$$

This is the clean necessary-and-sufficient answer: constant environments are sinks, but a time-dependent environment can rescue the predator exactly when the boundary Floquet exponent becomes positive.

### Local perturbation formula

If
$$
\sigma(t)=\sigma_0+\varepsilon\eta(t),
$$
with $\eta$ $T$-periodic and $\varepsilon$ small, and
$$
x_*(t)=K(\sigma_0)+\varepsilon\xi(t)+O(\varepsilon^2),
$$
then
$$
\Lambda[\sigma]=\lambda_{\mathrm{fr}}(\sigma_0)+\varepsilon\Lambda_1+O(\varepsilon^2),
$$
where $\Lambda_1$ is obtained by differentiating $h(x,\sigma)$ and the linearized prey equation. In particular, rescue is possible for small forcing only if the first-order correction $\Lambda_1$ is large enough to overcome the negative frozen exponent.

## 5. Fast and slow forcing limits

These limits are useful because they give explicit sufficient conditions in terms of $F$ and $G$.

### Slow forcing

If $\sigma(t)$ varies slowly, the prey tracks the moving frozen equilibrium:
$$
x_*(t)\approx K(\sigma(t)).
$$
Then
$$
\Lambda[\sigma] \approx \frac{1}{T}\int_0^T \lambda_{\mathrm{fr}}(\sigma(t))\,dt.
$$
So in the adiabatic regime, rescue requires the time-average of the frozen invasion exponent to be positive.

### Fast forcing

If $\sigma(t)$ varies rapidly, the prey sees the averaged vector field
$$
\bar f(x):=\frac{1}{T}\int_0^T f(x,\sigma(t))\,dt,
$$
with averaged equilibrium $\bar K$ solving $\bar f(\bar K)=0$.
Then
$$
\Lambda[\sigma]\approx \frac{1}{T}\int_0^T h(\bar K,\sigma(t))\,dt.
$$

So fast forcing rescues the predator exactly when the predator's averaged growth along the averaged prey state is positive.

## 6. Bonus: chaotic $\sigma(t)$

Now assume $\sigma(t)$ is chaotic rather than periodic.

### 6.1 Ergodic chaotic forcing

Suppose $\sigma(t)$ is generated by an ergodic chaotic flow, so that the coupled boundary dynamics admit an invariant measure. Then the predator's long-time exponent is the almost-sure Lyapunov exponent
$$
\Lambda_{\mathrm{chaos}}:=\lim_{T\to\infty}\frac{1}{T}\int_0^T h(x(t),\sigma(t))\,dt,
$$
whenever this limit exists.

By Birkhoff’s ergodic theorem, this equals the expectation of $h$ with respect to the invariant measure of the forced boundary system.

### Theorem 2 (Chaotic persistence criterion)

For ergodic chaotic forcing, the predator persists if and only if
$$
\Lambda_{\mathrm{chaos}}>0,
$$
and it goes extinct if and only if
$$
\Lambda_{\mathrm{chaos}}<0.
$$

Hence chaotic forcing can rescue the predator, but only if the invariant measure spends enough time in predator-favorable states to make the long-time average positive.

### 6.2 Non-ergodic chaotic forcing

If the chaotic forcing is not ergodic, there is no universal average. One must compute the top Lyapunov exponent on the specific attractor realized by the forcing. The criterion is unchanged in form:
$$
	\text{persistence} \iff \text{top exponent } > 0.
$$

## 7. Summary

The rigorous answer is:

1. **Frozen environments:** extinction occurs for every constant $\sigma$ if and only if
    $$
        \lambda_{\mathrm{fr}}(\sigma)=h(K(\sigma),\sigma)<0\quad\text{for all }\sigma
    $$


2. **Periodic forcing:** rescue occurs if and only if the boundary Floquet exponent is positive:
    $$
    	\Lambda[\sigma]=\frac{1}{T}\int_0^T h(x_*(t),\sigma(t))\,dt>0.
	$$

3. **Chaotic forcing:** rescue occurs if and only if the boundary Lyapunov exponent is positive:
    $$
    	\Lambda_{\mathrm{chaos}}>0.
	$$

So the phenomenon “dies for every constant environment but survives under varying $\sigma(t)$” is completely governed by the sign of the long-time boundary growth exponent induced by the forcing.
