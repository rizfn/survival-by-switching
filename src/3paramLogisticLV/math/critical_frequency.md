# Critical Switching Period $T_c$ for Predator Persistence

The right way to define the critical switching period is through the **invasion exponent** of the predator along the extinction boundary $y \to 0$. In that limit the prey decouples, so the sign of the long-term predator growth rate determines whether switching is slow enough to cause extinction or fast enough to permit persistence.

If the switching is periodic with dwell time $T/2$ in each environment, then the full period is $T$.

The critical period $T_c$ is the value, if it exists, where the predator invasion exponent changes sign:

$$
\lambda(T_c) = 0.
$$

---

## 1. Predator growth on the extinction boundary

Near $y=0$, your predator equation has the form

$$
\frac{dy}{dt} = y\,(x-\gamma_E),
$$

in environment $E \in \{A,B\}$, so the per-capita predator growth rate is simply $x-\gamma_E$.

Therefore the predator growth over one full switching cycle is

$$
\Delta \ln y
= \int_0^T \bigl(x_0(t)-\gamma(t)\bigr)\,dt,
$$

where $x_0(t)$ is the prey trajectory on the $y=0$ boundary and $\gamma(t)$ is the corresponding switching value, $\gamma_A$ or $\gamma_B$.

The corresponding invasion exponent is the cycle average:

$$
\lambda(T) = \frac{1}{T} \int_0^T \bigl(x_0(t)-\gamma(t)\bigr)\,dt.
$$

Persistence occurs when $\lambda(T) > 0$; extinction occurs when $\lambda(T) < 0$.

---

## 2. Exact dependence on the initial condition

If switching starts in environment $A$ at prey density $x_0$, then on the extinction boundary the prey obeys the logistic equations

$$
\frac{dx}{dt} = r_E x\left(1 - \frac{x}{K_E}\right), \qquad E \in \{A,B\}.
$$

To evaluate the invasion exponent, we first write this prey flow in closed form. Separation of variables gives

$$
\int \frac{dx}{x(1-x/K_E)} = \int r_E\,dt.
$$

Using partial fractions,

$$
\ln\!\left(\frac{x}{K_E-x}\right) = r_E t + C.
$$

Imposing the initial condition $x(0)=x_{\mathrm{in}}$ gives the exact flow in a fixed environment $E$:

$$
x_E(t; x_{\mathrm{in}})
=
\frac{K_E x_{\mathrm{in}} e^{r_E t}}{K_E + x_{\mathrm{in}}(e^{r_E t} - 1)}.
$$

If the system begins in $A$, define the prey density at the switch by

$$
x_1 = x_A\left(\frac{T}{2}; x_0\right),
$$

and then evolve through $B$ from $x_1$.

The per-cycle predator gain is the sum of the two dwell contributions,

$$
\Delta \ln y
=
\int_0^{T/2} \bigl(x_A(t; x_0) - \gamma_A\bigr)\,dt
+
\int_0^{T/2} \bigl(x_B(t; x_1) - \gamma_B\bigr)\,dt.
$$

Using the logistic flow, each integral has a closed form:

$$
\int_0^{T/2} x_E(t; x_{\mathrm{in}})\,dt
=
\frac{K_E T}{2} - \frac{K_E}{r_E}\ln\!\left(\frac{x_E\left(\frac{T}{2}; x_{\mathrm{in}}\right)}{x_{\mathrm{in}}}\right).
$$

Equivalently, the predator itself satisfies

$$
\frac{1}{y}\frac{dy}{dt} = x(t)-\gamma_E,
$$

so over one dwell interval

$$
 \ln\!\frac{y\left(\frac{T}{2}\right)}{y(0)} = \int_0^{T/2} \bigl(x_E(t;x_{\mathrm{in}})-\gamma_E\bigr)\,dt.
$$

Combining these pieces, the exact exponent over one full cycle is

$$
\lambda(T; x_0)
=
\frac{1}{T}
\left[
\sum_{E \in \{A,B\}} (K_E - \gamma_E)\frac{T}{2}
-
\sum_{E \in \{A,B\}} \frac{K_E}{r_E}
\ln\!\left(\frac{x_E\left(\frac{T}{2}; x_{\mathrm{in},E}\right)}{x_{\mathrm{in},E}}\right)
\right].
$$

where $x_{\mathrm{in},A} = x_0$ and $x_{\mathrm{in},B} = x_1$.

The critical period is therefore defined implicitly by

$$
\lambda(T_c; x_0) = 0.
$$

This is the exact answer for a prescribed initial condition and initial switching phase.

---

## 3. What depends on the initial condition?

There are two different notions of criticality:

1. **Finite-time threshold from a chosen initial state.**
	Here $T_c$ really does depend on $x_0$ and on whether the cycle starts in $A$ or $B$, because the first few dwell intervals affect the prey trajectory entering the exponent.

2. **Asymptotic Floquet threshold on the periodic orbit.**
	After transients, the prey converges to the periodic boundary orbit $x_*(t;T)$ on $y=0$. Then the invasion exponent becomes
	$$
	\lambda(T)
	= \frac{1}{T}\int_0^T \bigl(x_*(t;T) - \gamma(t)\bigr)\,dt,
	$$
	which depends on the parameters and the switching rule, but not on the initial condition.

So if you want the mathematically clean switching threshold, $T_c$ is a function of the environment parameters and the switching protocol. If you want the threshold experienced by a particular initial condition, keep $x_0$ explicit and solve the implicit equation above.

---

## 4. Slow and fast limits

The limit values bracket the answer:

$$
\lambda_{slow} = \frac{1}{2}\Big[(K_A - \gamma_A) + (K_B - \gamma_B)\Big],
$$

and, from the averaged fast-switching prey equilibrium,

$$
\bar{x}_{fast} = \frac{r_A + r_B}{\frac{r_A}{K_A} + \frac{r_B}{K_B}},
\qquad
\lambda_{fast} = \bar{x}_{fast} - \frac{\gamma_A + \gamma_B}{2}.
$$

If

$$
\lambda_{slow} < 0 < \lambda_{fast},
$$

then continuity implies at least one critical period $T_c \in (0,\infty)$.

---

## 5. Practical computation of $T_c$

There is generally no universal closed-form formula for $T_c$. The practical route is:

1. Choose the initial prey density $x_0$ and the initial environment.
2. Compute the boundary trajectory $x_0(t)$ for one switching cycle.
3. Evaluate $\lambda(T)$.
4. Solve $\lambda(T)=0$ numerically.

That gives the exact critical switching period for the chosen initial condition and parameter set.

---

## 6. Compact transcendental form

Starting from the previous equation,

$$
\lambda(T_c;x_0)=0,
$$

we multiply by $T_c$ and substitute the two half-cycle flows:

$$
0
=
\sum_{E\in\{A,B\}} (K_E-\gamma_E)\frac{T_c}{2}
-
\sum_{E\in\{A,B\}} \frac{K_E}{r_E}
\ln\!\left(\frac{x_E\left(\frac{T_c}{2};x_{\mathrm{in},E}\right)}{x_{\mathrm{in},E}}\right).
$$

For the $A$-dwell, write

$$
\frac{x_A\left(\frac{T_c}{2};x_0\right)}{x_0}
=
\frac{K_A e^{r_A T_c/2}}{K_A+x_0\bigl(e^{r_A T_c/2}-1\bigr)}
=
\frac{e^{r_A T_c/2}}{1+\frac{x_0}{K_A}\bigl(e^{r_A T_c/2}-1\bigr)}.
$$

Taking logs, the $r_A T_c/2$ term cancels against the linear contribution from the first sum, leaving

$$
\frac{K_A}{r_A}\ln\!\left(1+\frac{x_0}{K_A}\bigl(e^{r_A T_c/2}-1\bigr)\right).
$$

Doing the same for the $B$-dwell gives the analogous term with $x_1(T_c)=x_A(T_c/2;x_0)$. So the condition becomes the single scalar root problem

$$
F(T)
=
\frac{K_A}{r_A}\ln\!\Bigl(1+\frac{x_0}{K_A}\bigl(e^{r_A T/2}-1\bigr)\Bigr)
+
\frac{K_B}{r_B}\ln\!\Bigl(1+\frac{x_1(T)}{K_B}\bigl(e^{r_B T/2}-1\bigr)\Bigr)
-
\frac{\gamma_A+\gamma_B}{2}T,
$$

where

$$
x_1(T)=x_A\!\left(\frac{T}{2};x_0\right)
=
\frac{K_A x_0 e^{r_A T/2}}{K_A+x_0\bigl(e^{r_A T/2}-1\bigr)}.
$$

Then the critical period is given by the transcendental equation

$$
F(T_c)=0.
$$

This is usually the cleanest form for computation: it is scalar, explicit, and still captures the full switching dependence through the exponentials and logarithms.
