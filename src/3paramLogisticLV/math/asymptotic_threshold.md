# Asymptotic Switching Threshold $T_c$ for Predator Persistence

The asymptotic switching threshold is the value of the switching period where the predator changes from eventual decay to eventual persistence on the periodic boundary orbit. Unlike the finite-time threshold, it does not depend on a chosen initial prey density.

The clean definition is through the predator's Floquet exponent on the extinction boundary $y \to 0$.

---

## 1. Start from the model

Near the predator-free boundary, the predator equation has the form

$$
\frac{dy}{dt} = y\,(x-\gamma_E),
$$

where $E \in \{A,B\}$ denotes the current environment.

The prey equation on the boundary is logistic in each environment:

$$
\frac{dx}{dt} = r_E x\left(1-\frac{x}{K_E}\right).
$$

Because $y$ is small, the prey dynamics determine the predator's long-time growth rate.

---

## 2. The prey flow on the boundary

In a fixed environment $E$, the logistic equation is separable:

$$
\int \frac{dx}{x(1-x/K_E)} = \int r_E\,dt.
$$

Using partial fractions,

$$
\ln\!\left(\frac{x}{K_E-x}\right) = r_E t + C.
$$

So the exact boundary flow is

$$
 x_E(t; x_{\mathrm{in}})
=
\frac{K_E x_{\mathrm{in}} e^{r_E t}}{K_E + x_{\mathrm{in}}(e^{r_E t} - 1)}.
$$

For equal dwell times, the environment switches every $T/2$, so one full period is $T$.

---

## 3. Periodic boundary orbit

For the asymptotic threshold, we do not keep a specific initial condition in the final formula. Instead, after the transient has died out, the prey lies on the periodic boundary orbit $x_*(t;T)$, where

$$
x_*(t+T;T)=x_*(t;T).
$$

Let

$$
x_0(T)=x_*(0;T), 
\qquad
x_1(T)=x_*(T/2;T).
$$

Then the two half-period flows satisfy

$$
x_1(T)=x_A\!\left(\frac{T}{2};x_0(T)\right),
$$

$$
x_0(T)=x_B\!\left(\frac{T}{2};x_1(T)\right).
$$

These are the periodicity conditions.

---

## 4. Expand the Floquet exponent

Along the periodic boundary orbit, the predator satisfies

$$
\frac{1}{y}\frac{dy}{dt}=x_*(t;T)-\gamma(t),
$$

where $\gamma(t)=\gamma_A$ on $[0,T/2]$ and $\gamma(t)=\gamma_B$ on $[T/2,T]$.

So over one full cycle,

$$
\lambda(T)=\frac{1}{T}\int_0^T \bigl(x_*(t;T)-\gamma(t)\bigr)\,dt
=\frac{1}{T}\left[\int_0^{T/2}\bigl(x_A(t;x_0(T)) - \gamma_A\bigr)\,dt
+\int_0^{T/2}\bigl(x_B(t;x_1(T)) - \gamma_B\bigr)\,dt\right].
$$

Using the closed form for the logistic flow,

$$
\int_0^{T/2} x_E(t;x_{\mathrm{in}})\,dt
=
\frac{K_E T}{2}-\frac{K_E}{r_E}\ln\!\left(\frac{x_E\left(\frac{T}{2};x_{\mathrm{in}}\right)}{x_{\mathrm{in}}}\right),
$$

we get

$$
\lambda(T)
=\frac{1}{T}\Bigg[
\sum_{E\in\{A,B\}}(K_E-\gamma_E)\frac{T}{2}
-\frac{K_A}{r_A}\ln\!\left(\frac{x_1(T)}{x_0(T)}\right)
-\frac{K_B}{r_B}\ln\!\left(\frac{x_0(T)}{x_1(T)}\right)
\Bigg].
$$

Combining the logarithms gives the simplest form:

$$
\lambda(T)
=\frac{K_A+K_B-\gamma_A-\gamma_B}{2}
-\frac{1}{T}\left(\frac{K_A}{r_A}-\frac{K_B}{r_B}\right)
\ln\!\left(\frac{x_1(T)}{x_0(T)}\right).
$$

---

## 5. The transcendental equation for $T_c$

The asymptotic switching threshold is the root of the Floquet exponent:

$$
\lambda(T_c)=0.
$$

Substituting the simplified expression above gives the transcendental equation

$$
\frac{K_A+K_B-\gamma_A-\gamma_B}{2}
=
\frac{1}{T_c}\left(\frac{K_A}{r_A}-\frac{K_B}{r_B}\right)
\ln\!\left(\frac{x_1(T_c)}{x_0(T_c)}\right),
$$

with the periodicity conditions

$$
x_1(T_c)=x_A\!\left(\frac{T_c}{2};x_0(T_c)\right),
\qquad
x_0(T_c)=x_B\!\left(\frac{T_c}{2};x_1(T_c)\right).
$$

So the threshold is obtained by solving this scalar transcendental condition together with the boundary-orbit consistency relations. Starting far from the periodic orbit only changes the transient; it does not change the asymptotic equation itself.
