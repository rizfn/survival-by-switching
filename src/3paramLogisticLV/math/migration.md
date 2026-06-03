# Can diffusion save the predator when both isolated patches are sinks?

Consider the two-patch prey–predator system with conservative diffusion between patches. In each patch the local dynamics are logistic prey growth with Lotka-type predation; add linear exchange terms with rates $d_x$ (prey) and $d_y$ (predator). The full system is
$$
\begin{aligned}
\dot x_1 &= r_1 x_1\left(1-\frac{x_1}{K_1}\right) - x_1 y_1 + d_x(x_2-x_1),\\
\dot x_2 &= r_2 x_2\left(1-\frac{x_2}{K_2}\right) - x_2 y_2 + d_x(x_1-x_2),\\
\dot y_1 &= x_1 y_1 - \gamma_1 y_1 + d_y(y_2-y_1),\\
\dot y_2 &= x_2 y_2 - \gamma_2 y_2 + d_y(y_1-y_2).
\end{aligned}
$$

Assume both isolated patches are sinks (no diffusion):
$$
K_1<\gamma_1\quad\text{and}\quad K_2<\gamma_2.
$$

Let $(x_1^*,x_2^*)$ be the predator-free equilibrium of the prey subsystem (solve by setting $y_1=y_2=0$). The equilibrium equations are
$$
r_1 x_1^*\Bigl(1-\frac{x_1^*}{K_1}\Bigr) + d_x(x_2^*-x_1^*) = 0,\qquad
r_2 x_2^*\Bigl(1-\frac{x_2^*}{K_2}\Bigr) + d_x(x_1^*-x_2^*) = 0.
$$

Step 1 — bound on prey equilibria. Suppose for contradiction that for some index $i$ we have $x_i^*>K_i$. Then the local logistic term satisfies $r_i x_i^*(1-x_i^*/K_i)<0$, so the equilibrium equation for patch $i$ implies
$$
d_x(x_j^*-x_i^*) = -r_i x_i^*\Bigl(1-\frac{x_i^*}{K_i}\Bigr) > 0,
$$
hence $x_j^*>x_i^*$, contradicting the assumption that $x_i^*$ was a maximal component. Therefore no component can exceed $\max\{K_1,K_2\}$, i.e.
$$
\max\{x_1^*,x_2^*\}\le\max\{K_1,K_2\}.
$$

Step 2 — sign of predator linear growth. Under the hypothesis that both isolated patches are sinks,
$$
K_1<\gamma_1\quad\text{and}\quad K_2<\gamma_2,
$$
we obtain $x_i^*\le\max\{K_1,K_2\}<\gamma_i$ for each $i$. Define the linearised predator growth rates
$$
a_i := x_i^* - \gamma_i < 0,
$$
for $i=1,2$.

Step 3 — total predator decay. The linearised predator system is
$$
\dot y_1 = a_1 y_1 + d_y(y_2-y_1),\qquad \dot y_2 = a_2 y_2 + d_y(y_1-y_2).
$$
Summing gives the evolution of the total predator density $Y=y_1+y_2$:
$$
\dot Y = a_1 y_1 + a_2 y_2 \le \max(a_1,a_2) (y_1+y_2) = \max(a_1,a_2) Y < 0.
$$
Thus $Y$ decays monotonically to zero and extinction is unavoidable.

Answer (final): No. Under conservative diffusion of prey and predators, if both isolated environments are sinks (each would drive $y\to0$ without diffusion), then for any diffusion rates $d_x,d_y\ge0$ the predator-free prey equilibrium satisfies $x_i^*<\gamma_i$, the linearised predator growth rates are negative, and the total predator density decays. Therefore diffusion cannot save $y$ in that case.

## Time-varying symmetric migration

Suppose predator migration remains symmetric but varies in time: the per-capita exchange rate is $d_y(t)$ with the instantaneous flux terms $d_y(t)(y_2-y_1)$ and $d_y(t)(y_1-y_2)$ in the two equations. The linearised predator system becomes
$$
\dot y_1 = a_1 y_1 + d_y(t)(y_2-y_1),\qquad
\dot y_2 = a_2 y_2 + d_y(t)(y_1-y_2).
$$

Summing these two equations again eliminates the migration terms for every $t$:
$$
\dot Y = a_1 y_1 + a_2 y_2.
$$
Hence if the frozen-time linear growth rates satisfy $a_i<0$ for both patches at all times (for example when the prey subsystem remains such that $x_i^*(t)<\gamma_i$), then the same decay argument applies and $Y$ decreases monotonically regardless of the time dependence of $d_y(t)$.

Caveat — when time dependence can matter: if the prey densities (and thus $a_i$) themselves vary in time because of explicit temporal forcing (e.g. time-varying environmental parameters or time-dependent prey diffusion $d_x(t)$), then temporal averaging can change the sign of the long-term invasion exponent (this is the temporal-switching mechanism analysed elsewhere). By itself, a time-varying symmetric predator migration rate that does not alter prey dynamics cannot rescue the predator when the instantaneous $a_i$ remain negative.

Short answer: No — time-varying symmetric migration alone does not help. To obtain rescue you must allow movement to break conservation (e.g., transport prey/subsidy, biased/state-dependent dispersal, or correlated temporal forcing of prey that creates transient positive $a_i$ values).

## Predator taxis toward prey

Now replace simple diffusion of the predator by prey-oriented movement (taxis). A standard two-patch discrete form is
$$
\dot y_1 = x_1 y_1 - \gamma_1 y_1 + \tau\,y_1(x_2-x_1) - \tau\,y_2(x_2-x_1),
$$
$$
\dot y_2 = x_2 y_2 - \gamma_2 y_2 - \tau\,y_1(x_2-x_1) + \tau\,y_2(x_2-x_1),
$$
where $\tau\ge 0$ measures taxis strength and $x_2-x_1$ is the prey gradient. The extra terms represent predators leaving the lower-prey patch and entering the higher-prey patch. Importantly, these taxis terms are still conservative: when we add the two equations they cancel exactly.

Writing $Y=y_1+y_2$ and linearising the prey at the predator-free equilibrium $x_i^*$ gives
$$
\dot Y = (x_1^*-\gamma_1) y_1 + (x_2^* - \gamma_2) y_2.
$$
The taxis terms do not appear in the total balance. Therefore, if both isolated patches are sinks,
$$
x_1^*<\gamma_1\quad\text{and}\quad x_2^*<\gamma_2,
$$
then
$$
\dot Y \le \max\{x_1^*-\gamma_1,\,x_2^*-\gamma_2\}\,Y < 0.
$$
So prey-oriented taxis alone cannot rescue the predator when both isolated patches are sinks.

Why this is still useful: taxis can change the spatial distribution $y_1:y_2$ and can therefore matter when one patch is a source and the other is a sink, or when prey dynamics are also altered by movement. But under the present hypothesis (both isolated patches sinks), taxis only redistributes predator mass between two losing patches; it does not create a positive total growth rate.

Thus the mathematical answer is still no: any conservative movement law for the predator, even prey-following taxis, leaves the summed growth bounded above by a negative quantity when both local invasion rates are negative.
