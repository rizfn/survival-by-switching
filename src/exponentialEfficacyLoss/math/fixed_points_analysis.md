# Fixed Points and Bifurcation Analysis

For a single, unchanging environment, the system is described by the following ordinary differential equations coupling the prey ($x$) and predator ($y$):

$$
\frac{dx}{dt} = r x \left(1 - \frac{x}{K}\right) - x e^{-\alpha x} y
$$

$$
\frac{dy}{dt} = y \left(x e^{-\alpha x} - \gamma\right)
$$

Here, we use the Ricker-style consumption function $F_E(x) = x e^{-\alpha x}$.

## 1. Fixed Points

Fixed points $(x^*, y^*)$ occur where $\frac{dx}{dt} = 0$ and $\frac{dy}{dt} = 0$. 

From the $\frac{dy}{dt} = 0$ equation, we have two possibilities:
1. $y^* = 0$
2. $x^* e^{-\alpha x^*} = \gamma$

### The Trivial and Prey-Only Fixed Points ($y^* = 0$)
If $y^* = 0$, the $x$ equation simplifies to $r x \left(1 - \frac{x}{K}\right) = 0$. This gives us two fixed points:
- **Total Extinction:** $(0, 0)$
- **Prey Carrying Capacity:** $(K, 0)$

### The Coexistence Fixed Points ($y^* > 0$)
For coexistence, the prey population must satisfy the condition where predator consumption equals mortality:
$$
x^* e^{-\alpha x^*} = \gamma
$$
If a solution $x^*$ exists, the corresponding predator equilibrium is found by setting $\frac{dx}{dt} = 0$:
$$
y^* = \frac{r x^* \left(1 - \frac{x^*}{K}\right)}{x^* e^{-\alpha x^*}} = \frac{r x^*}{\gamma} \left(1 - \frac{x^*}{K}\right)
$$
For this fixed point to be biologically sound ($y^* > 0$), we require $x^* < K$.

The function $F_E(x) = x e^{-\alpha x}$ starts at $0$, increases to a maximum at $x_m = \frac{1}{\alpha}$ with a peak value of $F_E\left(\frac{1}{\alpha}\right) = \frac{1}{\alpha e}$, and then decays back to $0$ as $x \to \infty$.

Because of this hump shape, the equation $x^* e^{-\alpha x^*} = \gamma$ can have zero, one, or two solutions depending on the parameters:
- **No solutions:** If $\gamma > \frac{1}{\alpha e}$
- **One solution:** If $\gamma = \frac{1}{\alpha e}$
- **Two solutions:** If $\gamma < \frac{1}{\alpha e}$, let's call them $x_1^*$ and $x_2^*$ where $x_1^* < \frac{1}{\alpha} < x_2^*$.

## 2. Bifurcations

The system's qualitative behavior changes at critical values of $K$ and $\alpha$. Here are the primary bifurcations:

### Saddle-Node Bifurcation (Predator Existence)
The coexistence fixed points are born (or mutually annihilate) when the horizontal line $y = \gamma$ is exactly tangent to the peak of $F_E(x)$. 
This occurs when:
$$
\gamma = \frac{1}{\alpha e} \implies \alpha = \frac{1}{\gamma e}
$$
- For $\alpha > \frac{1}{\gamma e}$, the predator cannot harvest enough energy to counter mortality ($\gamma$), regardless of the prey density. The predator always goes extinct.
- As $\alpha$ decreases through $\frac{1}{\gamma e}$, a saddle-node bifurcation occurs, yielding a stable node and a saddle point for coexistence.

### Transcritical Bifurcation (Predator Invasion)
Even if coexistence fixed points *can* exist algebraically, they must be realistic (i.e., $x^* < K$). A transcritical bifurcation point occurs when a coexistence fixed point collides with the prey-only fixed point $(K,0)$. At this point, $x^* = K$.
We substitute $x^* = K$ into the coexistence condition:
$$
K e^{-\alpha K} = \gamma
$$
This condition defines the boundary of predator invasion. 
- If $K e^{-\alpha K} > \gamma$, a small introduction of predators into the $(K, 0)$ state will grow ($y$ can invade). 
- If $K e^{-\alpha K} < \gamma$, the $(K, 0)$ state is stable to predator invasions.
