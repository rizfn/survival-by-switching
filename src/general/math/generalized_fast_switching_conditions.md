# General Constraints for Fast-Switching Persistence

To orchestrate the switching survival paradox in an arbitrary continuous-time predator-prey system, we can abstract away the specific Lotka-Volterra parameters. Let the generalized dynamics in environment $i \in \{1, 2\}$ be defined by:

$$
\begin{align*}
\frac{dx}{dt} &= F_i(x,y) \\
\frac{dy}{dt} &= G_i(x,y) = y \cdot g_i(x,y)
\end{align*}
$$

Here, $g_i(x,y)$ represents the purely per-capita growth rate of the predator. 

To satisfy the paradox completely—that the predator perishes strictly in both isolated environments but thrives toward a stable interior fixed point in the fast-switching limit—any arbitrary functions $F$ and $G$ must satisfy three mathematical pillars.

---

## 1. The Isolated Extinction Condition

In any isolated environment $i$, the predator must be doomed. 

First, the prey must have a globally stable carrying capacity $K_i$ in the absence of predators:
$$ F_i(K_i, 0) = 0 \quad \text{with} \quad \frac{\partial F_i}{\partial x}(K_i, 0) < 0 $$

For the predator to strictly go extinct in this environment, it must not be able to invade when prey are at their maximum equilibrium:
$$ g_i(K_i, 0) < 0 $$

To ensure there are absolutely no "hidden" survival states (no interior fixed points), the nullclines of $F_i$ and $G_i$ must never cross in the positive quadrant:
$$ F_i(x,y) = 0 \text{ and } G_i(x,y) = 0 \text{ have no simultaneous solution for } (x>0, y>0) $$

---

## 2. The Averaged Invasion Condition

In the fast switching limit with duty cycle $p$ (e.g., $p=0.5$ for equal time), the populations react to the averaged vector field:
$$ F_{avg}(x,y) = p F_1(x,y) + (1-p) F_2(x,y) $$
$$ G_{avg}(x,y) = y \cdot \left[ p g_1(x,y) + (1-p) g_2(x,y) \right] $$

Let the **averaged prey-only equilibrium** be $\bar{x}$, defined as the root where $F_{avg}(\bar{x}, 0) = 0$. Note that $K_1 < \bar{x} < K_2$.

For the averaged environment to allow survival, the predator boundary must become a repeller (an invasion source). This directly requires the averaged per-capita growth rate at $\bar{x}$ to be strictly positive:
$$ g_{avg}(\bar{x}, 0) = p g_1(\bar{x}, 0) + (1-p) g_2(\bar{x}, 0) > 0 $$

### The "Paradox Engine" Implication:
Notice that $g_2(K_2, 0) < 0$. Because $\bar{x} < K_2$, and predator growth usually drops as prey decreases, $g_2(\bar{x}, 0)$ will typically be severely negative. 

Therefore, for $g_{avg} > 0$, **$g_1(\bar{x}, 0)$ must be overwhelmingly positive.**
Because $g_1(K_1, 0) < 0$, this implies $\bar{x}$ must be significantly larger than $K_1$. 

This reveals the two universal structural tricks required to engineer the paradox:
1.  **Asymmetric Prey Dominance ($F_2 \gg F_1$):** Environment 2's vector field must completely dominate the average along the x-axis, dragging $\bar{x}$ exceptionally close to $K_2$ and far away from $K_1$.
2.  **Over-Capacity Predator Spurt:** Environment 1's predator growth function ($g_1$) must explode positively when the prey density artificially exceeds $K_1$ ($x > K_1$).

---

## 3. The Interior Stability Condition (Persistence)

If $g_{avg}(\bar{x}, 0) > 0$, the boundary is unstable. Because biological systems naturally bound populations from growing to infinity (assuming $\partial_y g_i < 0$ and $\partial_x F_i < 0$ for large $x, y$), this local boundary instability guarantees the birth of at least one interior coexistence point $(x^*, y^*)$ or limit cycle.

To ensure it rests perfectly at a stable fixed point (like your canvas), we locate the interior root where $F_{avg}(x^*, y^*) = 0$ and $g_{avg}(x^*, y^*) = 0$. 

Its local stability requires the classic Jacobian trace/determinant conditions evaluated at the averaged field:
$$ \text{Det}(J_{avg}) > 0 \implies \left(\frac{\partial F_{avg}}{\partial x}\right) \left(y^* \frac{\partial g_{avg}}{\partial y}\right) - \left(\frac{\partial F_{avg}}{\partial y}\right) \left(y^* \frac{\partial g_{avg}}{\partial x}\right) > 0 $$

$$ \text{Tr}(J_{avg}) < 0 \implies \frac{\partial F_{avg}}{\partial x} + y^* \frac{\partial g_{avg}}{\partial y} < 0 $$

### Summary Checklist for generalized F and G:
Any pair of environments will exhibit this paradox if and only if they satisfy:
1. $g_1(K_1, 0) < 0$ and $g_2(K_2, 0) < 0$ *(Local environments are sinks)*
2. No intersection of local nullclines *(No isolated survival)*
3. $p g_1(\bar{x}, 0) + (1-p) g_2(\bar{x}, 0) > 0$ *(Averaged field boundary is an invasion source)*
4. $\text{Tr}(J_{avg}) < 0$ and $\text{Det}(J_{avg}) > 0$ *(Interior point is a stable node/spiral)*