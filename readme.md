## Survival only by switching

Let's say we have two species, $x$ and $y$, living in an ecosystem. The ecosystem goes through two different environments, $A$ and $B$.

Is it possible to construct a situation where $y$ dies out in environment $A$ and also in environment $B$, but can persist when we switch between the environments?

To be frankly honest: the idea is not my own, but that of my supervisor [Kunihiko Kaneko](https://scholar.google.com/citations?user=X9ffPWMAAAAJ). However, I was drawn to it, and spent a while trying to create something that works.

It's easy to do it with 3-species, but 2 is a bit more interesting. In addition, wanting to simplify things, I want to have as few parameters and rules as possible: while still being biologically realistic!

As a first step, the equations that work are

$$
\begin{align*}    
    \frac{\mathrm{d}x}{\mathrm{d}t} &= x\left(1-\frac{x}{k_E}\right) - x y e^{-\alpha_E x} \\
    \frac{\mathrm{d}y}{\mathrm{d}t} &= x y e^{-\alpha_E x} - y
\end{align*}
$$

with parameters:
$$
\begin{align*}
k_A &= 1 &&& k_B &= 2 \\
\alpha_A &= 0 &&& \alpha_B &= 0.35
\end{align*}
$$
