# Kelvin-Helmholtz Instability (2D)

Reference:
> A.S. Chamarthi, S.H. Frankel, A. Chintagunta, Implicit gradients based novel finite volume scheme for compressible single and multi-component flows, arXiv preprint arXiv:2106.01738 (2021).

**Key Parameters**

* $\gamma = 5/3$

**Initial Conditions**

* **Spatial Domain:** $(x, y) \in [0, 1] \times [0, 1]$

$$
p = 2.5, \rho(x, y) = \begin{cases}
2, \text{if } 0.25 < y \leq 0.75\\
1, \text{else, }\\
\end{cases}
$$

$$
u(x, y) = \begin{cases}
0.5, \text{if } 0.25 < y \leq 0.75\\
-0.5, \text{else, }
\end{cases}
$$

$$
v(x, y) = 0.1 \sin(4\pi x) \left[ \exp\left( -\frac{(y-0.75)^2}{2\sigma^2} \right) + \exp\left( -\frac{(y-0.25)^2}{2\sigma^2} \right) \right]
$$

where $\sigma = \frac{0.05}{\sqrt{2}}$.

### Initial State 
<img src="figure0.png" height="MAX_HEIGHT"/>

### Evolved State  
<img src="figure1.png" height="MAX_HEIGHT"/>
