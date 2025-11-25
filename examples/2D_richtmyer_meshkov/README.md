# Richtmyer-Meshkov Instability (2D)

Reference:
> A.S. Chamarthi, S.H. Frankel, A. Chintagunta, Implicit gradients based novel finite volume scheme for compressible single and multi-component flows, arXiv preprint arXiv:2106.01738 (2021)., see Example 4.18

**Key Parameters**

* $\lambda = 1.0$
* $Re = 1e4$

**Initial Conditions**

* **Spatial Domain:** $(x, y) \in [0, 16\lambda] \times [0, \lambda]$

The domain is divided into three regions: Post-shock Air, Pre-shock Air, and SF6. The shock position is given by $x_{shock}=0.7\lambda$. The interface position $x_{int}(y)$ is defined by the sinusoidal perturbation:

$$
x_{int}(y) = \lambda \left[ 0.4 - 0.1 \sin\left( 2\pi \left( \frac{y}{\lambda} + 0.25 \right) \right) \right]
$$

The initial state is given by:

$$
(\rho, u, v, p, \gamma) =
\begin{cases}
   (1.4112, 0.8787, 0, 1.6272/1.4, 1.4) & \text{if } x > x_{shock} \quad \text{(Post-shock Air)} \\
   (1.0, 1.24, 0, 1/1.4, 1.4) & \text{if } x_{int}(y) < x \leq x_{shock} \quad \text{(Pre-shock Air)} \\
   (5.04, 1.24, 0, 1/1.4, 1.093) & \text{if } x \leq x_{int}(y) \quad \text{(SF6)} \\
\end{cases}
$$

The initial perturbation is then smoothed by 
$$f_{sm} = \frac{1}{2}\left(1 + \text{erf}\left(\frac{\Delta D}{E_i\sqrt{\Delta x \Delta y}}\right)\right)$$

$$u = u_L(1 - f_{sm}) + u_Rf_{sm}$$

where $u$ is each primitive variable. $E_i$ is a thickness constant, and $\Delta D$ is the signed horizontal distance from $x_{int}(y)$. $u_L$ and $u_R$ are left and right interface conditions.


### Initial State
<img src="figure0.png" height="MAX_HEIGHT"/>

### Evolved State
<img src="figure1.png" height="MAX_HEIGHT"/>
