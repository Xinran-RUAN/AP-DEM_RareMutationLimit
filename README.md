# AP-DEM_RareMutationLimit
This is the joint work with Weijie Huang.

## The Model
We aim to solve the model:
$$
\left\{
\begin{aligned}
&\partial_t n_\varepsilon -D(\theta)\Delta_x n_\varepsilon -\varepsilon^2\Delta_\theta n_\varepsilon = n_\varepsilon(K(x)-\rho_\varepsilon(x)),\\
&\rho_\varepsilon = \int_0^1n_\varepsilon(x,\theta) d\theta
\end{aligned}
\right.
$$
with Neumann BC in $x$ and periodic BC in $\theta$. 

We assume that $K\in C^1(\Omega)$ and $D\in C^1([0,1))$ satisfying 
$$
K(x)\ge K_m>0, \quad D(\theta)\ge D(\theta_m):=D_m>0.
$$

### Theoretical results:
There are two limits of interest, namely $\varepsilon\to 0$ and $t\to \infty$. 

As $t\to\infty$, it will reach the steady state of  
$$
	\begin{cases}
	\begin{aligned}
&-D(\theta)\Delta_x n_\varepsilon -\varepsilon^2\Delta_\theta n_\varepsilon = n_\varepsilon(K(x)-\rho_\varepsilon(x)),\\
&\rho_\varepsilon = \int_0^1n_\varepsilon(x,\theta) d\theta
	\end{aligned}
	\end{cases}
$$

And then we let $\varepsilon\to0$, in the sense of distributions, it can be proved that 
$$
	n_\varepsilon(x,\theta) \to N_m(x) \delta(\theta-\theta_m), \quad \rho_\varepsilon \to N_m,
$$ 
where $N_m$ is the positive solution of the nonlinear Fisher-type stationary problem
$$
	\begin{cases}
	\begin{aligned}
&-D_m\Delta_x N_m  = N_m(K(x)-N_m),\\
& \frac{\partial}{\partial \nu} N_m = 0 \text{ on } \partial\Omega.
	\end{aligned}
	\end{cases}
$$
Moreover, 
$$
	\varepsilon \ln n_\varepsilon(x,\theta) \to u
$$
uniformly in $x$ and $\theta$, where $u$ is the unique 1-periodic **viscosity** solution to 
$$
\begin{cases}
\begin{aligned}
&-|\nabla_\theta u|^2 = -H(\theta, \rho), \\
& \max_\theta \, u(\theta) = 0, 
\end{aligned}
\end{cases}
$$
where $H(\theta,\rho)$ is the eigenvalue corresponding to the positive eigenfunction of 
$$
	\begin{cases}
	\begin{aligned}
&-D(\theta)\Delta_x N  = N(K(x)-\rho) + N H(\theta,\rho),\\
& \frac{\partial}{\partial \nu} N = 0 \text{ on } \partial\Omega.
	\end{aligned}
	\end{cases}
$$


## Design of an A-P scheme

### Reformulation via WKB

Set $n_\varepsilon(x,\theta) = W_\varepsilon(x,\theta) e^{u(\theta)/\varepsilon}$, then
$$
\begin{aligned}
&\partial_t n_\varepsilon = \partial_t W_\varepsilon e^{u/\varepsilon} + W_\varepsilon e^{u/\varepsilon}\dfrac{\partial_t u}{\varepsilon},\\
&\nabla_x n_\varepsilon = \nabla_x W_\varepsilon ~e^{u/\varepsilon},\\
&\Delta_x n_\varepsilon = \Delta_x W_\varepsilon ~e^{u/\varepsilon},\\
&\nabla_\theta n_\varepsilon = \nabla_\theta W_\varepsilon ~e^{u/\varepsilon} + W_\varepsilon ~e^{u/\varepsilon} \dfrac{\nabla_\theta u}{\varepsilon},\\
&\Delta_\theta n_\varepsilon = \Delta_\theta W_\varepsilon ~e^{u/\varepsilon} + 2\nabla_\theta W_\varepsilon~e^{u/\varepsilon}~\dfrac{\nabla_\theta u}{\varepsilon} + W_\varepsilon e^{u/\varepsilon}\dfrac{|\nabla_\theta u|^2}{\varepsilon^2} + W_\varepsilon e^{u/\varepsilon}\dfrac{\Delta_\theta u}{\varepsilon}.
\end{aligned}
$$
The model now reads as 
$$
\begin{aligned}
&\partial_t W_\varepsilon e^{u/\varepsilon} + \color{red}{W_\varepsilon e^{u/\varepsilon}\frac{\partial_t u}{\varepsilon}}\\
&-D(\theta)\Delta_x W_\varepsilon ~e^{u/\varepsilon} -\varepsilon^2\left(\Delta_\theta W_\varepsilon ~e^{u/\varepsilon} + 2e^{u/\varepsilon}~\dfrac{\nabla_\theta W_\varepsilon\cdot \nabla_\theta u}{\varepsilon} + \color{red}{W_\varepsilon e^{u/\varepsilon}\dfrac{|\nabla_\theta u|^2}{\varepsilon^2} + W_\varepsilon e^{u/\varepsilon}\dfrac{\Delta_\theta u}{\varepsilon}}\right)\\
& = W_\varepsilon e^{u/\varepsilon}(K(x)-\rho_\varepsilon(x)) \\
& = W_\varepsilon e^{u/\varepsilon}(K(x)-\rho_\varepsilon(x)) + W_\varepsilon e^{u/\varepsilon} H^*(\theta) \color{red}{- W_\varepsilon e^{u/\varepsilon} H^*(\theta)}.
\end{aligned}
$$
By properly separating the equation, we get
$$
	\partial_t u - \varepsilon|\nabla_\theta u|^2 - \varepsilon^2 \Delta_\theta u = -\varepsilon H^*(\theta),
$$
and 
$$
\partial_t W_\varepsilon-D(\theta)\Delta_x W_\varepsilon -\varepsilon^2\Delta_\theta W_\varepsilon - 2\varepsilon\nabla_\theta W_\varepsilon\cdot \nabla_\theta u = W_\varepsilon (K(x)-\rho_\varepsilon(x)) + H^*(\theta) W,
$$
where $H^*(\theta)=H(\theta,\rho^*)$ is the approximation of $H(\theta, \rho_\varepsilon)$ with $\rho_\varepsilon = \int_0^1n_\varepsilon(x,\theta) d\theta $.

### Idea for A-P scheme
For time evolution from $t_n$ to $t_{n+1}$:

1. To compute $\rho^n$ and $H^*$
2. To update $u^{n+1}$ and $W^{n+1}$

## Numerical treatment for conservation laws and Hamilton-Jacobi equation

Solutions of **H-J equations** 
$$
	\partial_t \phi + H(\nabla_x \phi) = 0, \quad \phi(x,0) = \phi_0(x),
$$
are **continuous** and, in the generic case, **form discontinuous derivatives in a finite time** even with smooth initial conditions.

Solutions with this kind of discontinuity are not unique.
Therefore, analogous to conservation laws, it is necessary to introduce the concept of the entropy-like condition to facilitate the selection of a unique solution, which leads to the so-called **viscosity solution**.

For **convex Hamiltonians**, the viscosity solution is characterized by a **semiconcave stability condition**. Such a viscosity solution coincides with the limit solution obtained by the **vanishing viscosity method**. For **general Hamiltonians**, the definition of the viscosity solution and the question of well-posedness (in L∞) were formulated and systematically studied by Crandall, Evans, Lions, Souganidis, and many others.

For one-dimensional H-J equations, there is a one-to-one correspondence with the conservation laws
$$
	\partial_t u + \sum_{i=1}^n \frac{\partial}{\partial x_i} f_i(u) = 0, \quad u(x,0) = u_0(x),
$$
by setting $u = \frac{\partial}{\partial x}\phi$. In the multidimensional case, however, this kind of one-to-one correspondence no longer exists. Instead,  satisfies a weakly hyperbolic systemof conservation laws [Kr, JiXi]. Instead, $\nabla_x\phi$ satisfies a weakly hyperbolic system of conservation laws. In view of these arguments, we can think of viscosity solutions of the H-J equations as primitives of entropy solutions for the conservation laws. Based on this idea, concepts used for conservation laws can be passed to H-J equations.

For the vanishing viscosity approximations with viscosity amplitude $\varepsilon$, finite difference, finite element, and finite volume solutions based on grid-cells of size $\Delta x=\varepsilon$.

For 1D conservation laws:  **The Nessyahu–Tadmor(NT) scheme**:
1. To reconstruct a **piecewise-linear (MUSCL-type) interpolant** from the known cell averages at time $t_n$: $\overline{w}_j^n \to w(x, t_n)$
2. To evolve the interpolant exactly in time, projected on the staggered cell-averages at the next time step, $t_{n+1}$, and approximated by the midpoint rule,
resulting with the two-step **predictor-corrector** form: $\overline{w}_j^n \to w_j^{n+\frac12}$ by Taylor's expansion and get $\overline{w}_{j+\frac12}^n$.

For 1D H-J equations: 


## References
1. *Rare Mutations Limit of a Steady State Dispersal Evolution Model*, B. Perthame and P. E. Souganidis, Math. Model. Nat. Phenom., 11(2016), 154-166.
2. *High-Resolution Nonoscillatory Central Schemes for Hamilton--Jacobi Equations*,  C.-T. Lin and E. Tadmor, SIAM J. Sci. Comput., 21(2000), 2163-2186.

