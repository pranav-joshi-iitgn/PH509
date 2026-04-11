# Problem 1: Time-dependent evolution of a Gaussian wavepacket in 1D

## Part (a): Numerical Method

### Split operator

#### In physical units

For a lot of systems, $V(x,t) = V(x)$. Using that, the hamiltonian operator operates entirely in the realm of position. 
For these systems, we have $\psi_t = e^{-it\hat H/\hbar}\psi_0$ . This is again similar to what we got for classical Hamiltonian systems. The full derivation for this is attached in the `Notes.md` file.  

Take a look at this equation :

$$
e^{i\frac{t}{2}\hat H/\hbar}\psi_t = \psi_{t/2}= e^{-i\frac{t}{2}\hat H/\hbar}\psi_0
$$

We can approximate this for small $t$ as :

$$
(1+i\frac{t}{2}\hat H/\hbar)\psi_t \approx (1 - i\frac{t}{2}\hat H/\hbar)\psi_0 
$$

Then, for small $\Delta x$, we can discretise the $\frac{\partial^2}{\partial x^2}$ operation and write :

$$
\psi_t(x+\Delta x) + (\frac{4m\Delta x^2}{t\hbar}i - \frac{2m\Delta x^2}{\hbar^2}V(x) -2)\psi_t(x) + \psi_t(x-\Delta x) \\ \approx \\
- \psi_0(x+\Delta x) + (\frac{4m\Delta x^2}{t\hbar}i + \frac{2m\Delta x^2}{\hbar^2}V(x) + 2)\psi_0(x) - \psi_0(x-\Delta x)
$$

Define 

$$
a(x) = \frac{4m\Delta x^2}{t\hbar}i - \frac{2m\Delta x^2}{\hbar^2}V(x) -2 \\
b(x) = \frac{4m\Delta x^2}{t\hbar}i + \frac{2m\Delta x^2}{\hbar^2}V(x) + 2
$$

Then, we can write :

$$
\begin{bmatrix}1 & a(x) & 1\end{bmatrix} \begin{bmatrix}\psi_t(x+\Delta x) \\ \psi_t(x) \\ \psi_t(x-\Delta x)\end{bmatrix} =
\begin{bmatrix}-1 & b(x) & -1\end{bmatrix} \begin{bmatrix}\psi_0(x+\Delta x) \\ \psi_0(x) \\ \psi_0(x-\Delta x)\end{bmatrix}
$$

From here, you can write the equation for the whole grid/segment using tridiagonal matrices $A,B$ as :

$$
A \vec \psi_t = B \vec \psi_0
$$

And then solve this using tridiagonal matrix solvers from `scipy` that can solve this in $O(N)$.

#### Atomic Units

For a single electron in first shell of hydrogen atom, we define:

$$
\alpha = \frac{1}{4\pi \epsilon_0}\frac{q_e^2}{\hbar c} \approx \frac{1}{137}\quad \text{ and} \quad V_1 = \frac{q_e^2}{4\pi \epsilon_0}\\
a_0 = \frac{\hbar^2}{m_eV_1}= \frac{4\pi \epsilon_0 \hbar^2}{m_e q_e^2} = \frac{\hbar}{\alpha m_e c} = 0.529\times 10^{-10} \text{ m} \\
E_0 = -V_0 = \frac{V_1}{a_0}= \frac{q_e^4 m_e}{(4\pi\epsilon_0)^2 \hbar^2} = \alpha^2m_ec^2 = 27.2 \text{eV}
$$

With these constants, we make these substitutions in the TDSE :

$$
x = \tilde x a_0 \\
V = \tilde V E_0 \\
t = \tilde t \hbar/E_0
$$

This gives us the A-dimensional Time Dependent Schrodinger Equation (ATDSE) :

$$
\boxed{
-\frac{1}{2}\frac{\partial^2 \psi}{\partial \tilde x^2} + \tilde V \psi = \tilde H \psi = i \frac{\partial \psi}{\partial \tilde t}
}
$$

Here $\tilde H = \frac{1}{E_0} \hat H = -\frac{1}{2}\frac{\partial^2}{\partial \tilde x^2} + \tilde V$

Once again, assuming $\tilde V(\tilde x,\tilde t) = \tilde V(\tilde x)$, we can write :

$$
\psi_t = e^{i\hat H t/\hbar}\psi_0 = e^{i \tilde H \tilde t }\psi_0
$$

The split operator approximation was:

$$
\psi_t(x+\Delta x) + (\frac{4m\Delta x^2}{t\hbar}i - \frac{2m\Delta x^2}{\hbar^2}V(x) -2)\psi_t(x) + \psi_t(x-\Delta x) \\ \approx \\
- \psi_0(x+\Delta x) + (\frac{4m\Delta x^2}{t\hbar}i + \frac{2m\Delta x^2}{\hbar^2}V(x) + 2)\psi_0(x) - \psi_0(x-\Delta x)
$$

This now becomes:

$$
\psi_{\tilde t}(\tilde x+\Delta \tilde x) + (\frac{4\Delta \tilde x^2}{\tilde t}i - 2\Delta \tilde x^2 \tilde V(\tilde x) -2)\psi_{\tilde t}(\tilde x) + \psi_{\tilde t}(\tilde x-\Delta \tilde x) \\ \approx \\
- \psi_0(\tilde x+\Delta \tilde x) + (\frac{4\Delta \tilde x^2}{\tilde t}i + 2\Delta \tilde x^2 \tilde V(\tilde x) + 2)\psi_0(\tilde x) - \psi_0(\tilde x-\Delta \tilde x)
$$

Similarly for the actual numerical method, now we define

$$
a(x) = \frac{4\Delta \tilde x^2}{\tilde  t}i - 2\Delta \tilde  x^2\tilde V(\tilde  x) -2 \\
b(x) = \frac{4\Delta \tilde x^2}{\tilde  t}i + 2\Delta \tilde  x^2\tilde V(\tilde  x) + 2
$$

Using that, we write :

$$
\begin{bmatrix}1 & a(\tilde x) & 1\end{bmatrix} \begin{bmatrix}\psi_t(\tilde  x+\Delta \tilde  x) \\ \psi_{\tilde  t}(\tilde  x) \\ \psi_{\tilde  t}(\tilde  x-\Delta \tilde  x)\end{bmatrix} \\ = \\
\begin{bmatrix}-1 & b(\tilde  x) & -1\end{bmatrix} \begin{bmatrix}\psi_0(\tilde  x+\Delta \tilde  x) \\ \psi_0(\tilde  x) \\ \psi_0(\tilde  x-\Delta \tilde x)\end{bmatrix}
$$

Using this, tri-diagonal matrices $A,B$ will be created, allowing us to do the implicit split-operator step by solving $A\psi_{\tilde t} = B \psi_0$ in $O(N)$ operations using `scipy`.

For further accuracy, we can even go to a penta-diagonal matrix, using 5-point evaluation of the $\frac{\partial^2}{\partial x^2}$ operator. 

### Space-Discretised Leap-Frog (SDLF)

The ATDSE is :

$$
- \frac{1}{2}\frac{\partial^2}{\partial \tilde  x^2} \psi + \tilde  V \psi= \tilde H \psi = i\frac{\partial }{\partial \tilde t}\psi
$$

Since the $\tilde H$ operator can be decomposed lineraly and maps real numbers to real numbers, we can write (after dropping the tilde) :

$$
\psi = R + iI \implies \\
H R + i H I= i \frac{\partial }{\partial t} R - \frac{\partial }{\partial t} I \implies \\
\frac{\partial}{\partial t} I  = H R \quad \text{and} \quad \frac{\partial}{\partial t} R = - H I
$$

On top of this, we have

$$
\int (I^2 + R^2) dx = 1
$$

When descritised, this equation becomes

$$
\sum_x (I_x^2 + R_x^2) = 1
$$

Similarly, the evolution equations are discretised in space as :

$$
\frac{\partial}{\partial t} I[x] = V[x]R[x] - \frac{1}{2}\frac{1}{\Delta x^2} (R[x+\Delta x]-2R[x]+R[x-\Delta x]) \\
\frac{\partial}{\partial t} R[x] = -V[x]I[x] + \frac{1}{2}\frac{1}{\Delta x^2} (I[x+\Delta x]-2I[x]+I[x-\Delta x])
$$

We can rewrite the first equation as :

$$
\frac{\partial}{\partial t} I[x] = - \frac{1}{2\Delta x^2} 
\begin{bmatrix}
1 & -2(1 + V[x]\Delta x^2) & 1
\end{bmatrix}
\begin{bmatrix}
R[x+\Delta x] \\
R[x] \\
R[x-\Delta x]
\end{bmatrix}
$$

Define $f[x] = -2(1+V[x]\Delta x^2)$ and then the tri-diagonal matrix

$$
F = \frac{1}{2\Delta x^2}\begin{bmatrix}
f[x_l] & 1 & 0 & 0 & \dots & 0 & 0  \\
1 & f[x_l + \Delta x] & 1 & 0 & \dots & 0 & 0 \\
0 & 1 & f[x_l + 2\Delta x] & 1 & \dots & 0 & 0 \\
\vdots & \vdots & \vdots & \vdots & & \vdots & \vdots \\
0 & 0 & 0  & 0 & \dots & 1 & f[x_r]
\end{bmatrix}
$$

Now, vectorising this, with $\vec R = [R_x]_x$ and $\vec I = [I_x]_x$, we can write :

$$
\frac{\partial}{\partial t}\vec R = F \vec I\\
\frac{\partial}{\partial t}\vec I = - F \vec R \\
||\vec R||^2 + ||\vec I||^2 = 1
$$

This is exactly how the hamiltonian systems we have seen evolve. The third equation is similar to the net energy being constant in a classical mechanics system.

So, just like before, we can use the position verlet LF integrator like so (for small step of size $t$):

$$
\vec R_{t/2} := \vec R_0 + \frac{t}{2} F \vec I_0 \\
\vec I_{t} := \vec I_0 - t F \vec R_{t/2} \\
\vec R_t := \vec R_{t/2} + \frac{t}{2} F\vec R_{t/2} 
$$

## Part (b): Normalization and Free Evolution

A gaussion packet with $\sigma=0.7$ centered at $x_0=-3$ was initialised and given the average initial momentum/wave-number of $k_0=2$. This was then evolved with a free potential. The evolution is visualised in the GIF animation shown. As you can see, the packet spreads in position-space as time passes, but doesn't spread in the wave-number/momentum domain. Here $\phi(k)$ is computed using FFT. The derivations and exact method can be found in `Notes.md`.  

![Free Potential Evolution of moving gaussian packet](media/P1/Part_B_Free_k0_2.gif)

To verify that the normalisation actually happens, a plot of  $N_x(t) = \braket{\psi|\psi}$ is shown in the (static) figure, along with the normalisation $N_k$ in the momentum domain.
You can see that the variation for both $N_x$ and $N_k$ is on the order of $10^{-7}$. 

![Free Potential Variation in normalisation](media/P1/Part_B_Free_k0_2_summary.png)

It can in-fact be proven analytically that the SDLF method preserves the discreate analogues of the $N_x$ and $N_k$ values. For the derivation, see the section "Leap-Frog like finite difference method" in `Notes.md`.

## Part (c): Harmonic Oscillator ($k_0 = 0$)

A gaussian packet with $\sigma=0.7$ was initialised at $x_0=-3$ and was given no average momentum, i.e. $k_0=0$. This was then evolved with the harmonic potential $V_0=\frac{1}{2}x^2$. You can see the periodic nature of the solution from the animation. Notice that here, the waveform in momentum space is also changing. This is because, unlike the free potential, with 0 classical force on the particle, here the classical force is $F(x) = -\frac{d}{dx}V = -2x$. You can also see vertical lines for 

![Animated oscilations of wave-packet in harmonic potential](media/P1/Part_C_Harmonic_k0_0.gif)

1. The most probable position $x_p(t)$ where $|\psi(x,t)|^2$ is highest
2. The mean position $\braket{x}(t)$, interpreting $|\psi|^2(x,t)$ as the probability density function. 
3. The classical solution $x_c(t)$ computed using (normal) Position-Verlet Leap-Frong integrator.

The same 3 quantities are also plotted as a function of time in a subplot of the static figure.

![Comparison of $\braket{x}$ with classical $x_c$](media/P1/Part_C_Harmonic_k0_0_summary.png)

## Part (d): Harmonic Oscillator ($k_0 = 2$)

Now, the wave-packet is given an initial momentum of $k_0=2$ rather than 0. This allows it to reach further to the right, as evident from the static plot for $\braket{x}$; similar to how $x_c$ goes further to the right than before. 
You will also notice that the average energy, namely $\braket{H}(t)$ is higher (see static figure) than what it was in the case when initial momentum was 0

![Animated oscilations of wave-packet in harmonic potential](media/P1/Part_D_Harmonic_k0_2.gif)
![Comparison of $\braket{x}$ with classical $x_c$](media/P1/Part_D_Harmonic_k0_2_summary.png)

## Part (e): Numerical Stability and Accuracy

Chapter 8 for Wang's textbook says that the SDLF integrator will give stable solutions when $\Delta t < \lambda \Delta x^2$ where $\lambda$ is a characteristic of the potential $V(x)$. It states that for a free potential, we have $\lambda=1$ and that for most potentials, $\lambda=0.5$ suffices. 
In the figure, I have plotted the trajectories for position, momentum, and the variation of $N_x,N_k$ and $\braket{H}$ with time for various $\Delta x, \Delta t$ as well as for the two methods discusses.  
I initialised the wave-function as a gaussion with $\sigma=0.7$, centered at $x_0 = -3$ and with initial average momentum of $k_0 =0$. Then I evolved the packet using the harmonic potential $V(x)=x^2/2$. 

![Effect of varying $\Delta x$ and $\Delta t$](media/P1/Part_E_Convergence.png)

You can see that for the cases when $\Delta x=0.5$, we get slight deviations of $\braket{x}$ from the classical solution $x_c$. This is expected, since the accuracy decreases on increasing $\Delta x$. 
You can also notice that the normalisation in position space, namely $N_x(t)=\braket{\psi|\psi}$ suddenly increases by a significant amount for the case with $\Delta x=0.1$ and $\Delta t=0.00825$ when using the SDLF method. If I ran the simulation for a little longer, the numerical solution with these settings would become unstable and deviate by large amounts away from the classical solution. If you look closely, you can already see the start of this instability in the curve for $\braket{x}$.
This is because we don't have $\Delta t < \lambda \Delta x^2$ in this case. The split-operator method on the other hand, doesn't face this issue and performs well even when $\Delta t$ is 0.05 (and $\Delta x =0.1$). 
Other than that, you will notice that the SDLF solution with $\Delta x=0.1$ is the closer to $x_c$ and $p_c$ curves than the split-operator solution, simply because it has a smaller $\Delta t$ (though not small enough for numerical stability). Thus, just like always, a smaller $\Delta t$ increases accuracy. Since these methods are 2nd order, the error scales as $O(\Delta t^2)$.

## Part (f): Linear Potential

I initialised the wave-packet with the same initial position and spread as all other parts, and 0 initial momentum. Then, I let it evolve under $V(x)=Fx =x$ for a longer duration so that it hits the boundary. This was done intentionally to observe the effects of the boundary. 

Since we assumed $\psi =0$ at the boundary at all times, we have effectively assumed that $V=\infty$ outside the grid/segment we are solving for. This makes the boundaries perfectly reflective, with no part of the wavefunction being transmitted forward, and all probability mass reflected backwards.
So, a wave moving towards left is reflected back. This reflected wave then interferes with the part of the original wave that is not yet reflected, giving us a wavy interference pattern. This can be seen towards the end of the animation. 

Before the particle gets close to the boundary, it has constant acceleration and thus follows the equation $\braket{x} \approx x_c = -3 -\frac{1}{2}Ft^2 = -3 -t^2/2$. This is a downward facing parabola, as evident from the subplot for position in the static figure. Note that even now, there is some effect from the boundary, but it isn't large enough for it to be visible.
Another thing to note is that the average momentum also changes linearly, similar to the the classical solution.
The wave-packet also spreads in position with time (see GIF animation), similar to what we saw in the case for the free particle. 

![Animation showing boundary effects and spreading for linear potential](media/P1/Part_F_Linear_k0_0.gif)

![Curves for $x_c(t)$ and $\braket{x}(t)$](media/P1/Part_F_Linear_k0_0_summary.png)

# Problem 2: Radial Schrodinger equation for the hydrogen atom

## Part (a): Effective Potential



## Part (b): Bound-State Energies
## Part (c): Radial Wavefunctions
## Part (d): Energy Degeneracy
## Part (f): Degeneracy Breaking

# Problem 3: 2D Quantum Scattering

## Part (a): Boundary Conditions and Normalization
## Part (b): Free Propagation
## Part (c): Repulsive Barrier
## Part (d): Parameter Variations
## Part (e): Angular Distribution
## Part (f): Attractive Potential
## Part (g): Convergence Tests
## Part (h): Comparison with Classical Expectation