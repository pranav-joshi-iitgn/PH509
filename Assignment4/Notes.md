# Chapter 8 : Time-Dependent Quantum mechanics

## TDSE

The one-dimensional time-dependent Schrodinger equation is given by :

$$
\frac{-\hbar^2}{2m} \frac{\partial^2}{\partial x^2} \psi(x,t) + V(x,t)\psi(x,t) = i\hbar \frac{\partial}{\partial t}\psi(x,t)
$$

## $V(x,t)=0$

In this case, we have :

$$
\begin{aligned}
\frac{-\hbar^2}{2m} \frac{\partial^2}{\partial x^2} \psi &= i\hbar \frac{\partial}{\partial t} \psi \\
\implies \frac{\partial}{\partial t}\psi= (\frac{i\hbar}{2m})\frac{\partial^2}{\partial x^2}\psi
\end{aligned}
$$

By separation of variables, i.e. $\psi(x,t) = a(x)b(t)$, we get :

$$
a(x)b'(t) = (\frac{i\hbar}{2m}) a''(x)b(t)
$$

Then, for all $(x,t)$ where $a\ne 0$ and $b \ne 0$, we have:

$$
\frac{b'(t)}{b(t)} = (\frac{i\hbar}{2m}) \frac{a''(x)}{a(x)} \quad [= f (x,t)]
$$

Since LHS is a constant w.r.t $x$ and RHS is constant w.r.t $t$, we can say set both equal to some _constant_ $-i \omega$. So, we get :

$$
b'(t) = -i\omega b(t) \\
\frac{\hbar}{2m} a''(x) = -\omega a(x)
$$

These equations are easy to solve. We get the solution

$$
b(t) = \beta e^{-i\omega t} \\
a(x) = \alpha_1 e^{i\sqrt{\frac{2m}{\hbar}\omega} x} + \alpha_2 e^{-i\sqrt{\frac{2m}{\hbar} \omega} x}
$$

Now, for stability and boundedness of these solutions, it is necessary that $\omega\in \mathbb{R}^{+}$. 

Of-course, this is only _a_ solution, and we can have as many as we want for all the different values of $\omega$. 

Therefore, we have a _basis_ of solutions, namely:

$$
\psi_k(x,t) = A_k e^{i(kx -\omega t)} \quad \text{where}\quad \omega=\frac{\hbar}{2m}k^2 \quad \text{and} \quad k\in \mathbb{R}
$$

This gives the general solution:

$$
\psi(x,t) = \int_{-\infty}^{\infty} [\frac{1}{\sqrt{2\pi}} \phi(k)] e^{i(kx-\omega t)} dk
$$

Set $t=0$ (which can be _any_ point you want really) and this looks like an inverse fourier transform.

From this, it's easy to see that:

$$
\phi(k) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} \psi(x,0) e^{-ikx}dx
$$

Now, since we did the whole setting $t=0$ think without realling caring about what point in time we are choosing as $t=0$, we effectively dissolved the $e^{-i\omega t}$ into $\phi(k)$. So, it's better to actually think of it that and not do the $t=0$ trick. This gives us :

$$
\psi(x,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} \phi(k,t) e^{ikx}dk = {\cal F}^{-1}[\phi] \\
\phi(k,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} \psi(x,t) e^{-ikx}dx = {\cal F}[\psi] \\
$$

If we use the 2nd equation to _define_ $\phi$, then these two equations are mathematically true in general, not just for $V(x,t)=0$

The only think special about the $V(x,t)=0$ case is that we can write 

$$
\phi(k,t) = \phi(k,0)e^{-i \frac{\hbar}{2m}k^2 t}
$$

That is to say, once the initial conditions are given, we have an _analytical_ solution for the evolution. 

## Momentum and Position

In quantum mechanics, any observable quantity $O(t)$ is represented by an Hermitian operator $\hat O$, i.e. one that has the property $\braket{ \hat O \psi_1 | \psi_2} = \braket{\psi_1 |\hat O \psi_2}$ where the "bra-ket" notation is defined as :

$$
\braket{\psi_1|\psi_2} = \int_{-\infty}^\infty \psi_1^*(x,t)\psi_2^*(x,t) dx
$$

In fact, we in particular have the expectation (mean) value:

$$
\braket{O}(t) = \bra {\psi(t)} \hat O \ket{\psi(t)} = \braket{\psi| \hat O \psi}
$$

For position, the operator is $\hat x = x$. Thus,

$$
\braket{x}(t) = \int_{-\infty}^\infty \psi^*(x,t) \,x\,\psi(x,t) \,dx = \int_{-\infty}^\infty x |\psi(x,t)|^2 dx = E[x](t)
$$

For momentum $p$, we have $\hat p =\frac{\hbar}{i}\frac{\partial}{\partial x}$. Thus, we get :

$$
\braket{p} = -i\hbar\int_{-\infty}^\infty \psi^*(x,t)\frac{\partial}{\partial x}\psi(x,t) dx = -i\hbar \bra{\psi}\ket{\frac{\partial}{\partial x}\psi}
$$

But remember that $\psi(x,t)$ can be written in terms of $\phi(k,t)$. This gives :

$$
\frac{\partial}{\partial x} \psi = \frac{1}{\sqrt{2\pi}}\int_{-\infty}^{\infty}\phi(k,t)\frac{\partial}{\partial x} e^{ikx}dk \\ = \frac{i}{\sqrt{2\pi}}\int_{-\infty}^{\infty} k \phi(k,t) e^{ikx}dk = i{\cal F}^{-1}[k\phi(k,t)]
$$

Substituting this back, we get :

$$
\braket{p}  = \hbar \braket{\psi |{\cal F}^{-1} [k \phi(k,t)]}
$$

Since the functions we are dealing with are all bounded are analytic, we can switch the order of integration like this :

$$
\braket{p}  = \hbar \frac{1}{\sqrt{2\pi}}\int_{-\infty}^\infty [\int_{-\infty}^{\infty}\psi^*(x,t) k \phi(k,t)e^{ikx}dx]dk \\
= \hbar \int_{-\infty}^\infty [\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{\infty}\psi^*(x,t) e^{ikx}dx\,]\; k \phi(k,t) dk \\ 
= \hbar \int_{-\infty}^\infty [\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{\infty}\psi^*(x,t) e^{ikx}dx\,]\; k \phi(k,t) dk \\
= \hbar \int_{-\infty}^\infty {\cal F}^{-1}[\psi^*](k,t)\; k\; \phi(k,t) \,dk
$$

Now, using the conjugation property of fourier transform, i.e. "fourier transform of complex conjugate is complex conjugate of fourier transform", we can write :

$$
\braket{p} = \hbar \int_{-\infty}^{\infty} ({\cal F}^{-1} [\psi])^* k \phi(k,t) dk \\ = \hbar \int_{-\infty}^{\infty} \phi^*(k,t) k \phi(k,t) dk \implies \\
\braket{p} = \hbar \int_{-\infty}^{\infty} k |\phi(k,t)|^2 dk
$$

Now, from modern physics, we know that for a particle moving with pomentum $p$, its wavelength is given by $\lambda = h/p$ . This gives the wave number $k = 2\pi/\lambda = 2\pi p/h = p/\hbar$. 
The solution for $V(x,t)=0$ is basically a superposition of many plane waves of the form $A_ke^{i(kx-\omega t)}$ . So of course, we can interpret the $``k"$ in our equations as the wave-number associated with the De-broglie wavelength. 

This leads us to saying that the momentum of each plane wave is given by $p = \hbar k$. Thus, we should rewrite the equation in a more suggestive form as :

$$
\braket{p} = \int_{-\infty}^\infty p |\tilde \psi(p,t)|^2 dp \quad \text{where} \quad p = \hbar k \quad \text{and}\quad \tilde \psi(p,t) = \frac{1}{\hbar} \phi(p/\hbar,t)
$$

Here, $\tilde \psi(p,t)$ is what is called the momentum wave function. 

Of course, this is more general and not just for a free particle, although the interpretation is most accurate for that case.

## Hamiltonian operator

Define the operator

$$
\hat H = \frac{\hat p^2}{2m} + V = -\frac{\hbar ^2}{2m} \frac{\partial^2}{\partial ^2 x} + V
$$

With this, we can write the TDSE as :

$$
\hat H \psi = i\hbar \frac{\partial}{\partial t}\psi \iff 

\frac{\partial}{\partial t}\psi = \frac{-i}{\hbar} \hat H \psi
$$

This is similar to the Liouville operator $D_H$ in classical mechanics.

## Evolution with time-independent potential

For a lot of systems, $V(x,t) = V(x)$. Using that, the hamiltonian operator operates entirely in the realm of position. Thus, we can write this :

$$
(\frac{\partial}{\partial t})(\hat H) = (\hat H)(\frac{\partial}{\partial t})
$$

Thus, we can write:

$$
\frac{\partial^n}{\partial t^n}\hat H = \hat H \frac{\partial^n}{\partial t^n} \implies \\

\frac{\partial^n}{\partial t^n} \psi = \frac{\partial^{n-1}}{\partial t^{n-1}} (\frac{\hat H}{i\hbar})\psi = (\frac{\hat H}{i\hbar}) \frac{\partial^{n-1}}{\partial t^{n-1}} \psi \implies 
\frac{\partial^n}{\partial t^n} \psi = (\frac{\hat H}{i\hbar})^n\psi
$$

Then, taylor expanding $\psi(x,t)$ in time, we have :

$$
\psi_t(x) = \psi(x,t) = \frac{1}{0!}\psi(x,0) + \frac{1}{1!}t[\frac{\partial }{\partial t}\psi](x,0) + \frac{1}{2!}t^2[\frac{\partial^2 }{\partial t^2}\psi](x,0) + \dots \\ = 
[\frac{1}{0!}\psi](x,0) + [\frac{1}{1!}(\frac{-i}{\hbar}{t\hat H})\psi](x,0) + [\frac{1}{2!}(\frac{-i}{\hbar}{t\hat H})^2\psi](x,0) + \dots \\ = 
[[\frac{1}{0!} + \frac{1}{1!}(\frac{-i}{\hbar}{t\hat H}) + \frac{1}{2!}(\frac{-i}{\hbar}{t\hat H})^2 + \dots] \psi](x,0) \\ 
= [e^{-i t\hat H/\hbar}\psi](x,0) = e^{-i t\hat H/\hbar}\psi_0
$$

Thus, we have $\psi_t = e^{-it\hat H/\hbar}\psi_0$ . This is again similar to what we got for classical Hamiltonian systems. 

## Conserved quantities

For any observable with operator $\hat O$, it is conserved if it commutes with $\hat H$, i.e. $\hat H\hat O = \hat O \hat H$, so that we may write $e^{-i\frac{\hat H}{\hbar} t}\hat O \psi_0 = \hat O e^{-i\frac{\hat H}{\hbar} t} \psi_0 = \hat O \psi_t$ and then $\bra{\psi_0}\hat O\ket{\psi_0} = \bra{\psi_t}\hat O\ket{\psi_t}$ . 

## Solution of TDSE for time-independent potentials

Now, since we have the case that $V(x,t) = V(x)$ and thus $\hat H$ is an operator in spatial domain, we can use separation of variables to do some nice stuff. First, consider _a_ solution to the TDSE of the form $\psi(x,t) = a(x)b(t)$ . We can write :

$$
\hat H \psi = b(t)(\hat H a)(x) = i\hbar \frac{\partial}{\partial}\psi = i\hbar a(x) b'(t)
$$

Then again, when $\psi(x,t)\ne 0$ we have :

$$
\frac{1}{a}\hat H a =i\hbar b'/b = E
$$

where $E$ is a constant. 
This gives :

$$
\hat H a = E a \\
b' = (-i\frac{E}{\hbar})b
$$

While the first equation doesn't have a general solution, the solutions of the second one is 

$$
b(t) = \beta e^{-it\frac{E}{\hbar}}
$$

Define $\psi(x) = \beta a(x)$, then we have :

$$ 
\psi(x,t) = \psi(x) e^{-itE/\hbar} \\
\hat H \psi(x) = E \psi(x)
$$

The second equation is called the time-independent Schrodinger Equation (TISE). Looking closely, the solution $\psi(x)$ for the TISE will be an _eigenfunction_ for the operator $\hat H$ with eigenvalue of $E$, which obviously is the Energy of that particular solution. Also notice that for the solution to have the property $\braket{\psi_t|\psi_t} =1 \;\forall t$, we need $E \in \mathbb{R}$. A complex value energy also doesn't make any physical sense. Fortunately, this is guaranteed since $\hat H$ is a hermitian operator and can only have real valued eigenvalues. Moreover since $\hat H$ is a hermitian operator, any two eigenfunctions must either have the same eigenvalue and inner product of 1, or must be orthogonal, i.e. $E_1 \ne E_2 \implies \braket{\psi_{E_1}|\psi_{E_2}}= 0 $ and $E_1=E_2 \implies \braket{\psi_{E_1}|\psi_{E_2}} = 1$.

Again, this will give us _a_ solution. The good thing is that most of th times, any initial condition, i.e. $\psi(x,0)$ can be expressed as a linear combination of the (orthonormal) eigenfunctions of $\hat H$. Thus, we simply have the general solution given as :

$$
\psi(x,t) = \sum_{-\infty}^{\infty} A_E \psi_E(x) e^{-itE/\hbar} \\
\text{where}\quad A_E = \braket{\psi_E(x)|\psi(x,0)} 
$$

Note that although I denote the basis functions as $\psi_E$, there can actually be mutiple functions for the same energy level. Using the properties of the eigenbasis, we can choose any function for the inner product in 2nd equation. The only thing compromised will be the uniqueness of the linear combination. 

## "Split operator" method

Take a loot at this equation :

$$
e^{i\frac{t}{2}\hat H/\hbar}\psi_t = \psi_{t/2}= e^{-i\frac{t}{2}\hat H/\hbar}\psi_0
$$

We can approximate this for small $t$ as :

$$
(1+i\frac{t}{2}\hat H/\hbar)\psi_t \approx (1 - i\frac{t}{2}\hat H/\hbar)\psi_0 
$$

Then, for small $x$, we have :

$$
[\frac{\partial ^2}{\partial x^2}\psi](0,t) \approx
\frac{1}{x^2}[\psi(x,t) - 2\psi(0,t) + \psi(-x,t)]
$$

Using this too, we get :

$$
(1 + i\frac{t}{2\hbar}V(0))\psi(0,t) - i\frac{t\hbar}{4m}\frac{1}{x^2}[\psi(x,t)-2\psi(0,t) + \psi(-x,t)] \\ = \\
(1 - i\frac{t}{2\hbar}V(0))\psi(0,0) + i\frac{t\hbar}{4m}\frac{1}{x^2}[\psi(x,0)-2\psi(0,0) + \psi(-x,0)] 
$$

Thus,

$$
(\frac{4mx^2}{t\hbar}i - \frac{2mx^2}{\hbar^2}V(0))\psi(0,t) + \psi(x,t) -2\psi(0,t) + \psi(-x,t) \\ \approx \\
(\frac{4mx^2}{t\hbar}i + \frac{2mx^2}{\hbar^2}V(0))\psi(0,0) - \psi(x,0) + 2\psi(0,0) - \psi(-x,0)
$$

Now, centering abount a non-zero $x$ value, we get :

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

We can also extend this to 2D by assuming $\Delta y = \Delta x$ and defining 

$$
a(x,y) = \frac{4m\Delta  x^2}{t\hbar}i - \frac{2m\Delta  x^2}{\hbar^2}V(x,y) -4 \\
b(x,y) = \frac{4m\Delta x^2}{t\hbar}i + \frac{2m\Delta x^2}{\hbar^2}V(x,y) + 4
$$

and writing

$$
\begin{bmatrix}
0 & 1 & 0 \\
1 & a(x) & 1 \\
0 & 1 & 0
\end{bmatrix} 
\circ
\begin{bmatrix}
\psi_t(x+\Delta x,y+\Delta y) & \psi_t(x,y+\Delta y) & \psi_t(x-\Delta x,y+\Delta y) \\
\psi_t(x+\Delta x,y) & \psi_t(x,y) & \psi_t(x-\Delta x,y) \\
\psi_t(x+\Delta x,y-\Delta y) & \psi_t(x,y-\Delta y) & \psi_t(x-\Delta x,y-\Delta y) \\
\end{bmatrix} \\
= \\
\begin{bmatrix}
0 & -1 & 0 \\
-1 & b(x) & -1 \\
0 & -1 & 0
\end{bmatrix} 
\circ
\begin{bmatrix}
\psi_0(x+\Delta x,y+\Delta y) & \psi_0(x,y+\Delta y) & \psi_0(x-\Delta x,y+\Delta y) \\
\psi_0(x+\Delta x,y) & \psi_0(x,y) & \psi_0(x-\Delta x,y) \\
\psi_0(x+\Delta x,y-\Delta y) & \psi_0(x,y-\Delta y) & \psi_0(x-\Delta x,y-\Delta y) \\
\end{bmatrix} 
$$

Solving this is going to be a nightmare.

So instead, we approximate even more, giving us what is known as the _Alternat Direction Implicit_ method. 

## ADI

Here, we split the hamiltonian operator as $\hat H = \hat H_x + \hat H_y$ where 

$$
\hat H_x = -\frac{\hbar^2 \hat p_x^2}{2m} + V_x(x,y) \\
\hat H_y = -\frac{\hbar^2 \hat p_y^2}{2m} + V_y(x,y)
$$

It it's possible to separate the potential as $V(x,y) = V_x(x) + V_y(y)$, we can simply split it. Otherwise, we can do $V_x = V_y = V/2$

Now, for the method, we discretise the wave-function over a grid, giving us the matrix $\Psi[y,x]$. Just like before, we create the tri-diagonal matrices (using $t/2$ rather than $t$), but now for each dimension separately. Then we solve :

$$
A_{x}(y)\Psi_{t/2}[y] = B_x(y)\Psi_0[y]
$$

for all $y$ and then

$$
A_y(x)\Psi_{t}^T[x] = B_y(x)\Psi_{t/2}^T[x]
$$

for all $x$.

This effectively takes $O(MN)$ operations each iteration, assuming the grid is $M\times N$ points big.

## Atomic Units

For a single electron in first shell of hydrogen atom, we define:

$$
\alpha = \frac{1}{4\pi \epsilon_0}\frac{q_e^2}{\hbar c} \approx \frac{1}{137}\quad \text{ and} \quad V_1 = \frac{q_e^2}{4\pi \epsilon_0}\\
a_0 = \frac{\hbar^2}{m_eV_1}= \frac{4\pi \epsilon_0 \hbar^2}{m_e q_e^2} = \frac{\hbar}{\alpha m_e c} = 0.529\times 10^{-10} \text{ m} \\
E_0 = -V_0 = \frac{V_1}{a_0}= \frac{q_e^4 m_e}{(4\pi\epsilon_0)^2 \hbar^2} = \alpha^2m_ec^2 = 27.2 \text{eV}
$$

Then, for an electron, we write the TDSE equation :

$$
-\frac{\hbar^2}{2m_e} \frac{\partial^2}{\partial x^2} \psi + V \psi = i\hbar \frac{\partial}{\partial t} V
$$

First we substitute $x = \tilde x a_0$ and $V = \tilde V E_0$, giving :

$$
-\frac{\hbar^2}{2m_ea_0^2} \frac{\partial^2}{\partial \tilde x^2} \psi + E_0\tilde V \psi = i\hbar \frac{\partial}{\partial t} \psi
$$

We use the fact that $E_0a_0^2 = a_0V_1 = \hbar^2/m_e$, giving us :

$$
-\frac{1}{2}\frac{\partial^2 \psi}{\partial \tilde x^2} + \tilde V \psi = i \frac{\partial \psi}{\partial t} (\frac{\hbar}{E_0})
$$

Then, we define $\tilde t = t E_0/\hbar$

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
