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

## Leap-Frog like finite difference method

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

Effectively what is happening is that

$$
\begin{bmatrix}
\vec R_{t/2} \\
\vec I_{t/2}
\end{bmatrix} = 
\begin{bmatrix}
\textbf{I}_N & \frac{t}{2} \textbf{F}\\
\textbf{0}_N & \textbf{I}_N
\end{bmatrix}
\begin{bmatrix}
\vec R_{0} \\
\vec I_{0}
\end{bmatrix} \\

\begin{bmatrix}
\vec R_{t} \\
\vec I_{t}
\end{bmatrix} = 
\begin{bmatrix}
\textbf{I}_N & \textbf{0}_N \\ 
- \frac{t}{2} \textbf{F}  & \textbf{I}_N
\end{bmatrix}
\begin{bmatrix}
\vec R_{t/2} \\
\vec I_{t/2}
\end{bmatrix}
$$

It's easy to show that both the matrices have a determinant of 1. Thus, the L2 norm is exactly preserved. 

But that isn't enough for stability. Instead, for von Neuman or BIBO stability, we require :

$$
t < \Lambda \Delta x^2
$$

Here, $\Lambda$ is a quantity dependent on $V(x)$. For free fall, i.e. $V(x)=0$, we have $\Lambda =1$ . For most other potentials we will encounter, we will have $\Lambda \ge 1/2$ . 

Thus, as long as $t \ll \frac{1}{2}\Delta x^2$, the simulation will give stable outputs. 

Obviously, we won't do the iteration for just the $[0,t]$ interval in time, but rather till $T = t\times N_f$ where $N_f$ is number of frames. 

The time complexity using either band-matrix multiplication of vector dot products after slicing and padding is $O(NN_f)$.
If $N_f \gg N$, we can optimise by doing fast matrix exponentiation of the matrix

$$
M := 
\begin{bmatrix}
\textbf{I}_N & \textbf{0}_N \\ 
- \frac{t}{2} \textbf{F}  & \textbf{I}_N
\end{bmatrix}
\begin{bmatrix}
\textbf{I}_N & \frac{t}{2} \textbf{F}\\
\textbf{0}_N & \textbf{I}_N
\end{bmatrix}
$$

to get $M^{N_f}$ in $O(N^3 \log N_f)$ time. 
For example, if we have $(x_l,x_r) = -10,10$ and $\Delta x = 0.1$, giving us $N = 201$, but have $t = 0.001 < \frac{1}{2}\Delta x^2$ and we evolve till $T_1=1$, giving us $N_f = 1000$. It may seem that $N^3 \gg N N_f$ and so this is a waste of time, but pre-computing $M^{N_f}$ and repeatedly using it allows us to speed up the simulation by computing effect of $N_f$ iterations via just $O(N^2)$ operations, rather than $O(NN_f)$. So, if the simulation has to be done for a longer time, say $T=1000$, we can just use the matrix for $T=10$ and use that to do 100 of 160000-operation matrix-vector multiplications, rather than doing $2\times 3 \times 201000000$ operations. 

This also allows us to make the $t$ value really small. 

Eventually, there will come a point when

$$
M^{T/t} = (I+ \frac{t}{T}(\frac{M-\textbf{I}_n}{t/T}))^{T/t} \approx \exp(\frac{T}{t}(M-\textbf{I}_n))
$$

This quantity can be computed using `scipy.linalg.expm` in $O(N^3)$ and is basically the space discretised analogue of the $e^{-iHt}$ operator. Notice I said "space discretised". In the time dimension, it's effectively the solution to the equation

$$
\frac{\partial}{\partial t} \begin{bmatrix}
\vec R \\
\vec I
\end{bmatrix} = \begin{bmatrix}
0 & F \\
-F & 0
\end{bmatrix} \begin{bmatrix}
\vec R \\
\vec I
\end{bmatrix}
$$

i.e.

$$
\begin{bmatrix}
\vec R \\
\vec I
\end{bmatrix}(T) = \exp(T G) (\begin{bmatrix}
\vec R \\
\vec I
\end{bmatrix}(0)) \\
\text{where}\quad  G = \begin{bmatrix}
0 & F \\
-F & 0
\end{bmatrix}
$$

It can be shown that as $t\to 0$, the $\frac{T}{t}(M-I)$ matrix become $G$.
This then raises the obvious question : why bother with leap-frog at all where we can get the time-analytical solution for the space discretised scenario directly ? 

Since this has become a very different method from space-discretised leap-frog (SDLP) after all the optimisations and analysis, I will give it a new name: the "Real-Imaginary Split Evolution" (RISE) method. 

The algorithm for RISE given $\Delta x, \Delta t$ where $\Delta t$ is large, and $\psi_0, V$ is :

* Compute $f[x]$ from $V[x],\Delta x, \Delta t$
* Compute $F$ and $G$ matrices
* Compute the transition matrix $P = \exp(G\Delta t)$
* Initialise $R = \Re (\psi_0)$ and $I = \Im(\psi_0)$
* for each frame:
    * Compute any quantities from $\psi = R + iI$ that are to be tracked (such as $\braket{x},\braket{p},\braket{E}, \phi$, etc.)
    * Update $R,I$ as $[R\;I]^T := P[R\;I]^T$

Although this method gives analytically correct time evolution, the $O(N^3)$ bottleneck means that the spatial resolution is poor for 2D and 3D cases, and thus the solution will also be in-accurate.
Thus, this method isn't used professionally. 

But, for the toy examples we have, I can certainly use it to verify correctness. 