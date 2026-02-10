# Symplecitc Integrators

## Hamiltonian systems

Setting for Hamiltonian mechanics is a $2n$ dimensional coordinate system $\vec z = (\vec q,\vec p) = (q_1,q_2\dots q_n,p_1\dots p_n)$
The Hamiltonian $H(q,p)$ of the system is specified before-hand. 
As the system evolves, $\vec p (t), \vec q(t)$ trace out a curve with the tangent in the direction of $(\nabla_{\vec p} H, -\nabla_{\vec q} H)$. Essentially, we have : 

$$
\frac{d\vec q}{dt} = \nabla_{\vec p} H \\
\frac{d\vec p}{dt} = -\nabla_{\vec q} H
$$

This means, that we have :

$$
\frac{d}{dt}H = \frac{d\vec q}{dt}\cdot \nabla_{\vec q} H + \frac{d\vec p}{dt} \cdot \nabla_{\vec p} H \\
= (\nabla_{\vec p} H)\cdot (\nabla_{\vec q} H) + (-\nabla_{\vec q} H) \cdot (\nabla_{\vec p} H) = 0
$$

Now, in terms of $\vec z$, we can write the evolution equations as :

$$
\frac{d\vec z}{dt} = J \nabla_{\vec z} H
$$

where $J$ is the matrix constructed as :

$$
J = \begin{bmatrix}
\boldsymbol{0} & I_n \\
- I_n & \boldsymbol{0}
\end{bmatrix}
$$

Notice that $J$ is skew symmetric. 

Suppose we have $\vec v = (\vec v_1,\vec v_2)$ and $\vec w = (\vec w_1,\vec w_2)$, then we define the "linear symplectic structure" as :

$$
[\vec v, \vec w] = (J\vec v)\cdot \vec w = \vec v_2 \cdot \vec w_1 - \vec v_1\cdot \vec w_2
$$

One can rewrite this as :

$$
[\vec v, \vec w] = \sum_{i=1}^n (v_{2,i}w_{1,i} - v_{1,i}w_{2,i}) = \sum_{i=1}^n \begin{vmatrix}w_{1,i} & v_{1,i} \\ w_{2,i} & v_{2,i}\end{vmatrix}
$$

Thus, the $[\cdot, \cdot]$ operation gives a kind of "area". 

## Phase flow

When one solves the Hamiltonian equations with a given initial value $z_0$, we get a "phase flow function" $\psi(t,\vec z_0)$ or $\psi_t(\vec z_0)$ that gives the value of $z$ at any given $t$. 

If we define $\Psi'_t(\vec z_0) = \frac{d \vec \psi_t(\vec z_0)}{d\vec z_0}$ as the Jacobian, then we can write : 

$$
[\boldsymbol{\Psi_t'}(\vec z_0) \vec v, \boldsymbol{\Psi_t'}(\vec z_0) \vec w] = [\vec v,\vec w] \;\forall \vec v,\vec w \in \mathbb{R}^{2n}
$$

This can be equivalently written as :

$$
[\Psi_t'(\vec z_0)]^T J [\Psi_t'(\vec z_0)] = J
$$

Any matrix/operator that satisfies the equation $M^TJM = J$ is called a symplectic operator. 

## Poisson Bracket

For scalar functions $f(\vec z), g(\vec z)$, the poisson braket is defined as :

$$
\{f,g\} = [\nabla_{\vec z} \vec g, \nabla_{\vec z} f] = \sum_{i=1}^n (\frac{f}{\vec q_i}\cdot \frac{g}{\vec p_i} - \frac{f}{\vec p_i}\cdot \frac{g}{\vec q_i})
$$

For a given scalar function $f(\vec z)$ and the hamiltonian $H(\vec z)$, we can define the Liouville operator as :

$$
D_H f = \{f H\} = [\nabla H, \nabla f] = (J \nabla H)\cdot \nabla f 
$$

Now, let's evaluate $\frac{d}{dt} f(\vec z(t))$ using the chain rule. We have :

$$
\frac{d}{dt}f = \nabla_{\vec z} f \cdot \frac{d}{dt}\vec z = \nabla_{\vec z}f \cdot (J\nabla_{\vec z}H) = D_H f
$$

Now, since $D_H f$ is also just another function of $\vec z$, one can say 

$$
\frac{d^2}{dt^2}f = \frac{d}{dt}(\frac{df}{dt}) = \frac{d}{dt}(D_H f) = D_H (D_H f) = D_H^2 f
$$

and then by induction, 

$$
\frac{d^k}{dt^k}f(\vec z) = D_H^k f
$$

For a vector, we just apply this operator to each component point-wise. That returs a vector of the same size as the input. 

Thus, we have the identity

$$
\boxed{
\frac{d^k}{dt^k}f(\vec z(t)) = D_H^k \vec f
}
$$

Of course, this also applies to the function $\vec f(\vec z) = \vec z$ which gives us $\frac{d^k}{dt^k}\vec z(t) = D_H^k \vec z$

As a reminder, the operator $D_H$ is dependent on the current value of $\vec z$ . 

Now, one can write the taylor series solution

$$
\vec z(t) = \vec z (0) + t(\frac{d}{dt}\vec z)(0) +  t^2(\frac{1}{2!}\frac{d^2}{dt^2}\vec z)(0) + t^3(\frac{1}{3!}\frac{d^3}{dt^3}\vec z)(0) + \dots \\ = 
\vec z(0) + t(D_H\vec z)(0) + \frac{1}{2!}t^2(D_H^2 \vec z)(0)  + \dots \\
= ((1 + (tD_H) + \frac{1}{2}(tD_H)^2 + \dots)\vec z)_{\vec z = \vec z_0} \\= \exp(tD_H) z_0 
$$

Trying to wrap your head around the $\exp(tD_H)$ notation isn't worth it. Let's just accept this monstrousity and move on. 
What we have effectively done is "solved" the hamiltonian differential equation with the solution

$$
\vec \psi(t,\vec z_0) = [\exp(tD_H) \vec z]_{\vec z=\vec z_0}
$$

## Assumption on form of $H$

Consider hamiltonians of the form

$$
H(\vec z) = V(\vec q) + T(\vec p)
$$

This allows us to write :

$$
D_H f = (\nabla_{\vec p}T(\vec p))\cdot  \nabla_{\vec q} f - (\nabla_{\vec q}V(\vec q)) \cdot \nabla_{\vec p} f
$$

Or with better notation, 

$$
D_H = (\nabla_{\vec p}T(\vec p))^T \nabla_{\vec q} - (\nabla_{\vec q}V(\vec q))^T \nabla_{\vec p}
$$

Now, let's define

$$
D_T = (\nabla_{\vec p}T(\vec p))^T \nabla_{\vec q} \\
D_V = (\nabla_{\vec q}V(\vec q))^T \nabla_{\vec p}
$$

which allows us to write $D_H = D_T - D_V$

Now, unfortunatly $D_T D_V \ne D_V D_T$ and so we cannot use $\exp(tD_H) = \exp(tD_T)\exp(-tD_V)$.

## Small step approximation

But the good thing is that 

$$
\exp(t(D_T - D_V)) = \exp(tD_T)\exp(-tD_V) + O(t^2)
$$

So, up-to a first order, we can basically apply the solution ($\exp(-tD_V)$) of $\frac{d}{dt}\vec z = -D_V \vec z$ first and then the solution $\exp(tD_T)$ of $\frac{d}{dt}\vec z = D_T \vec z$ . Now since both $\exp(tD_T)$ and $\exp(-tD_V)$ are symplectic operators, applying one after another produces another symplectic operator. 

This is enough to get a first order method. Since you'll be having a $O(t^2)$ of error anyway, it makes sense to evaluate the taylor series for $\exp(tD_V)$ and $\exp(tD_T)$ only to the first order terms. 


This leads to the Euler Cromer method. 

Then, if you want an $O(t^2)$ method, we can use 

$$
\exp(tD_H) = 1 + tD_H + t^2 D_H^2/2 + O(t^3) \\ = 1 + tD_T - tD_V + \frac{1}{2}t^2(D_T^2 - D_TD_V - D_VD_T + D_V^2) + O(t^3) \\
\exp(tD_T/2)\exp(-tD_V)\exp(tD_T/2) = \\ 1 + (tD_T - tD_V) + (\frac{t^2}{4}D_T^2 + \frac{t^2}{2}D_V^2) + (- t^2D_TD_V/2 - t^2D_VD_T/2 + t^2D_T^2/4) + O(t^3) \\
\implies \exp(tD_H) = \exp(tD_T/2)\exp(-tD_V)\exp(tD_T/2) + O(t^3)
$$

Now, simulating this evolution using a second order method for each individual evolutions gives us a symplecitic second order method, namely the Leap Frog method. 

Now, all this is good, but it's a bit of overkill.

The equations we are interested in are :

$$
D_T\vec z = (\nabla_{\vec p}T(\vec p))^T \nabla_{\vec q}\vec z = \begin{bmatrix}\nabla_{\vec p} T \\ \vec 0 \end{bmatrix}\\
D_V\vec z = (\nabla_{\vec q}V(\vec q))^T \nabla_{\vec p}\vec z =  \begin{bmatrix} \vec 0 \\ \nabla_{\vec q} V\end{bmatrix}\\
$$

Going a bit further:

$$
D_T^2 \vec z = \vec 0_{2n} \\
D_V^2\vec z = \vec 0_{2n}
$$

What this means is that the separated equations actually have very simple solutions, namely:

$$
\vec z_T(t) = \vec z(0) + t\begin{bmatrix}(\nabla_{\vec p} T)(\vec p_{0}) \\ \vec 0 \end{bmatrix} \\
\vec z_V(t) = \vec z(0) + t\begin{bmatrix}\vec 0 \\(-\nabla_{\vec q} V)(\vec q_{0}) \end{bmatrix}
$$

So, we can very well just use Euler's method for simulating the evolution $\exp(-D_Vt)$ and $\exp(D_T t/2)$

Note the the final goal is still to simulate $\exp(D_H t)$ which, doing it this way will still give us an $O(t^3)$ error. 

Note : When I say $O(t^k)$, it should actually be interpreted as $O(\Delta t^k)$ . And when I say $\vec z_0$ and $\vec z(t)$, it should be interpreted as $\vec z_n$ and $\vec z_{n+1}$ in the context of the algorithm. 

## Energy conservation

Finally, we wanted methods that conserve energy on the system.
While symplectic solvers try to preserve the mathematical properties of the system, they do not necessarily preserve the enrgy at all times. 

This is because every symplectic solver is a combination of exact solutions applied for small steps for systems with a _different_ hamiltonian ($T$ and $-V$ rather than $H$). These exact steps preserve the modified hamiltonian function values, but they do not preserve the original value. 

## $T=||\vec p||_2^2/2, \; \nabla_{\vec q}V = -\vec a(\vec q)$

This is the usual kinematric case. The evolution steps are now :

$$
\vec z_T(t+\Delta t) = \vec z(t) + \begin{bmatrix} \vec p(t) \\ \vec 0\end{bmatrix}\Delta t \\
\vec z_V(t+\Delta t) = \vec z_(t) + \begin{bmatrix} \vec 0 \\ \vec a(\vec q(t)) \end{bmatrix}\Delta t
$$

For the position verlet leap-frog method, we do : 

$$
\vec z_{T} = \vec z_n + \begin{bmatrix} \vec p_n \\ \vec 0\end{bmatrix}(\Delta t/2) \\
\vec z_{VT} = \vec z_T + \begin{bmatrix} \vec 0 \\ \vec a(\vec q_{T})\end{bmatrix}\Delta t \\
\vec z_{n+1}= \vec z_{TVT} = \vec z_{VT} + \begin{bmatrix} \vec p_{VT} \\ \vec 0\end{bmatrix}(\Delta t /2)
$$

The same thing, in simpler notation is :

$$
\vec q_{n+1/2} = \vec q_{VT} = \vec q_T =  \vec q_n + \vec p_n (\Delta t/2) \\
\vec p_{n+1} = \vec p_{TVT} = \vec p_{VT} =  \vec p_n + \vec a(\vec q_{n+1/2})\Delta t \\
\vec q_{n+1} = \vec q_{TVT}= \vec q_{n+1/2} + \vec p_{n+1} (\Delta t/2) 
$$

## Simple Harmonic Oscilator

For the potential $V(q) = ||\vec q^2||/2$ , we have $\vec a = - \vec q$ . This then gives us our familiar system of equations :

$$
\frac{d}{dt} \vec z =  \frac{d}{dt} \begin{bmatrix}\vec q \\ \vec p\end{bmatrix}= \begin{bmatrix}\vec p \\ -\vec q\end{bmatrix} = D_H \vec z
$$

Here, we have $D_H = J$, which simplifies a lot of things. 

In any case, the leap-frog method now becomes :

$$
\vec q_{n+1/2}  =  \vec q_n + \vec p_n (\Delta t/2) \\
\vec p_{n+1} =  \vec p_n - \vec q_{n+1/2}\Delta t \\
\vec q_{n+1} =  \vec q_{n+1/2} + \vec p_{n+1} (\Delta t/2) 
$$

Plugging the value of $q_{n+1/2}$ into the the other 2 equations, we get :

$$
\vec p_{n+1} = \vec p_n (1-\Delta t^2/2) -\vec q_n\Delta t \\
\vec q_{n+1} = \vec q_n + \vec p_{n+1}(\Delta t/2) + p_n(\Delta t/2)
$$

Once again, substituting the value of $\vec p_{n+1}$ , we get :

$$
\vec p_{n+1} = \vec p_n (1-\Delta t^2/2) -\vec q_n\Delta t \\
\vec q_{n+1} = \vec q_n (1-\Delta t^2/2) + \vec p_{n}(\Delta t - \Delta t^3/4)
$$

which can be written as :

$$
\vec z_{n+1} =
\begin{bmatrix}
\vec q_{n+1} \\
\vec p_{n+1}
\end{bmatrix}
=
M \vec z_n \\
\text{where } M = \begin{bmatrix}

(1- \Delta t^2/2)I_n & (\Delta t - \Delta t^3/4)I_n \\
-\Delta tI_n &  (1-\Delta t^2/2)I_n

\end{bmatrix}
$$

Now, $M$ is a symplectic operator since it has $M^TJM = J$ . 

## Area preservation

Any symplectic operator (by definition) satisfies $M^T J M = J$ (assuming $M^T$ can be defined) or more generally 

$$
[M \vec u, M \vec v] = [\vec u, \vec v] = (J \vec u)\cdot \vec v = \sum_{i=1}^n \begin{vmatrix} v_i & u_i \\ v_{n+i} & u_{n+i} \end{vmatrix} \;\forall \vec u,\vec v \in \mathbb{R}^{2n}
$$

This is a sum of 2D areas of parallelograms. For the simple case where $n=1$ and we have $u = (\theta_1,\omega_1)$ and $v = (\theta_2,\omega_2)$ , the quatity $[\vec u,\vec v]$ is the signed area of the parallelogram with $\vec u,\vec v$ as edges. The operator (such as $D_H,D_T,D_V,\Psi_t$, etc.), namely $M$, preserves this area in the phase space, if it is symplectic.

While in class, we said that symplectic methods "preserve energy" since they preserve area in phase space, the area enclosed by the trajectory of a system in the phase space is NOT the energy, but rather, the "action". 
One can also verify that the transformation for the step in Euler-Cromer method is area-preserving. But as we saw, there are oscilations in energy $E(t)$. Thus, the Euler-Cromer update step does not preserve energy. 

I am going to be honest and say that I don't understan lagrangian mechanics which is where I last saw the concept of "action" mentioned. 
But from what I read on searching, a constant action means that the process can be reversed. 


## Aharmonic oscilator

Consiider the oscilator described by $\vec z(\tau) = (\theta, \omega)$ that follow these differential equations :

$$
\frac{d\theta}{d\tau} = \omega \\
\frac{d\omega}{d\tau} = a(\theta) = -\sin \theta
$$

From the first equation, we can define $T(\omega) = \frac{1}{2}\omega^2$ and from the second equation, $V(\theta) = 1 -\cos \theta$ . Note that I partly _chose_ these specific functions. One can add constants to these functions and it would still be valid. 

Now, the position Verlet update step is :

$$
\theta_{n+1/2} =  \theta_n + \omega_n \frac{\Delta t}{2} \\
\omega_{n+1} =  \omega_n + a(\theta_{n+1/2})\Delta t \\
\vec q_{n+1} = \theta_{n+1/2} + \omega_{n+1} \frac{\Delta t}{2} 
$$

These are the same equations I derived for the general case with $T = p^2/2$ . In this case, we have $\omega$ and $\theta$ rather than $\vec p$ and $\vec q$.

By plugging in $a(\cdot)$, we get :

$$
\theta_{n+1/2} =  \theta_n + \omega_n \frac{\Delta t}{2} \\
\omega_{n+1} =  \omega_n - \sin(\theta_{n+1/2})\Delta t \\
\vec q_{n+1} = \theta_{n+1/2} + \omega_{n+1} \frac{\Delta t}{2} 
$$

Renaming $\Delta t = h$, we get the same equations as those given in part (b) of problem 3. 

## (g) : Time step as a phase space coordinate transformation

### Part 1

The Euler method time step is the map $F_E : (\theta,\omega) \to (\Theta,\Omega) = (\theta + h \omega, \omega - h\sin\theta)$ .

This gives us the partial derivatives 

$$
\frac{\partial \Theta}{\partial \theta} = 1 \\
\frac{\partial \Theta}{\partial \omega} = h \\
\frac{\partial \Omega}{\partial \theta} = -h\cos\theta \\
\frac{\partial \Omega}{\partial \omega} = 1 \\
$$

Thus, we have the jacobian 

$$
J_E = \begin{bmatrix}
1 & h \\
-h\cos\theta & 1
\end{bmatrix}
$$

This has the determinant

$$
\det J_E (\theta,\omega) = \begin{vmatrix}
1 & h \\
-h\cos\theta & 1
\end{vmatrix} = 1\cdot1 - h\cdot (-h\cos\theta) = 1 - h^2\cos\theta
$$

This means that $\lnot (\det J_E (\theta,\omega) = 1 \forall \theta,\omega)$. 

Since this transformation doesn't necessarily preserve the phase space area locally, we cannot say that it will preserve the are enclosed by the a constant energy trajectory in the phase space. 

Now, since the energy of the oscilator is given by $H(\theta, \omega) = \frac{1}{2}\omega^2 + (1-\cos\theta)$, the energy of the transformed points will be :

$$
H(\Theta,\Omega) = \frac{1}{2}\Omega^2 + (1-\cos\Theta)
$$

Since it will be hard to evaluate this for large $h$, I am opting to compue the derivative with respect to $h$. If that turns out to be non-zero, we'll know that the transformation doesn't preserve energy. 

$$
\frac{d}{dh}H(\Theta,\Omega) = [\Omega(-\sin\theta) + \sin\Theta (\omega)] 
$$

Since clearly, evaluating it at $h=0$ gives us 0, it means that we need to go for the second order derivative. 

Doing that gives us :

$$
\frac{d^2}{dh^2}H(\Theta,\Omega) = [(-\sin\theta)^2 + \cos\Theta (\omega)^2] 
$$

Now evaluating the second derivative at $h = 0$ doesn't give us 0 for all $\theta,\omega$ . What this means is that while the linear term in the taylor expansion of $H(\Theta(h),\Omega(h))$ is 0 and the constant term matches the current energy, since the quadratic term is not 0, thus for large enough $h$, this transformation will alter the energy. 

### Part 2

The Euler Cromer method has the transformation 

$$
\Theta = \theta + h\omega - h^2 \sin\theta \\
\Omega = \omega - h\sin\theta \\ 
$$

This gives us the Jacobian 

$$
J_\text{EC}(\theta,\omega) = 
\begin{bmatrix}
1 - h^2\cos\theta &  h \\
-h\cos\theta & 1
\end{bmatrix}
$$

This has the determinant

$$
\det J_\text{EC} = 1 - h^2\cos\theta + h^2\cos\theta = 1\;\forall \theta,\omega
$$

Thus it preserves area in phase space. 

### Part 3

The determinant being 1 is the discrete analogue for the Lioville operator $D_H$ and the flow function $\Psi_t = \exp(t D_H)$. For the continous evolution, we had 

$$
\vec z(t) = \Psi_t (\vec z_0) = \exp(t D_H) \vec z_0
$$

The flow function had the area preserving property too. Namely, we had $[\Psi_t \vec z_1, \Psi_t \vec z_2] = [\vec z_1,\vec z_2] = (J\vec z_1)\cdot \vec z_2$ . 

Now, rather than $\vec z(t) = \Psi_t \vec z(0)$, we have $\vec z_n = M^n \vec z_0$ where $M$ is the map that the update steps implement. When the algorithm is symplectic, these maps also preserve area. The maps preserving the area in phase space is equivalent to the jacobian of the transformation/map being 1 everywhere (for all $\vec z$).

This also shows in the phase space trajectories of the numerical solutions given by Euler's method and Euler-Cromer (EC) method. The EC method, being symplectic, preserves the area enclosed, whereas with Euler's method the trajectory first starts getting farther from the origina and eventually deviates from what we expect. This is because the area (as we would get for the full, but transformed trajectory) keeps changing. 