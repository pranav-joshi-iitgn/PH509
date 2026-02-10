# Planetary Motion

### Python2 (?) code for basic earth-sun simulation

```python
import ode,numpy as np
import visual as vp

def earth(id, r, v, t):
    global GM
    if id==0 : return v
    s = vp.mag(r) # magnitude
    return - GM*r/s**3

def go():
    r = np.array([1.017, 0.0])
    v = np.array([0.0,6.179]) 
    scene = vp.display(
        title="Planetary motion",
        background=(0.2,0.5,1),
        forward =(0.2,-1),
        planet = vp.sphere(
            pos=r, radius=0.1, 
            make_trail=True, 
            material =vp.materials.earth, 
            up=(0,0,1)
            ),
        sun = vp.sphere(
            pos=(0,0),
            radius = 0.2,
            color=vp.color.yellow,
            material = vp.materials.emissive
            ),
        sunlight = vp.local_light(pos=(0,0),color=vp.color.yellow)
    )
    t = 0.0
    h = 0.001
    while True:
        vp.rate(200)
        r,v = ode.leapdrog(earth, r, v, t, h)
        planet.pos = r
        if (scene.kb.keys): 
            scene.kb.getkey(),scene.kb.getkey()

GM = 4* (np.pi**2)
go()
```

### Terminology

- Aphelion : Position of Earth farthest from Sun
- Perihelion : Position of Earth closes to Sun
- 1 AU (astronomical usint) = $1.496 \times 10^{11}$ m
- 1 year = $3.156 \times 10^7$ s
- 1 AU/year = 4740 m/s
- For the sun, $GM = 4\pi^2 \text{AU}^3\text{year}^{-2}$


### Reduction to one body problem

If we have bodies with masses $m,M$ and positions $\vec r_1 (t), \vec r_2 (t)$ following gravity, then we can rewrite the system as :

$$
\vec R = (\frac{m}{m+M})\vec r_1 + (\frac{M}{m+M})\vec r_2 \quad [\text{CoM}] \\
\vec r = \vec r_1 - \vec r_2 \\
m\ddot{\vec r_1} = \frac{GMm}{r^2} (-\hat r) \quad \text{and} \quad
M\ddot{\vec r_2} = \frac{GMm}{r^2} (\hat r) \\\implies

\ddot{\vec R} = 0 \quad [\text{addition}] \\
\text{and } \ddot{\vec r} = \frac{G(M+m)}{r^2}\hat r
$$

Just solve for $\vec r$ and you get $\vec r_1,\vec r_2$ easily. Since this is exactly the body orbiting-fixed-mass problem, we know the solutions are conic section. 

### Angular Momentum and Effective potential

Define $\omega = \dot{\theta}$ and $v_r = \dot{r}$ as angular and radial velocites. Note that $v_r$ is different from $\vec v = \frac{d}{dt}\vec r$ . 

Define angular momentum as $\vec L = \vec r \times (m \vec v)$ . 

We are assuming movement in x-y plane only, so $\vec L = L \hat k$ . 

Since the torque on the mass $m$ (and the full system for reduced two body case) is 0, we have $\frac{d}{dt}\vec L = \vec 0$ . 

This causes Kepler's second law (equal areas in equal times)

We want our numerical method to keel angular momentum constant too. 

Another thing to keep constant is Energy (Hamiltonian), given by

$$
E(\vec r,\vec v) = \frac{1}{2}mv^2 + V(r)
$$

Here, 

$$
V(\vec r) = -\frac{GMm}{r}
$$

is our usual gravitational potential. 

Now, rewriting things in polar coordinates, 

$$
E(\theta, r, \omega, v_r) = \frac{1}{2}mv_r^2 + \frac{1}{2}m(r\omega)^2 + V(r)
$$

But, due to constant $L$, we have :

$$
mr^2\omega = L \implies r^2\omega^2 = \frac{L^2}{m^2r^2}
$$

which is a function of $\vec r$, and not $\omega$. 

Thus, we have :

$$
E(r, v_r) = \frac{1}{2}mv_r^2 + V_\text{eff}(r) \\
\text{where } V_\text{eff}(r) = \frac{1}{2}\frac{L^2}{mr^2} - \frac{GMm}{r^2}
$$

The $\frac{L^2}{2mr^2}$ term is the "centrifugal potential energy". This arises because we have effectively shiften to a rotating reference frame. 

The equation is now entirely one-dimensional (since $v_r = \dot{r}$) .

Usually, for bounded motion (elliptical trajectories), the values $r_\text{max}$ and $r_\text{min}$ are called "turning points" , which is also where the energy reaches a maximum or minimum value. 

### Kepler Orbits

In polar coordinates, the solution to the two body probem is :

$$
r(\theta) = \frac{L^2}{mk}[\frac{1}{1+e\cos(\theta + \theta_0)}] \\
\text{ where } e = \sqrt{1 + \frac{2EL^2}{mk^2}} \quad [\text{ eccentricity}]
$$

Here $k=GMm$. 

The shape of the orbit is determined as :

|Energy | eccentricity      | shape    |
|-------|-------------------|----------|
| +ve   | > 1               | hyperbola|
| 0     | 1                 | parabola |
| -ve   | < 1               | ellipse  |
| $-\frac{mk^2}{2L^2}$ | 0  | circle   |

For elliptical orbits, we have :

$$
a = \frac{-k}{2E} \\
b = a\sqrt{1-e^2} = \frac{-k}{2E} \sqrt{\frac{-2EL^2}{mk^2}} = \sqrt{-\frac{L^2}{2mE}} = \frac{L}{\sqrt{-2mE}} \\
r(\theta) = \frac{a(1-e^2)}{1-e\cos(\theta + \theta_0)}
$$

Setting $\theta_0 = \pi$ (body starts at perihelion), we get :

$$
r(\theta) = \frac{a(1-e^2)}{1-e\cos\theta}
$$

The tangential velocity components at aphelion ($r_\text{max}$) and perihelion ($r_\text{min}$) are :

$$
v_\text{ap} = \sqrt{\frac{k}{ma} [\frac{1-e}{1+e}]} \\
v_\text{peri} = \sqrt{\frac{k}{ma} [\frac{1+e}{1-e}]}
$$

Since the angulare momentum is fixed, the area is sweeped at a rate of $\frac{L}{2m} = b\sqrt{-2mE}/(2m) = \frac{b\sqrt{-E}}{\sqrt{2m}} = b \sqrt{\frac{k}{4am}}$ . The total area is $\pi a b$ . Thus, the time period is $T = \pi a \sqrt{\frac{4am}{k}} = \sqrt{4\pi^2 a^3 \frac{1}{GM}}$ . This gives Kepler's 3rd law :

$$
T^2 = \frac{4\pi^2 a^3}{GM}
$$

Using this for the sun-earth system, we have :

$$
(1\text{ year})^2 = 4\pi^2 (\text{AU})^3 \frac{1}{GM_s}
$$

Thus, $GM_s = 4\pi^2$ using AU and years. 

### Stable orbits

If by changing the energy by a small amount, the orbits doesn't change qualitatively (is still closed), it's a stable orbit.

A closed orbit has the ratio of number of oscilations in $r$ (namely $n_r (t)$) and oscilations in $\theta$ (namely $n_\theta (t)$) to be a constant ($n_r/n_\theta = \text{constant}$) .

Only forces with a $1/r^2$ or $r$ proportionality can produce closed orbits. Other kind of forces/potentials cause orbital precession. 

### Precesion of platenary orbits

While classic gravity _does_ follow the inverse square law, the planets in the solar system still undergo precession. This has two reasons. 

The first is that other planets also tug and pull at a given planet. This can be simulated numerically simply using newton's theory of gravity. 
For example, doing this for mercury give a precession of 574 arc-seconds per century (1 arc-second is 1/3600 of a degree), or 0.16 degrees per century. This rate is $\frac{d\theta_0}{dt}$ where $\theta_0$ is the aphelion. 

Even with this, there is a leftover precession of 43 arc-seconds, which is due to the second reason. 

The second reason is that, according to general relativity, gravity isn't entirely an inverse-square force. Instead, we have :

$$
F = \frac{GMm}{r^2} + \lambda \frac{GMm}{r^4}
$$ 

### Rung-Lenz vector

We define it as 

$$
\vec A = \vec p \times \vec L - mk\hat r
$$

This is 0 for circular orbits and is a constant vector pointing to the perihelion (negative x axis direction) in case of closed elliptical orbits. 

The rate of change of the angle that $\vec A$ makes with x-axis is the rate of precession. 

### Leapfrog with time transformation

We scale the time step $h$ based on the distance $r$ with a scaling factor $\Omega(r)$. Basically, rather than using $t$ as the path parameter, we use $s(t)$ with the property that $\frac{ds}{dt} = \Omega(r)$ is large when $r$ is small. So now, we need to solve for $\vec v(s), \vec r(s), t(s)$ rather than $\vec v(t), \vec r(t)$ . 

This means :

$$
\frac{d}{ds} \vec r = \vec v / \Omega(r) \\
\frac{d}{ds} \vec v = \vec a(\vec r) /\Omega(r) \\
\frac{d}{ds} t = 1/\Omega(r)
$$


In position-verlet Leapfrog, the way a single step works is through

$$
\vec z_{n+1} = \exp(D_Th/2)\exp(D_V h)\exp(D_Th/2)\vec z_n
$$

and thus we require the time step $h$ to be constant. 

Having time step $\Omega (r_{n}) h$ in the first sub-step and $\Omega (r_{n+1/2}) h$ in other two will make it wrong. 

So, we use the same step scale, call it $W$, which is computed at the very start. So our equations become :

$$
\frac{d}{ds} \vec r = \vec v / W \\
\frac{d}{ds} t = 1/W \\
\frac{d}{ds} \vec v = \vec a(\vec r) /\Omega(r) \\
\frac{d}{ds} W = \frac{d}{ds} \Omega(r) = \nabla \Omega(r) \cdot \frac{d \vec r}{ds} = \nabla \Omega(r) \cdot \vec v / \Omega(r) 
$$

The reason for still writing $\Omega(r)$ is so that this system can be directly put into a leapfrog solver function by thinking of $\vec q = (\vec r, t)$ as position and $\vec  p = (\vec v. W)$ as velocoty vector.

This is not actually necessary (doubt) and doesn't even make physical sense since the $\vec q, \vec p$ don't correspond to a real $T(\vec p),V(\vec q)$ pair. We could just compute $W = \Omega(r)$ directly before doing the leap-frog step and update $t$ separately. This is a more sane approach IMO. 

Their choise of $\Omega(r)$ is $1/r$. 

In their update equations, they don't actually have the same scaling factor for the substeps. I really don't understand what's going on at this point. 

The updates are :

$$
\vec r_{1/2} = \vec r_0 + \vec v_0 \frac{h}{2}\frac{1}{W_0} \\
t_{1/2} = t_0 + \frac{h}{2}\frac{1}{W_0} \\
\vec v_1 = \vec a (\vec r_{1/2}) r_{1/2} h \\
W_1 = W_0 - \vec r_{1/2}\cdot (\frac{\vec v_0 + \vec  v_1}{2})h \frac{1}{r_{1/2}^2} \\
\vec r_1 = \vec r_{1/2} + \vec v_1 \frac{h}{2}\frac{1}{W_1} \\
t_1 = t_{1/2} + \frac{h}{2}\frac{1}{W_1}
$$

### Other planets' effect

To get the net effect of other planet's gravitational feild over centuries, we can do a central field apprroximation. Essentially, over centuries, the other planet has bee at all positions of a (hypothetical) ring with mass distributed over regions in proportion to the time of a "year" that the planet spends in that region. 

We can do something even simpler and just consider a circular orbit, making the ring uniform. 

Then, the non-relativistic gravitational potential for such a planet with orbit radius $a$ when mercury is at distance $r = \rho a$ from sun, will be :

$$
V(r) = - \frac{G m}{a} \frac{M}{2\pi a}\int_{0}^{2\pi} \frac{d\theta}{\sqrt{1 + \rho^2 -2\rho \cos\theta}}
$$

Since $\rho < 1$ for mercury and any other planet, we can (somehow) write this as a power series : 

$$
V(r) = -\frac{GM_p m}{a} [1 + \frac{1}{4}\rho^2 + \frac{9}{64}\rho^4 + \frac{25}{256}\rho^6 + O(\rho^8)]
$$

Then, we take the negative gradient to get the force and then the acceleration :

$$
\vec a_p = \frac{GM_p}{2a^3}\vec r [1 + \frac{9}{8}\rho^2 + \frac{75}{64}\rho^4 + O(\rho^6)]
$$

Note the sign. It's away from the sun.

### Book

- "Computational Physics" by Giordano, Naganishi
- Show that Rung Lenz is constant both computationally and analytically. 
- Show that Rung Lenz vector is related to $SO(4)$