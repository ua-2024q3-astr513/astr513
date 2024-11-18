---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Integration of Ordinary Differential Equations. II

+++

## Derivation of the Classical 4th-order Runge-Kutta Method

In the last lecture, we learned the very popular and robust classical 4th-order Runge-Kutta method.
Here, we will try to derive it.
Recalling we want to solve an ODE:
\begin{align}
\frac{dx}{dt} = f(x)
\end{align}
The solution $x(t)$ evaluated at $t_{n+1} = (n+1) \Delta t$ can be written in terms of $x(t)$ evaluated at $t_n$ in Taylor series.
I.e.,
\begin{align}
x_{n+1} = x_n + x'_n \Delta t + \frac{1}{2}x''_n \Delta t + \frac{1}{3!} x'''_n \Delta t^3 + \frac{1}{4!} x''''_n \Delta t^4 + \cdots
\end{align}
Here, we use $x_n \equiv x(t_n) \equiv x(n \Delta t)$.
Substituting the ODE into the Taylor series, we obtain
\begin{align}
x_{n+1} = x_n + f_n \Delta t + \frac{1}{2}f'_n \Delta t + \frac{1}{3!} f''_n \Delta t^3 + \frac{1}{4!} f'''_n \Delta t^4 + \cdots
\end{align}

+++

To construct Runge-Kutta method, we consider a formulation
\begin{align}
x_{n+1} = x_n + a_1 \Delta_1 x_n + a_2 \Delta_2 x_n + \cdots + a_s \Delta_s x_n
\end{align}
for some $s$, where
\begin{align}
\Delta_1 x_n &\equiv f(x(t_n)) \Delta t \\
\Delta_2 x_n &\equiv f(x(t_n + b_2 \Delta t)) \Delta t = (f_n + f'_n b_2 \Delta t + \frac{1}{2} f''_n b_2^2 \Delta t^2 + \cdots)\Delta t\\
\cdots \\
\Delta_s x_n &\equiv f(x(t_n + b_s \Delta t)) \Delta t = (f_n + f'_n b_s \Delta t + \frac{1}{2} f''_n b_s^2 \Delta t^2 + \cdots)\Delta t
\end{align}

+++

Substitute, we obtain
\begin{align}
x_{n+1} = x_n 
&+ a_1 f_n \Delta t \\
&+ a_2 (f_n \Delta t + f'_n b_2 \Delta t^2 + \frac{1}{2} f''_n b_2^2 \Delta t^3 + \cdots) \\
&+ \cdots \\
&+ a_s (f_n \Delta t + f'_n b_s \Delta t^2 + \frac{1}{2} f''_n b_s^2 \Delta t^3 + \cdots)
\end{align}

For 4th-order scheme, collecting the terms and require all terms up to $\Delta t^4$ match, we obtain the conditions
\begin{align}
a_1       + a_2       + a_3       + a_4       &=       1     \\
            a_2 b_2   + a_3 b_3   + a_4 b_4   &= \frac{1}{2} \\
            a_2 b_2^2 + a_3 b_3^2 + a_4 b_4^2 &= \frac{1}{3} \\
            a_2 b_2^3 + a_3 b_3^3 + a_4 b_4^3 &= \frac{1}{4}
\end{align}

+++

For the classical 4th-order Runge-Kutta scheme, we have already decided $b_1 = 0$, $b_2 = b_3 = 1/2$, and $b_4 = 1$.
Therefore, the system of coefficients read:
\begin{align}
a_1       + a_2       + a_3       + a_4 &=       1     \\
\frac{1}{2} a_2 + \frac{1}{2} a_3 + a_4 &= \frac{1}{2} \\
\frac{1}{4} a_2 + \frac{1}{4} a_3 + a_4 &= \frac{1}{3} \\
\frac{1}{8} a_2 + \frac{1}{8} a_3 + a_4 &= \frac{1}{4}
\end{align}
It is then easy to verify that 
\begin{align}
(a_1, a_2, a_3, a_4) = \left(\frac{1}{6},\frac{1}{3},\frac{1}{3},\frac{1}{6}\right)
\end{align}
is the solution.

+++

Note that, if one choose $b_1 = 0$, $b_2 = 1/3$, $b_3 = 2/3$, and $b_4 = 1$, then the solutoin is 
\begin{align}
(a_1, a_2, a_3, a_4) = \left(\frac{1}{8},\frac{3}{8},\frac{3}{8},\frac{1}{8}\right).
\end{align}
This is Wilhelm Kutta (1901)'s "3/8 method".

+++

This suggests that Runge-Kutta methods are really a "family", where many different choices can be used to construct numerical schemes with the same order.
The perform of the numerical scheme, nevertheless, depends on the number of oeprations as well as the properties of the ODEs being solved.

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
```

```{code-cell} ipython3
def RK4(f, x, t, dt, n):
    T = np.array(t)
    X = np.array(x)
    
    for i in range(n):
        k1 = dt * np.array(f(*(x         )))
        k2 = dt * np.array(f(*(x + 0.5*k1)))
        k3 = dt * np.array(f(*(x + 0.5*k2)))
        k4 = dt * np.array(f(*(x +     k3)))
        
        t += dt
        x += k1/6 + k2/3 + k3/3 + k4/6
        
        T = np.append( T, t)
        X = np.vstack((X, x))
        
    return T, X
```

```{code-cell} ipython3
def RK38(f, x, t, dt, n):
    T = np.array(t)
    X = np.array(x)
    
    for i in range(n):
        k1 = dt * np.array(f(*(x             )))
        k2 = dt * np.array(f(*(x +       k1/3)))
        k3 = dt * np.array(f(*(x +    k2-k1/3)))
        k4 = dt * np.array(f(*(x + k3-k2+k1  )))
        
        t += dt
        x += (k1 + 3*k2 + 3*k3 + k4)/8
        
        T = np.append( T, t)
        X = np.vstack((X, x))
        
    return T, X
```

```{code-cell} ipython3
def f(theta, omega):
    return omega, -theta
```

```{code-cell} ipython3
def error_RK4(N=100):
    T, X = RK4(f, (0, 0.01), 0, 10/N, N)
    Theta  = X[:,0]
    Thetap = 0.01 * np.sin(T)
    return np.max(abs(Theta - Thetap))

def error_RK38(N=100):
    T, X = RK38(f, (0, 0.01), 0, 10/N, N)
    Theta  = X[:,0]
    Thetap = 0.01 * np.sin(T)
    return np.max(abs(Theta - Thetap))

N     = np.array([64, 128, 256, 512, 1024])
ERK4  = np.array([error_RK4(n)  for n in N])
ERK38 = np.array([error_RK38(n) for n in N])

plt.loglog(N, 1/N**4,      label='1/N^4')
plt.loglog(N, ERK4,  'o-', label='RK4')
plt.loglog(N, ERK38, 'o:', label='RK38')
plt.xlabel('N')
plt.ylabel(r'$\text{err} = max|x_\text{numeric} - x|$')
plt.legend()
```

## The Double Pendulum Problem

The double pendulum is a well known example of a non-linear, chaotic system in classical mechanics.
It consists of a pendulum with another pendulum attached to its end, resulting in a system with two degrees of freedom.
This configuration leads to highly complex and sensitive-dependent dynamics, making the double pendulum an excellent subject for studying chaos theory and non-linear dynamics.
Because it is not possible to construct analytical solutions, it is also a great example to numerical integrators.

![Double pendulum](https://upload.wikimedia.org/wikipedia/commons/c/c9/Double-compound-pendulum-dimensioned.svg)

+++

To setup the equations of motion, we assume:

* The two arms of the pendulums have the same length $l$.

* The mass of each arm is $m$.

* The angle between the first and second pendulums, with respect to the vertical axis, are denoted by $\theta_1$ and $\theta_2$.

Newton's second law suggests that we will need to solve a system of two second-order ordinary differential equations (ODEs).
Using the methods we learn in the lecture, we can cast the problem into a system of four first-order ODEs.
\begin{align}
\frac{d\theta_1}{dt} &=
\frac{6}{m l^2}\frac{2 p_1 - 3 \cos(\theta_1 - \theta_2) p_2}{16 - 9 \cos^2(\theta_1 - \theta_2)}\\
\frac{d\theta_2}{dt} &=
\frac{6}{m l^2}\frac{8 p_2 - 3 \cos(\theta_1 - \theta_2) p_1}{16 - 9 \cos^2(\theta_1 - \theta_2)}\\
\frac{dp_1}{dt} &=
-\frac{1}{2} m l^2 \left(\frac{d\theta_1}{dt} \frac{d\theta_2}{dt}\sin(\theta_1 - \theta_2) +
                           3\frac{g}{l}\sin\theta_1\right)\\
\frac{dp_2}{dt} &=
-\frac{1}{2} m l^2 \left(-\frac{d\theta_1}{dt} \frac{d\theta_2}{dt}\sin(\theta_1 - \theta_2) +
                            \frac{g}{l}\sin\theta_2\right)
\end{align}
where $p_1$ and $p_2$ are called the generalized momenta.
(There might be typos in the equation.
Please [double check](https://en.wikipedia.org/wiki/Double_pendulum).)

```{code-cell} ipython3
def f(th1, th2, p1, p2):
    m  = 1
    l  = 1
    g  = 1
    
    u1 = m * l * l
    u2 = g / l
    f  = 6 / (u1 * (16 - 9 * np.cos(th1 - th2)**2))
    
    dth1 = f * (2 * p1 - 3 * np.cos(th1 - th2) * p2)
    dth2 = f * (8 * p2 - 3 * np.cos(th1 - th2) * p1)
    
    dp1  = - 0.5 * u1 * (  dth1 * dth2 * np.sin(th1 - th2) + 3 * u2 * np.sin(th1))
    dp2  = - 0.5 * u1 * (- dth1 * dth2 * np.sin(th1 - th2) +     u2 * np.sin(th2))
    
    return dth1, dth2, dp1, dp2
```

```{code-cell} ipython3
T = 100
N = 1000
T, X = RK4(f, (np.pi/2, np.pi/2, 0.0, 0.0), 0, T/N, N)
```

```{code-cell} ipython3
plt.plot(T, X[:,0])
plt.plot(T, X[:,1])
```

We may improve the visualization and create a movie:

```{code-cell} ipython3
%%capture

# Step 6. Improve the visualization

from matplotlib import animation
from IPython.display import HTML

fig = plt.figure(figsize=(8,8))
ax  = plt.axes(xlim=(-2.5, 2.5), ylim=(-2.5, 2.5))
ax.set_aspect('equal')

line, = ax.plot([], [], 'o-', lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    th1 = X[i,0]
    th2 = X[i,1]
    
    x1 =   np.sin(th1)
    y1 = - np.cos(th1)
    
    x2 =   np.sin(th2)
    y2 = - np.cos(th2)
    
    line.set_data([0, x1, x1+x2], [0, y1, y1+y2])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=N, interval=20, blit=True)
```

```{code-cell} ipython3
# Uncommon the following when run locally 
#HTML(anim.to_html5_video())
```

```{code-cell} ipython3
# Uncommon the following when run locally 
#anim.save('double_pendulum.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
```

However, because the problem is so non-linear, the solution changes significantly not just when we change the initial conditions, but also change the time step.
How do we know if the time step is chosen propertly?
How do we know if the solution is accurate enough?

```{code-cell} ipython3
for i, n in enumerate([1000,2000,4000,8000,16000]):
    T, X = RK4(f, (np.pi/2, np.pi/2, 0.0, 0.0), 0, 100/n, n)
    plt.plot(T, X[:,0], '-',  color=f'C{i}', label=f'n={n}')
    plt.plot(T, X[:,1], '--', color=f'C{i}')
plt.legend()
```

## Numerical Stability of Integrators

Numerical Stability in the context of ODE solvers refers to the ability of a numerical method to control the growth of errors introduced during the iterative process of approximation.
A method is stable if the numerical errors do not amplify uncontrollably as computations proceed.
* Definition: A numerical method is stable for a given problem if the errors (from truncation or round-off) do not grow exponentially with time steps.
* Importance: Stability ensures that the numerical solution behaves similarly to the true solution, especially over long integration intervals.

Stability is a different concept than accuracy.
* Stability: Concerns the boundedness of errors over time.
* Accuracy: Pertains to how closely the numerical solution approximates the exact solution.

A method can be stable but not necessarily accurate, and vice versa.
However, both properties are requried for reliable numerical solutions.

+++

### Stability Analysis Using the Linear Test Equation

To analyze the stability of numerical integrators, we commonly use a linear test equation from last lecture:
\begin{align}
\frac{dx}{dt} = \lambda x
\end{align}
where $\lambda \in \mathbb{C}$.
The exact solution is:
\begin{align}
x(t) = x_0 e^{\lambda t}.
\end{align}

+++

Consider a one-step numerical method applied to the test equation.
The update can generally be expressed as:
\begin{align}
x_{n+1} = R(z) x_n
\end{align}
where $R(z)$ is the amplification factor and $z = \lambda \Delta t$ is the stability parameter.

For the numerical method to be stable, the magnitude of the amplification factor must satisfy:
\begin{align}
|R(z)| \leq 1
\end{align}
This **stability condition** ensures that errors do not grow exponentially with each step.

+++

In [ODE I](ode1.md), we introduced the Forward Euler method as the simplest explicit numerical integrator for solving ODEs.
Let's revisit its stability properties using the linear test equation.

The forward Euler update formula can be rewritten as:
\begin{align}
x_{n+1} = x_n + \Delta t \cdot f(x_n, t_n) = (1 + \lambda\Delta t) x_n.
\end{align}

The amplification factor is therefore
\begin{align}
R(z) = 1 + \lambda\Delta t.
\end{align}
The stability condition is
\begin{align}
|R(z)| = |1 + \lambda\Delta t| \leq 1.
\end{align}

Graphically, the stability region is a circle in the complex plane centered at $(-1, 0)$ with a radius of 1, as shown in the following figure.

```{code-cell} ipython3
# Define the grid for the complex plane
Re = np.linspace(-3, 3, 601)
Im = np.linspace(-2, 2, 401)
Re, Im = np.meshgrid(Re, Im)
Z = Re + 1j * Im

# Forward Euler stability condition |1 + Z| <= 1
abs_R_fE = np.abs(1 + Z)

# Plotting
plt.figure(figsize=(8, 6))
plt.contourf(Re, Im, abs_R_fE, levels=[0, 1], colors=['red'])

plt.title('Stability Regions for Euler Method')
plt.xlabel(r'Re($\lambda \Delta t$)')
plt.ylabel(r'Im($\lambda \Delta t$)')
plt.gca().set_aspect('equal')

## Adding legend manually
import matplotlib.patches as mpatches
blue_patch = mpatches.Patch(color='red', label='Forward Euler')
plt.legend(handles=[blue_patch])
```

Assuming $\lambda > 0 \in \mathbb{R}$.
As the forward Euler method requires $\Delta t > 0$, the forward Euler is **unconditionally unstable**.
Although we have used forward Euler to solve $dx/dt = \lambda x$ earlier, the error of the solution is unbounded and forward Euler is actually not useful in solving this problem!

+++

### Stability Analysis Using Simple Harmonic Oscillator

From the previous lecture [ODE I](ode1.md), we also solved the simple harmonic oscillator:
\begin{align}
\frac{d}{dt}\begin{bmatrix} \theta(t) \\ \Omega(t) \end{bmatrix} = 
\begin{bmatrix} \Omega(t) \\ -\frac{g}{l} \theta(t) \end{bmatrix} = 
\begin{bmatrix} 0 & 1 \\ -\frac{g}{l} & 0 \end{bmatrix} 
\begin{bmatrix} \theta(t) \\ \Omega(t) \end{bmatrix}
\end{align}
using the forward Euler method:
\begin{align}
\begin{bmatrix} \theta_{n+1} \\ \Omega_{n+1} \end{bmatrix} = 
\begin{bmatrix} \theta_{n} \\ \Omega_{n} \end{bmatrix} +
\begin{bmatrix} 0 & 1 \\ -\frac{g}{l} & 0 \end{bmatrix} 
\begin{bmatrix} \theta_{n} \\ \Omega_{n} \end{bmatrix} = 
\begin{bmatrix} 1 & \Delta t \\ -\frac{g}{l}\Delta t & 1 \end{bmatrix} 
\begin{bmatrix} \theta_{n} \\ \Omega_{n} \end{bmatrix}
\end{align}
The "amplification factor" is no longer a scalar but a matrix.
The stability condition $|R| \le 1$ become
\begin{align}
\det \begin{bmatrix} 1 & \Delta t \\ -\frac{g}{l}\Delta t & 1 \end{bmatrix} \le 1.
\end{align}
Hence,
\begin{align}
\frac{g}{l}\Delta t^2 \le 0
\end{align}
which canno be satisfied.
The forward Euler method is therefore again unconditional unstable.

+++

Similar, we may compute the stability regions for different numerical schemes.
Given a numerical scheme with order $n$ should agree with the analytical solution in its Taylor series up to order $n$, we have:

```{code-cell} ipython3
# Define the grid for the complex plane
Re = np.linspace(-5, 3, 600)  # Real axis from -5 to 3
Im = np.linspace(-3, 3, 600)  # Imaginary axis from -3 to 3
Re, Im = np.meshgrid(Re, Im)
Z = Re + 1j * Im  # Complex grid

def R_fE(z):
    """Stability function for forward Euler"""
    return 1 + z

def R_bE(z):
    """Stability function for backward Euler"""
    return 1 - z

def R_RK2(z):
    """Stability function for RK2 (Heun's method)"""
    return 1 + z + 0.5 * z**2

def R_RK4(z):
    """Stability function for RK4 (classical Runge-Kutta)"""
    return 1 + z + 0.5 * z**2 + (1/6) * z**3 + (1/24) * z**4

# Compute |R(z)| for each method
abs_R_fE  = np.abs(R_fE(Z))
abs_R_bE  = np.abs(R_bE(Z))
abs_R_RK2 = np.abs(R_RK2(Z))
abs_R_RK4 = np.abs(R_RK4(Z))

# Define a list of methods and their corresponding data
methods = {
    'forward Euler' :abs_R_fE,
    'backward Euler':abs_R_bE,
    'RK2'           :abs_R_RK2,
    'RK4'           :abs_R_RK4,
}

plt.figure(figsize=(8, 6))

# Plot stability regions for each method
patches = []
for c, (title, abs_R) in enumerate(list(methods.items())[::-1]):
    # Contour where |R(z)| = 1
    plt.contourf(Re, Im, abs_R, levels=[0, 1], colors=[f'C{c}'])
    patches.append(mpatches.Patch(color=f'C{c}', label=title))

plt.xlabel(r'Re($\lambda \Delta t$)')
plt.ylabel(r'Im($\lambda \Delta t$)')
plt.gca().set_aspect('equal')

## Adding legend manually
plt.legend(handles=patches)
```
