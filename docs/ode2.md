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
def f_sh(theta, omega):
    return omega, -theta
```

```{code-cell} ipython3
def error_RK4(N=100):
    T, X = RK4(f_sh, (0, 0.01), 0, 10/N, N)
    Theta  = X[:,0]
    Thetap = 0.01 * np.sin(T)
    return np.max(abs(Theta - Thetap))

def error_RK38(N=100):
    T, X = RK38(f_sh, (0, 0.01), 0, 10/N, N)
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
def f_dp(th1, th2, p1, p2):
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
T, X = RK4(f_dp, (np.pi/2, np.pi/2, 0.0, 0.0), 0, T/N, N)
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
How do we know if the time step is chosen properly?
How do we know if the solution is accurate enough?

```{code-cell} ipython3
for i, n in enumerate([1000,2000,4000,8000,16000]):
    T, X = RK4(f_dp, (np.pi/2, np.pi/2, 0.0, 0.0), 0, 100/n, n)
    plt.plot(T, X[:,0], '-',  color=f'C{i}', label=f'n={n}')
    plt.plot(T, X[:,1], '--', color=f'C{i}')
plt.legend()
```

## Adaptive Stepsize Control

The double pendulum problem demostrate the challenges posed by highly non-linear and chaotic systems, where even small variations in initial conditions or time step sizes can lead to drastically different outcomes.
This sensitivity highlights the critical importance of selecting an appropriate time step when numerically solving such complex systems.
To determine whether a chosen time step is suitable, one must assess whether it sufficiently captures the rapid changes and intricate dynamics inherent to the system without introducing significant numerical errors.
Additionally, ensuring the solution's accuracy necessitates implementing mechanisms that can evaluate and control these errors effectively.
Adaptive time step control emerges as a vital strategy in this context, as it dynamically adjusts the integration step based on real-time error estimates.
By doing so, it maintains the balance between computational efficiency and the precision required to accurately model the chaotic behavior of systems like the double pendulum, thereby addressing both the proper selection of time steps and the assurance of solution accuracy.

+++

### Error Estimate

It is obvious that there are multiple ways to advance an ODE system with $2\Delta t$:
1. Step the ODE system with a single step $\Delta t' = 2\Delta t$.
2. Step the ODE system with two steps, each step is $\Delta t$.
The error of these two approaches are different.
For a 4th-order algorithm, they are:
\begin{align}
x(t + 2\Delta t) &= x_1 + (2\Delta t)^5 \phi + \mathcal{O}(\Delta t^6) + \dots \\
x(t + 2\Delta t) &= x_2 +  2\Delta t^5  \phi + \mathcal{O}(\Delta t^6) + \dots
\end{align}
where, to order $\Delta t^5$, the value $\phi$ remains constant over the step.

The difference between the two numerical estimates is a convenient indicator of truncation error,
\begin{align}
\Delta = x_2 - x_1.
\end{align}
It is this difference that we shall endeavor to keep to a desired degree of accuracy, neither too large nor too small. We do this by adjusting $\Delta t$.

+++

It might also occur to you that, ignoring terms of order $\Delta t^6$ and higher, we can
solve the two equations to improve our numerical estimate of the true
solution $x(t + 2 \Delta t)$, namely,
\begin{align}
x(t + 2\Delta t) = x_2 + \frac{\Delta}{15} + \mathcal{O}(\Delta t^6)
\end{align}
This estimate is accurate to fifth order, one order higher than the original Runge-Kutta steps.
This is nothing but Richardson extrapolation.

However, we can't have our cake and eat it too.
The above equation may be fifth-order accurate, but we have no way of monitoring its truncation error. Higher order is not always higher accuracy!

+++

### Embedded Runge-Kutta Formulas: The Dormand–Prince Method

Step doubling is a simple adaptive step method.
However, it is now mostly replaced by a more efficient stepsize adjustment algorithm based on embedded Runge-Kutta formulas, originally invented by Merson and popularized in a method of Fehlberg.

An interesting fact about Runge-Kutta formulas is that for orders $M$ higher than four, more than $M$ function evaluations are required.
This accounts for the popularity of the classical fourth-order method: It seems to give the most bang for the buck. However, Fehlberg discovered a fifth-order method with six function evaluations where another combination of the six functions gives a fourth-order method.

The difference between the two estimates of $y(x + \Delta t)$ can then be used as an estimate of the truncation error to adjust the stepsize.
Since Fehlberg's original formula, many other embedded Runge-Kutta formulas have been found.

+++

The Dormand–Prince method is a 5th-order Runge-Kutta method with an embedded 4th-order solution.
This pairing enables accurate error estimation, which is essential for adaptive step size control.
The DP method is particularly favored for its robustness and efficiency, making it a staple in many numerical computing libraries, including MATLAB's ode45.

The Dormand–Prince method is characterized by its specific set of coefficients, which dictate how intermediate slopes are calculated and combined to produce the final solution estimates.
These coefficients are meticulously chosen to minimize error and optimize computational performance.

+++

The general form of a fifth-order Runge-Kutta formula is
\begin{align}
k_1 &= \Delta t f(x_n, t_n)\\
k_2 &= \Delta t f(x_n + a_{21}k_1, t_n + c_2 \Delta t)\\
\cdots\\
k_6 &= \Delta t f(x_n + a_{61}k_1 + \cdots + a_{65}k_5, t_n + c_6 \Delta t)\\
x_{n+1} &= x_n + b_1 k_1 + b_2 k_2 + \cdots + b_6 k_6 + \mathcal{O}(\Delta t^6)
\end{align}

The Dormand–Prince method employs a set of coefficients $a$, $b$, and $c$ to compute intermediate stages and combine them into higher and lower-order solutions.
Below are the coefficients for the DP method:

```{code-cell} ipython3
a = [
    [],
    [1/5],
    [3/40, 9/40],
    [44/45, -56/15, 32/9],
    [19372/6561, -25360/2187, 64448/6561, -212/729],
    [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656],
    [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84],
]
b_high = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0] # Fifth-order accurate solution estimate
b_low  = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40] # Fourth-order accurate solution estimate
c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]
```

A DP45 step is therefore

```{code-cell} ipython3
def DP45_step(f, x, t, dt):
        # Compute intermediate k1 to k7
    k1 = np.array(f(*x))
    k2 = np.array(f(*(x + dt*(a[1][0]*k1))))
    k3 = np.array(f(*(x + dt*(a[2][0]*k1 + a[2][1]*k2))))
    k4 = np.array(f(*(x + dt*(a[3][0]*k1 + a[3][1]*k2 + a[3][2]*k3))))
    k5 = np.array(f(*(x + dt*(a[4][0]*k1 + a[4][1]*k2 + a[4][2]*k3 + a[4][3]*k4))))
    k6 = np.array(f(*(x + dt*(a[5][0]*k1 + a[5][1]*k2 + a[5][2]*k3 + a[5][3]*k4 + a[5][4]*k5))))
    k7 = np.array(f(*(x + dt*(a[6][0]*k1 + a[6][1]*k2 + a[6][2]*k3 + a[6][3]*k4 + a[6][4]*k5 + a[6][5]*k6))))

    ks = [k1, k2, k3, k4, k5, k6, k7]

    # Compute high and low order estimates
    x_high = x + dt * np.dot(b_high, ks)
    x_low  = x + dt * np.dot(b_low,  ks)

    return x_high, x_low, ks
```

### Proportional-Integral Step Size Control

Once we have an embedded Runge-Kutta method like Dormand–Prince in place, the next step is to implement a mechanism to adjust the step size based on the estimated local error.
The PI (Proportional-Integral) controller is a widely-used strategy for this purpose, combining proportional and integral components to achieve stable and efficient step size adjustments.

The PI controller adjusts the step size $\Delta t$ based on the current error $E$ and past errors.
The objective is to maintain the error within specified tolerances while preventing drastic changes in step size that could lead to instability or inefficiency.

The general formula for updating the step size is:
\begin{align}
h_{\text{new}} = h \cdot \min\left(\text{fac}{\text{max}}, \max\left(\text{fac}{\text{min}}, \text{fac} \cdot \left(\frac{\text{tol}}{E}\right)^{\alpha}\right)\right)
\end{align}
where $\text{fac}$ is a scaling factor (typically around 0.9) to provide a safety margin;
$\text{fac}{\text{min}}$ and $\text{fac}{\text{max}}$ set the minimum and maximum allowable step size multipliers to prevent excessive changes; and
$\alpha$ is an exponent that determines the responsiveness of the step size adjustment.

```{code-cell} ipython3
def dt_update(dt, error, tol, fac=0.9, fac_min=0.1, fac_max=4.0, alpha=0.2):
    if error == 0:
        s = fac_max
    else:
        s = fac * (tol / error) ** alpha
    s = min(fac_max, max(fac_min, s))
    dt_new = dt * s
    return dt_new
```

Combining the single embedded step and the step controller, we obtain the DP45 algorithm:

```{code-cell} ipython3
def DP45(f, x, t, T, dt, atol, rtol):

    Ts = [t]
    Xs = [np.array(x)]

    while t < T:
        if t + dt > T:
            dt = T - t  # Adjust step size to end exactly at tf

        # Perform a single Dormand–Prince step
        x_high, x_low, _ = DP45_step(f, x, t, dt)

        # Compute the error estimate
        error = np.linalg.norm(x_high - x_low, ord=np.inf)

        # Compute the tolerance
        tol = atol + rtol * np.linalg.norm(x_high, ord=np.inf)

        # Check if the step is acceptable
        if error <= tol:
            # Accept the step
            t += dt
            x = x_high
            Ts.append(t)
            Xs.append(x)

        # Compute the new step size
        dt = dt_update(dt, error, tol)

    return np.array(Ts), np.array(Xs)
```

We can now apply it to the ODEs.

Applying to the simple harmonic oscillator, we may specifically ask for the accuracy:

```{code-cell} ipython3
def error_DP45(tol):
    T, X = DP45(f_sh, (0, 0.01), 0, 10, 0.1, tol, tol)
    Theta  = X[:,0]
    Thetap = 0.01 * np.sin(T)
    return np.max(abs(Theta - Thetap))

N     = np.array([64, 128, 256, 512, 1024])
EDP45 = np.array([error_DP45(tol) for tol in 2.0**(-22-4*np.arange(5))])

plt.loglog(N, 1/N**4,      label='1/N^4')
plt.loglog(N, ERK4,  'o-', label='RK4')
plt.loglog(N, ERK38, 'o:', label='RK38')
plt.loglog(N, EDP45, 'o--', label='RK38')
plt.xlabel('N')
plt.ylabel(r'$\text{err} = max|x_\text{numeric} - x|$')
plt.legend()
```

For non-linear problems, we can compare different accuracies:

```{code-cell} ipython3
for i, atol in enumerate([1e-3,1e-6,1e-9,1e-12]):
    T, X = DP45(f_dp, (np.pi/2, np.pi/2, 0.0, 0.0), 0, 10, 0.1, atol, 0)
    plt.plot(T[::10], X[::10,0], '.-',  color=f'C{i}', label=f'atol={atol}')
    plt.plot(T[::10], X[::10,1], '.--', color=f'C{i}')
plt.legend()
```

Up to this point, we have explored the advanced concepts of Adaptive Step Size Control in Runge-Kutta methods, focusing on the Dormand–Prince (DP) method.
We began by understanding the significance of embedded Runge-Kutta formulas, which enable simultaneous computation of solutions of different orders for effective error estimation.
This foundation allowed us to implement a PI controller that dynamically adjusts the integration step size based on local error estimates, ensuring that the numerical solution remains within desired accuracy bounds while optimizing computational efficiency.
