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
import numpy as np
import matplotlib.pyplot as plt

# Define the grid for the complex plane
Re = np.linspace(-3, 3, 601)
Im = np.linspace(-2, 2, 401)
Re, Im = np.meshgrid(Re, Im)
Z = Re + 1j * Im

# Forward Euler stability condition |1 + Z| <= 1
Euler_stable = np.abs(1 + Z) <= 1

# Plotting
plt.figure(figsize=(8, 6))
plt.contour(Re, Im, Euler_stable , levels=[0.5],
            colors='blue', linestyles='-')

plt.title('Stability Regions for Euler Method')
plt.xlabel(r'Re($\lambda \Delta t$)')
plt.ylabel(r'Im($\lambda \Delta t$)')
plt.gca().set_aspect('equal')

## Adding legend manually
import matplotlib.patches as mpatches
blue_patch = mpatches.Patch(color='blue', label='Forward Euler')
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
