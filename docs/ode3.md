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

# Integration of Ordinary Differential Equations. III

+++

## Numerical Stability of Integrators

Numerical Stability in the context of Ordinary Differential Equation (ODE) solvers refers to the ability of a numerical method to control the growth of errors introduced during the iterative process of approximation.
A method is considered stable if the numerical errors do not amplify uncontrollably as computations proceed.
This concept is crucial for ensuring that the numerical solution remains meaningful and behaves similarly to the true solution, especially over long integration intervals.

+++

### Definition and Importance

A numerical method is stable for a given problem if the errors---whether from truncation or round-off---do not grow exponentially with each time step.
Stability ensures that the numerical solution does not diverge from the true solution due to the accumulation of numerical errors.
This is particularly important in long-term simulations where even small errors can accumulate to produce significant deviations from the true behavior of the system.

+++

### Stability vs. Accuracy

It's essential to distinguish between stability and accuracy:
* Stability: Pertains to the boundedness of errors over time.
  A stable method ensures that errors remain controlled and do not grow exponentially, preventing the solution from becoming meaningless.
* Accuracy: Refers to how closely the numerical solution approximates the exact solution.
  An accurate method minimizes the difference between the numerical and true solutions.

A method can be stable but not accurate, meaning it controls error growth but doesn't closely follow the exact solution.
Conversely, a method can be accurate but unstable, producing precise results initially but diverging over time due to uncontrolled error growth.
For reliable numerical solutions, both stability and accuracy are required.

+++

### Stability Regions for Explicit Methods

To analyze the stability of numerical integrators, we commonly use the linear test equation introduced in previous lectures:
\begin{align}
\frac{dx}{dt} = \lambda x
\end{align}
where $\lambda \in \mathbb{C}$.
The exact solution to this ODE is:
\begin{align}
x(t) = x_0 e^{\lambda t}.
\end{align}
For a numerical method to be stable when applied to this equation, it must ensure that the numerical errors do not grow exponentially.
This is evaluated using the concept of the **stability region**.

+++

### Forward Euler Method Stability

In [ODE I](ode1.md), we learned that the Forward Euler method is the simplest explicit numerical integrator for solving ODEs.
Its update formula for the linear test equation is derived as follows:

Starting with the Forward Euler update:
\begin{align}
x_{n+1} = x_n + \Delta t \cdot f(x_n, t_n)
\end{align}
For the linear test equation $f(x, t) = \lambda x$, this becomes:
\begin{align}
x_{n+1} = x_n + \Delta t \cdot \lambda x_n = (1 + \lambda \Delta t) x_n
\end{align}
Here, the amplification factor $R(z)$ is defined as:
\begin{align}
R(z) = 1 + z \quad \text{where} \quad z = \lambda \Delta t
\end{align}
The stability condition requires that:
\begin{align}
|R(z)| \leq 1
\end{align}
Substituting the amplification factor:
\begin{align}
|1 + z| \leq 1
\end{align}
This condition defines the stability region for the Forward Euler method.

+++

To better understand the stability regions, we can visualize them in the complex plane.
The stability region for Forward Euler is a circle centered at $(-1, 0)$ with a radius of 1.

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
```

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

plt.title('Stability Region for Forward Euler Method')
plt.xlabel(r'Re($\lambda \Delta t$)')
plt.ylabel(r'Im($\lambda \Delta t$)')
plt.gca().set_aspect('equal')

# Adding legend manually
import matplotlib.patches as mpatches
blue_patch = mpatches.Patch(color='red', label='Forward Euler')
plt.legend(handles=[blue_patch])
```

The stability region plot leads to:
* The Forward Euler method is conditionally stable, meaning it is only stable within a specific region of the complex plane.
* Specifically, it is stable for values of $\lambda \Delta t$ that lie within the circle centered at $(-1, 0)$ with a radius of 1.
* For real positive $\lambda$ (i.e., $\lambda > 0$), the Forward Euler method becomes **unconditionally unstable** because such values fall outside the stability region.
  This implies that even with small time steps, the method cannot control error growth for these cases.

+++

The implication is:
* When dealing with ODEs where $\lambda > 0$, especially in stiff equations, the Forward Euler method may lead to solutions that diverge from the true behavior due to uncontrolled error amplification.
* This limitation necessitates the use of more stable methods, such as implicit integrators, which will be discussed in subsequent sections.

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
