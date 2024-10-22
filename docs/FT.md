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

# Fourier Transform and Spectral Analyses

## Introduction

The Fourier Transform is a fundamental tool in computational astrophysics and many other fields of science and engineering.
It allows for the decomposition of complex functions or signals into sums of sinusoidal components, facilitating the analysis of their frequency content.
This capability is essential for understanding various phenomena, from the distribution of matter in the universe to the behavior of signals received from distant astronomical objects.

The development of the Fast Fourier Transform (FFT) algorithm significantly advanced computational methods in the 20th century.
The FFT reduces the computational complexity of performing Fourier Transforms from $\mathcal{O}(N^2)$ to $\mathcal{O}(N \log N)$, enabling efficient processing of large datasets typical in astronomical observations and simulations.
Its applications include:
* Communication Systems:
  * Signal Processing: Fundamental in digital communication for encoding and decoding signals, allowing efficient data transmission over various frequencies without interference.
  * Compression: Techniques like MP3 (audio) and JPEG (image) rely on Fourier Transform methods to compress data by transforming signals into frequency components and discarding less significant ones.
* Radar and Sonar:
  * Object Detection: Helps interpret reflected electromagnetic or sound waves to determine the distance, speed, and characteristics of objects, crucial for navigation, meteorology, and defense.
  * Doppler Shifts: Analyzing frequency shifts in returned signals enables calculation of the velocity of moving objects.
* Astronomy:
  * Spectroscopy: Used to analyze light spectra from celestial objects, determining their composition, temperature, velocity, and other properties.
  * Imaging: Employed in reconstructing images from interferometric data, enhancing resolution beyond the limitations of individual telescopes.
  * Very Long Baseline Interferometry (VLBI): Combines data from multiple telescopes across vast distances to simulate a telescope the size of Earth, relying on Fourier analysis to synthesize high-resolution images of distant astronomical sources.

+++

### Historical Context and the Heat Equation

The origins of Fourier analysis trace back to the early 19th century with the work of Jean-Baptiste Joseph Fourier (1768–1830), a French mathematician and physicist.
While studying heat flow, Fourier introduced the revolutionary idea that any periodic function could be expressed as an infinite sum of sine and cosine functions—known today as the Fourier Series.

The one-dimensional heat equation is:
\begin{align}
\frac{\partial u(x,t)}{\partial t} = \alpha \frac{\partial^2 u(x,t)}{\partial x^2},
\end{align}
where $u(x,t)$ is the temperature distribution along a rod at position $x$ and time $t$, and
$\alpha$ is the thermal diffusivity constant of the material.

+++

### Solution Using Separation of Variables

Here, we use the method of separation of variables, which is taught in undergraduate level Electricity and Magnetism, to motivate the Fourier series.
Assuming the solution can be written as a product of functions, each depending only on a single variable:
each depending only on a single variable:
\begin{align}
u(x,t) = X(x) T(t).
\end{align}
Substituting into the heat equation:
\begin{align}
X(x) \frac{dT(t)}{dt} = \alpha T(t) \frac{d^2 X(x)}{dx^2}.
\end{align}
Dividing both sides by $X(x) T(t)$, we obtain
\begin{align}
\frac{1}{T(t)} \frac{dT(t)}{dt} = \alpha \frac{1}{X(x)} \frac{d^2 X(x)}{dx^2}.
\end{align}

+++

Since the left side depends only on $t$ and the right side only on $x$, both sides must equal a constant $-\lambda$:
\begin{align}
\frac{1}{T(t)} \frac{dT(t)}{dt} = -\lambda = \alpha \frac{1}{X(x)} \frac{d^2 X(x)}{dx^2}
\end{align}
This yields two ordinary differential equations (ODEs):

1. Temporal Equation:
\begin{align}
\frac{dT(t)}{dt} + \lambda T(t) = 0
\end{align}

2. Spatial Equation:
\begin{align}
\frac{d^2 X(x)}{dx^2} + \frac{\lambda}{\alpha} X(x) = 0
\end{align}

+++

The general solution to the spatial equation is:
\begin{align}
X(x) = A \sin(kx) + B \cos(kx)
\end{align}
where $k^2 = \lambda/\alpha$.
Note that this form of solution is originated from the second-order derivative.

Assuming Dirichlet boundary conditions for a rod of length $L$
\begin{align}
u(0,t) = u(L,t) = 0,
\end{align}
At $x = 0$, $X(0) = 0 \implies B = 0$.
At $x = L$, all non-trivial solutions ($A \neq 0$) require:
\begin{align}
\sin(kL) = 0 \implies kL = n\pi, \quad n = 1,2,3,\dots
\end{align}
Thus, the "eigenvalues" are
\begin{align}
k_n = \frac{n\pi}{L}
\end{align}
and the eigenfunctions are 
\begin{align}
X_n(x) = \sin\left( \frac{n\pi x}{L} \right).
\end{align}

+++

With $\lambda_n = \alpha k_n^2$, the temporal ODE becomes
\begin{align}
\frac{dT_n(t)}{dt} + \alpha \left( \frac{n\pi}{L} \right)^2 T_n(t) = 0.
\end{align}
with solution
\begin{align}
T_n(t) = C_n \exp\left[ -\alpha \left( \frac{n\pi}{L} \right)^2 t \right].
\end{align}

+++

Combining spatial and temporal parts and realizing the heat equation is linear, the general solution is the sum of all solutions:
\begin{align}
u(x,t) = \sum_{n=1}^\infty C_n \sin\left( \frac{n\pi x}{L} \right) \exp\left[ -\alpha \left( \frac{n\pi}{L} \right)^2 t \right]
\end{align}
and the coefficients $C_n$ are determined from the initial condition $u(x,0) = f(x)$:
\begin{align}
C_n = \frac{2}{L} \int_0^L f(x) \sin\left( \frac{n\pi x}{L} \right) dx.
\end{align}
This represents the Fourier sine series expansion of $f(x)$.

+++

## Fourier Series: Foundation and Interpretation

## Implementing Fourier Series in Python

## Transition to Fourier Transform

## Sampling Theory and DFT

## FFT and Computational Efficiency

## Fourier Transform and the Heat Equation

## Spectral Derivatives

## Convolution and Parseval's Theorem

## Denoising and Signal Processing Applications

## Time-Frequency Analysis

## Astrophysical Applications and VLBI

## Conclusion and Further Resources

```{code-cell} ipython3

```

```{code-cell} ipython3

```
