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

### Definition and Motivation

In many physical systems, including those encountered in astrophysics, functions describing phenomena are often periodic or can be approximated as such over certain intervals.
The Fourier Series provides a powerful method to represent these periodic functions as infinite sums of simpler sinusoidal functions---sines and cosines.
This enables us to break down complex periodic functions and simplifies analysis and computation.
Since nany physical systems can be described by harmonic components, making Fourier Series naturally enables some physical interpretation.
In astrophysics, signals received from observations can be decomposed to analyze their frequency content.

A function $f(x)$ is said to be periodic with period $L$ if:
\begin{align}
    f(x + L) = f(x) \quad \text{for all } x.
\end{align}
Common examples include trigonometric functions like $\sin(x)$ and $\cos(x)$, which have a period of $2\pi$.

Any reasonably well-behaved, periodic function $f(x)$ with period $L$ can be expressed as a Fourier Series:
\begin{align}
f(x) = \frac{A_0}{2} + \sum_{n=1}^\infty \left[ A_n \cos\left( \frac{2n\pi x}{L} \right) + B_n \sin\left( \frac{2n\pi x}{L} \right) \right].
\end{align}
$A_n$ and $B_n$ are the Fourier coefficients, which quantify the contribution of each harmonic component.

The coefficients are calculated using the orthogonality properties of sine and cosine functions:
\begin{align}
A_n &= \frac{2}{L} \int_{-L/2}^{L/2} f(x) \cos\left( \frac{2n\pi x}{L} \right) dx, \quad n = 0,1,2,\dots \\
B_n &= \frac{2}{L} \int_{-L/2}^{L/2} f(x) \sin\left( \frac{2n\pi x}{L} \right) dx, \quad n = 1,2,3,\dots
\end{align}

+++

### Example

Consider a square wave function defined over the interval $[-L/2, L/2)$:
\begin{align}
  f(x) =
  \begin{cases}
     1, &  0   < x < L/2, \\
    -1, & -L/2 < x < 0.
  \end{cases}
\end{align}
The Fourier coefficients for this function can be computed using the integrals above:
\begin{align}
  B_n = \frac{4}{n\pi}
  \begin{cases}
    0, & \text{even } n, \\
    1, & \text{odd }  n.
  \end{cases}
\end{align}

```{code-cell} ipython3
:jp-MarkdownHeadingCollapsed: true

import numpy as np
import matplotlib.pyplot as plt

def f(x, L):
    return np.where((x % L) < L/2, 1, -1)

def f_approx(x, L, N):
    fsum = np.zeros_like(x)
    for n in range(1, N + 1, 2):  # Sum over odd n
        B = 4 / (n * np.pi)
        fsum += B * np.sin(2 * n * np.pi * x / L)
    return fsum
```

```{code-cell} ipython3
L  = 2 * np.pi  # Period is 2L
xi = np.linspace(-L/2, L/2, 1000)

# Plotting
plt.figure(figsize=(10, 6))

# Original function
fi = f(xi, L)
plt.plot(xi, fi, label='Square Wave', color='k')

# Fourier series approximation
N = list(range(5,100,5))
for n in N[:5]:
    fi_n = f_approx(xi, L, n)
    plt.plot(xi, fi_n, label=f'Fourier Series Approximation (N={n})')

plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)
```

Note that, near points of discontinuity in the function (e.g., the jumps in a square wave), the Fourier series overshoots the function's value.
This overshoot does not diminish as more terms are added;
instead, the maximum overshoot approaches a finite limit ($\sim 9\%$ of the jump's magnitude).
This so called **Gibbs Phenomenon** states that, while the Fourier series converges pointwise almost everywhere, it does so non-uniformly at discontinuities.

+++

## Implementing Fourier Series in Python

We will implement the Fourier Series representation of functions using Python.
We will compute the Fourier coefficients numerically and reconstruct the original function from these coefficients.

```{code-cell} ipython3
def An(f, x, L, n):
    I = f(x) * np.cos(2 * n * np.pi * x / L)
    return (2 / L) * np.trapezoid(I, x)

def Bn(f, x, L, n):
    I = f(x) * np.sin(2 * n * np.pi * x / L)
    return (2 / L) * np.trapezoid(I, x)
```

```{code-cell} ipython3
def Fourier_coefficients(f, x, L, N):
    A = [An(f, x, L, n) for n in range(0, N)]
    B = [Bn(f, x, L, n) for n in range(0, N)]
    return A, B
```

```{code-cell} ipython3
def Fourier_series(A, B, L, x):
    fsum = (A[0]/2) * np.ones_like(x)
    for n, An in enumerate(A[1:],1):
        fsum += An * np.cos(2 * n * np.pi * x / L)
    for n, Bn in enumerate(B[1:],1):
        fsum += Bn * np.sin(2 * n * np.pi * x / L)
    return fsum    
```

Now, we can obtain Fourier series numerically using arbitrary functions.

```{code-cell} ipython3
L = 2 * np.pi

def f(x): # closure on L
    a = -L/2
    b =  L/2
    x = (x - a) % (b - a) + a
    return np.exp(-x*x*32)

xi = np.linspace(-L/2, L/2, 10_000)
A, B = Fourier_coefficients(f, xi, L, 100)
```

```{code-cell} ipython3
plt.figure(figsize=(10, 6))

# Original function
xi = np.linspace(-L, L, 20_000)
fi = f(xi)
plt.plot(xi, fi, color='k', label='Original function')

# Fourier series approximation
N    = list(range(5,100,5))
fi_N = [Fourier_series(A[:n], B[:n], L, xi) for n in N]
for n, fi_n in list(zip(N, fi_N))[:5]:
    plt.plot(xi, fi_n, label=f'Fourier Series Approximation (N={n})')

plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)
```

### Error Analysis

We can quantify how the approximation improves with $N$ by calculating the Mean Squared Error (MSE) between the original function and its Fourier series approximation.

```{code-cell} ipython3
def mse(f, f_N):
    df = (f - f_N)
    return np.sqrt(np.mean(df * df))

errs = []
for fi_n in fi_N:
    errs.append(mse(fi, fi_n))

plt.loglog(N, errs, label='Error')
plt.loglog(N, 1/np.array(N), label=r'$N^{-1}$')
plt.legend()
```

Try adjust different functions and observe how the errors behave.

+++

## Complex Fourier Series

It is often convenient to combine the sine and cosine in the Fourier series.
Recalling Euler's formula is:
\begin{align}
e^{i\theta} = \cos(\theta) + i\sin(\theta).
\end{align}
Therefore,
\begin{align}
\cos(\theta) &= \frac{e^{i\theta} + e^{-i\theta}}{2}, \\
\sin(\theta) &= \frac{e^{i\theta} - e^{-i\theta}}{2i}.
\end{align}
Substituting these into the definition of Fourier series, we obtain the Complex Fourier Series:
\begin{align}
f(x) = \sum_{n=-\infty}^{\infty} C_n e^{i n \omega_1 x},
\end{align}
where $\omega_1 = 2\pi/L$ is the fundamental angular frequency.

The complex coefficients $C_n$ are given by:
\begin{align}
C_n = \frac{1}{L} \int_{-L/2}^{L/2} f(x) e^{-i n \omega_1 x} dx.
\end{align}

What are the relationship between the complex Fourier coefficients $C_n$ and the Fourier series coefficients $A_n$ and $B_n$?
What is a special property of $C_n$ is the function $f(x)$ if purely real or purely imaginary?

+++

## Transition to Fourier Transform

In the previous sections, we explored how periodic functions can be represented as sums of sines and cosines using Fourier Series.
However, many functions of interest in physics and engineering are not periodic or are defined over an infinite domain.
To analyze such functions, we need to extend the concept of Fourier Series to the Fourier Transform.

### From Discrete to Continuous Spectrum

Recalling for a function $f(x)$ with period $L$, the (Complex) Fourier Series has discrete frequencies $\omega_n = 2n\pi/L \equiv n \omega_1$.

When the period $L$ becomes infinitely large ($L \rightarrow \infty$), the spacing between the frequencies in the Fourier Series becomes infinitesimally small, turning the sum into an integral.

\begin{align}
\sum_{n=\infty}^{\infty} \rightarrow
\int_{-\infty}^{\infty} \frac{d\omega}{\omega_1}
\end{align}

The discrete set of frequencies $\omega_n$ becomes a continuous variable $\omega$.
This limit leads us to the Fourier Transform, which represents non-periodic functions defined over the entire real line.

### Complex Fourier Transform Definitions

The Fourier Transform $F(\omega)$ of a function $f(x)$ is defined as:
\begin{align}
F(\omega) = \int_{-\infty}^{\infty} f(x) e^{-i \omega x} dx.
\end{align}
The Inverse Fourier Transform is:
\begin{align}
f(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} F(\omega) e^{i \omega x} d\omega.
\end{align}
These equations allow us to transform a time-domain function $f(x)$ into its frequency-domain representation $F(\omega)$, and vice versa.

What is a special property of $F(\omega)$ if the function $f(x)$ is purely real or purely imaginary?

Because computers have finite memory, we almost never use the true Fourier transform in computational science.
Nevertheless, it is extremely useful in proofs, study numerical properties of different algorithms, and derive corrections to algorithms.

+++

## Sampling Theory and DFT

This section focuses on the implications of sampling, introduces the Nyquist Sampling Theorem, explores the phenomenon of aliasing, and explains the role of the Discrete Fourier Transform (DFT) in analyzing discrete signals.
We will switch our notation from $f(x)$ to $f(t)$.

### Sampling Continuous Signals

Continuous vs. Discrete Signals:
* Continuous Signal: A function $f(t)$ defined for all real values of $t$.
* Discrete Signal: A sequence of values $f_n = f(n T_\text{s})$, where $T_\text{s}$ is the sampling interval/time and $n$ is an integer index.
We call $f_s = 1/\Delta t$ and $\omega_\text{s} = 2\pi f_\text{s}$ the sampling frequency and sampling angular frequency, respectively.

The Nyquist Sampling Theorem states:

A band-limited continuous-time signal with maximum frequency $f_\max$ can be perfectly reconstructed from its samples if the sampling frequency $f_\text{s}$ is greater than twice the maximum frequency of the signal:
\begin{align}
f_\text{s} > 2 f_\max.
\end{align}
The minimum sampling rate $f_{\text{Nyquist}} = 2 f_\max$ is called the Nyquist Rate.
And half the sampling frequency $f_{\text{Nyquist}} = f_\text{s}/2 is called the 
Nyquist Frequency.

Sampling above the Nyquist rate ensures that high-frequency components do not overlap and cause distortion, preventing **aliasing**.
Under ideal conditions, a band-limited signal can be reconstructed exactly from its samples, resulting perfect reconstruction.

+++

### Aliasing Errors

Aliasing occurs when a continuous signal is sampled at a rate insufficient to capture its changes accurately, causing different frequency components to become indistinguishable in the sampled signal.
Frequencies higher than the Nyquist frequency are "folded" back into lowe frequencies, making them appeare as low frequencies.

Let's illustrate aliasing with an example:

```{code-cell} ipython3
# Parameters
f_signal = 3  # Frequency of the signal in Hz
f_conf   = 1  # Frequency of the confused signal in Hz
f_s      = 4  # Low sampling frequency in Hz (below Nyquist rate)
t_end    = 2   # Duration in seconds

# Time arrays
t_cont = np.linspace(0, t_end, 10_000)
t_s    = np.arange(0, t_end, 1 / f_s)

# Signal
x_cont =  np.sin(2 * np.pi * f_signal * t_cont)
x_conf = -np.sin(2 * np.pi * f_conf   * t_cont)
x_s    =  np.sin(2 * np.pi * f_signal * t_s)

# Plotting
plt.figure(figsize=(12, 6))
plt.plot(t_cont, x_cont, label='Original Signal', color='k')
plt.plot(t_cont, x_conf, label='Confused Signal', color='b')
plt.scatter(t_s, x_s, label=f'Sampled at {f_s} Hz', color='red')

plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)
```

### Discrete Fourier Transform (DFT)

Fourier transforming discrete samples $f_n$ of function $f(t)$ at $t = n T_\text{s}$ results the Discrete Fourier Transform (DFT).
It transforms a sequence of $N$ (complex) numbers $f_n$ into another sequence of  $N$ complex numbers $F_k$:
\begin{align}
F_k = \sum_{n=0}^{N-1} f_n e^{-i 2\pi k n / N}, \quad k = 0, 1, …, N-1.
\end{align}
The inverse DFT is given by:
\begin{align}
f_n = \frac{1}{N} \sum_{k=0}^{N-1} F_k e^{i 2\pi k n / N}, \quad n = 0, 1, \dots, N-1.
\end{align}

The Discrete Fourier Transform possesses several key properties that make it an essential tool in signal processing:

DFT is periodic.
Both the input sequence $f_n$ and the output sequence $F_k$ are periodic with period $N$.
This means that $f_{n+N} = f_n$ and $F_{k+N} = F_k$ for all integers $n$ and $k$.
The periodicity arises because the DFT assumes that the finite-length sequence represents one period of an infinitely periodic signal.
As a result, the DFT effectively analyzes the frequency content within this repeating structure.

Another important property is the orthogonality of the basis functions used in the DFT.
The complex exponentials $\exp(-i 2\pi k n / N)$ serve as the basis functions, and they are orthogonal over the discrete set of points $n = 0, 1, \dots, N-1$.
Orthogonality implies that the inner product (sum over $n$) of different basis functions is zero:
\begin{align}
\sum_{n=0}^{N-1} e^{-i 2\pi k n / N} e^{i 2\pi k' n / N} = 0, \quad \text{for } k \ne k'.
\end{align}
This property ensures that each frequency component represented by $F_k$ is independent of the others, allowing for the separation and analysis of individual frequency components within the signal.

Furthermore, the DFT provides an exact transformation of the discrete signal into the frequency domain.
Unlike numerical approximations of the continuous Fourier Transform, the DFT computes the precise frequency components of a finite-length discrete signal, given the assumption of periodicity.
This means that no information is lost or approximated during the transformation, and the original signal can be perfectly reconstructed from its DFT coefficients using the inverse DFT.

+++

### Implementing the DFT in Python

Although there are many Fourier transform packages (notiable `numpy.fft`), implementing the DFT directly helps to understand its underlying mechanics.

A direct implementation of the DFT in Python involves computing the summations in the DFT definition for each value of $k$:

```{code-cell} ipython3
def DFT(f):
    N = len(f)
    F = np.zeros(N, dtype=complex)
    for k in range(N):
        for n in range(N):
            F[k] += f[n] * np.exp(-2j * np.pi * k * n / N)
    return F
```

This implementation has a computational complexity of $\mathcal{O}(N^2)$, making it inefficient for large $N$.
However, it is valuable for educational purposes to see how each frequency component $F_k$ is calculated from the input sequence $f_n$.

Consider a simple sine wave signal, where the DFT accurately identifies the frequency component at 5 Hz, as expected.
The magnitude spectrum shows a peak at this frequency, demonstrating that the DFT can extract the frequency information from the sampled signal.

```{code-cell} ipython3
f_signal = 5 # Frequency in Hz
f_s = 50     # Sampling frequency in Hz
N = 100      # Number of samples

t = np.arange(N) / f_s
f = np.sin(2 * np.pi * f_signal * t)

# Compute DFT
F = DFT(f)
freq = np.fft.fftfreq(N, d=1/f_s)

# Plotting
plt.figure(figsize=(10, 6))
plt.stem(freq[:N // 2], np.abs(F[:N // 2])*(2/N), basefmt=" ")
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.grid(True)
```

### The DFT Matrix Representation

The DFT can also be expressed in matrix form, which provides a different perspective on the transformation. The DFT of a sequence $f_n$ can be represented as a matrix multiplication:
\begin{align}
\mathbf{F} = \mathbf{W} \mathbf{f},
\end{align}
where
$\mathbf{f}$ is the column vector of input samples $f_n$,
$\mathbf{F}$ is the column vector of DFT coefficients $F_k$, and
$\mathbf{W}$ is the DFT matrix with elements $W_{k n} = \exp(-2\pi i k n / N)$.
This formulation highlights that the DFT is a linear transformation from the time domain to the frequency domain.

Implementing the DFT using matrix multiplication in Python:

```{code-cell} ipython3
def DFT_matrix(N):
    n = np.arange(N)
    k = n.reshape((N, 1))
    W = np.exp(-2j * np.pi * k * n / N)
    return W

# Compute DFT using matrix multiplication
W = DFT_matrix(N)
```

```{code-cell} ipython3
F_matrix = np.dot(W, f)

# Plotting
plt.figure(figsize=(10, 6))
plt.stem(freq[:N // 2], np.abs(F_matrix[:N // 2])*(2/N), basefmt=" ")
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.grid(True)
```

This method produces the same result as the direct computation but is also computationally intensive with  $\mathcal{O}(N^2)$ complexity.
The matrix representation is useful for theoretical analysis and understanding properties of the DFT, such as its eigenvalues and eigenvectors.

+++

### Aliasing in the DFT

The periodicity inherent in the DFT can lead to aliasing when analyzing signals that contain frequency components higher than the Nyquist frequency.
Since the DFT assumes the input signal is one period of a periodic signal, any frequencies above $f_\text{s}/2$ will be indistinguishable from lower frequencies due to the wrapping around in the frequency domain.

In practical applications, the maximum frequency that can be accurately represented without aliasing is  $f_\max = f_\text{s}/2$.
If the signal contains components beyond this limit, aliasing will occur, and the DFT will incorrectly attribute these high frequencies to lower ones.
Therefore, it is essential to ensure that the sampling rate is sufficiently high or that the signal is filtered to remove high-frequency components before sampling.

+++

### Preventing Aliasing

To prevent aliasing, an anti-aliasing filter is applied to the continuous signal before sampling.
This filter is typically a low-pass filter that attenuates frequency components above the Nyquist frequency, ensuring that the sampled signal does not contain frequencies that could cause aliasing.
In hardware, analog filters are used for this purpose, conditioning the signal prior to digitization.

Another approach to preventing aliasing is to sample at a rate significantly higher than the Nyquist rate, a practice known as oversampling.
Oversampling provides a buffer against aliasing by increasing the Nyquist frequency.
Additionally, adaptive sampling techniques can adjust the sampling rate based on the signal's frequency content, sampling more rapidly when high-frequency components are present.

+++

## FFT and Computational Efficiency

## Fourier Transform and the Heat Equation

## Spectral Derivatives

## Convolution and Parseval's Theorem

## Denoising and Signal Processing Applications

## Time-Frequency Analysis

## Astrophysical Applications and VLBI

## Conclusion and Further Resources
