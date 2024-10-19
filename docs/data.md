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

# Data Representation and Errors

```{epigraph}
The decimal system has been established, somewhat foolishly to be
sure, according to man's custom, not from a natural necessity as most
people think.

-- Blaise Pascal (1623--1662)
```

+++

## A Moment of ZEN

Q: What is the output of this simple code?
What is the last value of `x` it will print?
```c
/* To run me, save into "demo.c" and then compile with `gcc demo.c -o demo` */

#include <stdio.h>

int
main(int argc, char *argv[])
{
	float x;

	for(x = 0.0; x <= 1.0; x += 0.1)
		printf("x = %f, f(x) = %f\n", x, x * x);

	return 0;
}
```

+++

To see the problem, change the `printf` statement to
```c
printf("x = %10.8f, f(x) = %10.8f\n", x, x * x);
```

+++

Then, change the increment to `0.125`.

+++

# Rounding Errors Can Have Fatal Consequences

On February 25, 1991, during Operation Desert Storm, a Scud missile fired by the Iraqi army struck a U.S. army barracks at Dhahran Air Base in Saudi Arabia, killing 28 soldiers.
Although the base was protected by a Patriot Air Defense System, it failed to track and intercept the incoming missile.
A subsequent investigation by the U.S. Army revealed that the failure was caused by an accumulation of rounding errors in the Patriot's tracking software.

Each Patriot battery consists of a ground-based radar and eight missile launchers.
The radar detects airborne objects and tracks their movement, while the system's weapons control computer calculates the expected trajectory of any detected missile.
If an object follows the predicted ballistic trajectory, the system launches a Patriot missile to intercept it.
To perform these calculations, the computer relies on an internal clock, which tracks time in increments of 0.1 seconds, stored in 24-bit registers.

However, in binary form, 0.1 seconds is represented as an infinite repeating decimal, $(0.1)_{10} = (0.0001100011…)_2$, and must be rounded during computations.
After approximately 300 hours of continuous operation, the accumulated rounding errors became significant enough that the system failed to accurately predict the Scud missile's trajectory.
As a result, the Patriot system did not identify the missile as a threat, and no interception occurred.

+++

## Another Moment of ZEN

We all learned in high school that the solutions (roots) to the qudratic equation $a x^2 + b x + c = 0$ is
$$
x = \frac{-b \pm \sqrt{b^2 - 4 a c}}{2a}
$$

Q: Why one of the roots become zero when solving the qudratic equation with $b = 1$ and $a = c = 10^{-9}$?

```{code-cell} ipython3
a = 1e-9
b = 1
c = 1e-9

x1 = (-b + (b*b - 4*a*c)**(1/2)) / (2*a)
x2 = (-b - (b*b - 4*a*c)**(1/2)) / (2*a)

print(f'{x1:.16f}, {x2:.16f}')
```

It is straightforward to show in the limit $a, c \ll b$, the roots are
$$
x \approx -\frac{b}{a} \mbox{ or } -\frac{c}{b}
$$
Is it possible to recover the small root $-c/b$?

+++

When $b > 0$, a catastropic cancellation (see below) happens only in the "+" equation.
We may replace the first qudratic equation by its "conjugate" form
$$
x = \frac{2c}{-b \mp \sqrt{b^2 - 4 a c}}
$$

```{code-cell} ipython3
x1 = (2*c) / (-b - (b*b - 4*a*c)**(1/2))
x2 = (-b - (b*b - 4*a*c)**(1/2)) / (2*a)

print(f'{x1:.16f}, {x2:.16f}')
```

Equivalently, we may use the "numerically stable form",
\begin{align}
x_1 &= \frac{-b - \mathrm{sign}(b)\sqrt{b^2 - 4 a c}}{2a} \\
x_2 &= \frac{c}{a x_1}
\end{align}
as used by
[GSL](https://git.savannah.gnu.org/cgit/gsl.git/tree/poly/solve_quadratic.c#n57) and 
[fadge](https://github.com/adxsrc/fadge/blob/main/mod/fadge/utils.py#L25).

+++

## Catastropic Cancellation
