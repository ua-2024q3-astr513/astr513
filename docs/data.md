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

## Rounding Errors Can Have Fatal Consequences

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

## The Challenge of Representing Numbers

As we saw with the missile story, even small errors in representing numbers can accumulate with devastating consequences.
The issue of how to efficiently and accurately represent numbers is not trivial.
In fact, the introduction of the Arabic numeral system was a major advancement over earlier methods like Roman numerals---and even more primitive systems like unary numbers.

The **unary number system** is the most basic of all, where each number is represented by a corresponding number of identical marks.
For example, the number 5 is written as "|||||".
While simple, unary numbers are extremely inefficient for large quantities.
Representing the number 888, for instance, would require writing 888 marks!
Unary numbers are rarely used in practice because they require a lot of space to represent even moderately sized numbers.

Next, consider **Roman numerals***, which advanced by using specific symbols (I, V, X, L, C, D, M) to represent numbers like 1, 5, 10, 50, 100, 500, and 1000.
To represent a number like 888 in Roman numerals, you'd need 12 symbols: DCCCLXXXVIII.
While more efficient than unary, it's still cumbersome for large numbers.

In contrast, the **Arabic number system** allows any number from 1 to 999 to be represented using at most 3 digits.
While the notation is trivial for us, this was a major breakthrough in number representation, as it reduced the complexity of writing numbers while greatly increasing the efficiency and clarity.

+++

## Positional Notation Systems

The Arabic numeral system is an example of a positional notation system, in which the value of a digit is determined by both the digit itself and its position within the number.
This is different from systems like Roman numerals or unary numbers, where the position of a symbol doesn’t change its value.
In positional notation, the place of a digit corresponds to a specific power of the base of the system.

In any positional system, to represent a number:
* we first decide on the base (or radix) $b$.
* then identify the notation for the digits
* inally write $\pm (\dots d_3 d_2 d_1 d_0 . d_{-1} d_{-2} d_{-3} \dots)$ to represent $\pm (\dots + d_3 b^3 + d_2 b^2 + d_1 b^1 + d_0 b^0 + d_{-1} b^{-1} + d_{-2} b^{-2} + d_{-3} b^{-3} + \dots)$.

To convert from a numeral system of base $b$ to the decimal one, we simply use the definition 

$$
\pm (\dots + d_3 b^3 + d_2 b^2 + d_1 b^1 + d_0 b^0 + d_{-1} b^{-1} + d_{-2} b^{-2} + d_{-3} b^{-3} + \dots).
$$

Example:

$$
(256.4)_8 = 2\times8^2 + 5\times8^1 + 6 + 4\times8^{-1} = (174.5)_{10}
$$

+++

## Binary Numbers

* Base: $b = 2$
* Digits: $0$, $1$

The binary system has been used in various forms long before the age of computers.
Invented by merchants in medieval England, the units of liquid measure were based on the binary system. For example:
* 1 gallon = 2 pottles;
* 1 pottle = 2 quarts;
* 1 quart = 2 pints;
* 1 pint = 2 cups; etc.

![measure](figures/measure.png)

Similarly, the binary system is used in music to define note durations, i.e., whole note, half note, quarter note, eighth note, sixteenth note, etc.
These everyday examples reflect the fundamental nature of the binary system, which underpins much of modern computing.

In the binary system, only two digits are used: 0 and 1.
The position of each digit in a binary number corresponds to a power of 2, just as the position of a digit in the decimal system corresponds to a power of 10. 
For example, the binary number $1011_2$ represents: $1 \times 2^3 + 0 \times 2^2 + 1 \times 2^1 + 1 \times 2^0$.
This gives the decimal value: $1 \times 8 + 0 \times 4 + 1 \times 2 + 1 \times 1 = 11$.

![binary-shirt](figures/binary-shirt.jpg)


+++

## Another Moment of ZEN

We all learned in high school that the solutions (roots) to the qudratic equation $a x^2 + b x + c = 0$ is

$$
x = \frac{-b \pm \sqrt{b^2 - 4 a c}}{2a}.
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

$$
x_1 = \frac{-b - \mathrm{sign}(b)\sqrt{b^2 - 4 a c}}{2a} 
$$

$$
x_2 = \frac{c}{a x_1}
$$

as used by
[GSL](https://git.savannah.gnu.org/cgit/gsl.git/tree/poly/solve_quadratic.c#n57) and 
[fadge](https://github.com/adxsrc/fadge/blob/main/mod/fadge/utils.py#L25).

+++

## Catastropic Cancellation
