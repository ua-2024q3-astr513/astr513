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

```{note}
### Rounding Errors Can Have Fatal Consequences

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
```

+++

## Representing Numbers

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

### Positional Notation Systems

The Arabic numeral system is an example of a positional notation system, in which the value of a digit is determined by both the digit itself and its position within the number.
This is different from systems like Roman numerals or unary numbers, where the position of a symbol doesn’t change its value.
In positional notation, the place of a digit corresponds to a specific power of the base of the system.

In any positional system, to represent a number:
* we first decide on the base (or radix) $b$.
* then identify the notation for the digits
* inally write $\pm (\dots d_3 d_2 d_1 d_0 . d_{-1} d_{-2} d_{-3} \dots)$ to represent $\pm (\dots + d_3 b^3 + d_2 b^2 + d_1 b^1 + d_0 b^0 + d_{-1} b^{-1} + d_{-2} b^{-2} + d_{-3} b^{-3} + \dots)$.

To convert from a numeral system of base $b$ to the decimal one, we simply use the definition 
\begin{align}
\pm (\dots + d_3 b^3 + d_2 b^2 + d_1 b^1 + d_0 b^0 + d_{-1} b^{-1} + d_{-2} b^{-2} + d_{-3} b^{-3} + \dots).
\end{align}
Example:
\begin{align}
(256.4)_8 = 2\times8^2 + 5\times8^1 + 6 + 4\times8^{-1} = (174.5)_{10}
\end{align}

+++

### Binary Numbers

* Base: $b = 2$
* Digits: 0, 1

```{figure} figures/measure.png
---
height: 320px
align: right
---
Binary system was invented by merchants in medieval England.
```

The binary system has been used in various forms long before the age of computers.
Invented by merchants in medieval England, the units of liquid measure were based on the binary system. For example:
* 1 gallon = 2 pottles;
* 1 pottle = 2 quarts;
* 1 quart = 2 pints;
* 1 pint = 2 cups; etc.

Similarly, the binary system is used in music to define note durations, i.e., whole note, half note, quarter note, eighth note, sixteenth note, etc.
These everyday examples reflect the fundamental nature of the binary system, which underpins much of modern computing.

```{figure} figures/binary-shirt.jpg
---
height: 320px
---
There are 10 types of people in the world...
```

In the binary system, only two digits are used: 0 and 1.
The position of each digit in a binary number corresponds to a power of 2, just as the position of a digit in the decimal system corresponds to a power of 10. 
For example, the binary number $1011_2$ represents: $1 \times 2^3 + 0 \times 2^2 + 1 \times 2^1 + 1 \times 2^0$.
This gives the decimal value: $1 \times 8 + 0 \times 4 + 1 \times 2 + 1 \times 1 = 11$.

+++

### The Hexadecimal System

* Base: $b = 16$
* Digits: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F

The hexadecimal system allows for writing a binary number in a very compact notation.

```{figure} figures/color.png
---
height: 320px
---
Hex numbers are used to select colors.
```

+++

```{note}

### Quantum Computing: A New Level of Number Representation

Quantum computing takes number representation to a whole new level.
In quantum mechanics, data is represented using quantum bits, or qubits.
Unlike classical bits, a single qubit can exist in a superposition of two states, $|0\rangle$ and $|1\rangle$, simultaneously.
This means that (ignoring normalization for now) one qubit can represent two numbers $C_0$ and $C_1$ as in $C_0 |0\rangle + C_1 |1\rangle$.

With two qubits, the system can represent four possible states: $|00\rangle$, $|01\rangle$, $|10\rangle$, and $|11\rangle$, again, in superposition.
The number of possible states grows exponentially as more qubits are added:
* Three qubits can represent eight states: $|000\rangle$, $|001\rangle$, $|010\rangle$, ..., $|111\rangle$.
* Four qubits represent 16 states, and so on.

In general, $n$ qubits can represent $2^n$ states at once.
This exponential scaling is what gives quantum computers their enormous potential to perform certain types of calculations far more efficiently than classical computers.

For example, IBM's 53-qubit quantum computer can represent $2^{53}$ states simultaneously.
In terms of classical information, this is equivalent to storing approximately 1 petabyte (PB) of data, which is comparable to the amount of memory available in some of the world's largest supercomputers today.
```

+++

```{note}

### Superposition Hypothesis in Large Language Models

In large language models (LLMs), the
[Superposition Hypothesis](https://transformer-circuits.pub/2022/toy_model/index.html) suggests that individual neurons or parameters can store multiple overlapping features.
This is similar to how, in quantum computing, qubits exist in superposition.
Instead of each neuron being dedicated to a single concept, neurons in LLMs can represent different features depending on the context.
For example, a neuron might activate for both "dog" and "pet", depending on the input, allowing the model to efficiently reuse its parameters for different purposes.

This efficiency is further enhanced by the use of high-dimensional embeddings, where words and phrases are mapped to vectors in a space with many dimensions.
The [Johnson–Lindenstrauss Lemma](https://en.wikipedia.org/wiki/Johnson%E2%80%93Lindenstrauss_lemma) shows that in such high-dimensional spaces, exponentially many almost orthogonal directions can exist.
This means that as the dimensionality of the space increases, the model can represent and distinguish between a vast number of concepts without interference, as these nearly orthogonal vectors can encode distinct meanings.

The Superposition Hypothesis may help explain the remarkable power of LLMs, particularly in how they manage to store and process such vast amounts of information efficiently.
However, while this hypothesis offers promising insights, it is not yet conclusive, and ongoing research is required to fully understand the underlying mechanisms that drive the success of these models.
```

+++

### Addition in the Binary System

In the binary system, addition follows similar rules to decimal addition, but with only two digits: 0 and 1.
The key rules are:
```{code}
0 + 0 = 0
0 + 1 = 1
1 + 0 = 1
1 + 1 = 0 + 1 carryover
```

We start from the rightmost bit, adding the bits with carry when needed:
```{code}
    C = 111
    A = 0101 (  A = 5)
    B = 0011 (  B = 3)
-------------
A + B = 1000 (A+B = 8)
```

+++

### Multiplication in the Bindary System

A binary multiplication is simply shift and add:
```{code}
       11001010
     x 01001001
     ----------
       11001010
      00000000
     00000000
    11001010
   00000000
  00000000
 11001010
00000000
---------------
011100110011010
```

+++

## Hardware Implementations

### CMOS: Complementary Metal-Oxide Semiconductor

```{figure} figures/cmos.png
---
height: 160px
---
N-type and P-type MOS have different properties.
```

| PMOS | NMOS |
| --- | --- |
| $g=0 \Rightarrow s=0$ | $g=0 \Rightarrow s=1$ |
| $g=1 \Rightarrow s=1$ | $g=1 \Rightarrow s=0$ |

We may combine PMOS and NMOS gates to create more complex logic gates, e.g., an XOR gate.

```{figure} figures/xor.png
---
height: 160px
---
PMOS and NMOS can be combined to create, e.g., XOR gate.
```

```{code}
0 XOR 0 = 0
0 XOR 1 = 1
1 XOR 0 = 1
1 XOR 1 = 0
```

We can combine complex logic gates to make "half" adders.
"Sum" is implemented by XOR, "carry" is implemented by AND.

```{figure} figures/half.png
---
height: 160px
---
"Half" adder.
```

```{code}
A + B = S , C
-------------
0 + 0 = 0 , 0
0 + 1 = 1 , 0
1 + 0 = 1 , 0
1 + 1 = 0 , 1
```

We can combine two "half" adders into a "full" adder.

```{figure} figures/full.png
---
height: 160px
---
"Full" adder.
```

Then, we can combine "full" adders to make multi-digit adders.

```{figure} figures/multi.png
---
height: 160px
---
Multi-digit adder.
```

+++

### Moore's Law

```{figure} figures/Moores_Law.png
---
height: 640px
---
Moore's Law.
```

+++

### UA HPC

```{figure} figures/elgato.jpg
---
height: 240px
---
Elgato
```

```{figure} figures/ocelote.jpg
---
height: 240px
---
Ocelote
```

+++

### What Determines Computational Performance?

```{figure} figures/42-years-processor-trend.png
---
height: 640px
---
Origin data up to year 2010 collected and plotted by Horowitz et al.; new data collected and plotted by [Rupp](https://www.karlrupp.net/2018/02/42-years-of-microprocessor-trend-data/).
```

+++

## How Do Computers Handle Non-Integer Numbers?

### Floating Point Representation

The easiest way to describe floating-point representation is through an example.
Consider the result of the mathematical expression $e^6 \approx 403.42879$.
To express this in normalized floating-point notation, we first write the number in scientific notation:
\begin{align}
e^6 = 4.0342879 \times 10^2
\end{align}
In scientific notation, the number is written such that the significand (or mantissa) is always smaller than the base (in this case, 10).
To represent the number in floating-point format, we store the following components:
* The sign of the number,
* The exponent (the power of 10),
* The significand (the string of significant digits).

For example, the floating-point representation of $e^6$ with 4 significant digits is:
\begin{align}
e^6 = (+, 2, 4034).
\end{align}
And with 8 significant digits:
\begin{align}
e^6 = (+, 2, 40342879).
\end{align}

+++

### Single-Precision Floating Point

```{figure} figures/float.svg
---
height: 80px
---
The value $+0.15625 = (-1)^{(0)_2} \times 2^{(01111100)_2-127} \times (1.01...0)_2$ stored as single precision float.
```

In the IEEE 754 standard for floating-point arithmetic, used by most modern computers, special rules are applied to store numbers.
In single precision (32-bit), the significand is stored with an implicit leading bit to the left of the binary point, which is always assumed to be 1 in normalized numbers.
This means that, instead of storing the full significand, only the fractional part (digits to the right of the binary point) is stored.
For example, the number `1.101` in binary would be stored as just `101` in the significand, with the leading 1 implied.
This optimization increases precision, as it effectively adds an extra bit without needing more storage.

In addition, the exponent is stored with a "bias" of 127 in single precision.
The actual exponent is offset by this bias, allowing both positive and negative exponents to be represented.
The smallest exponent that can be stored is 0 (representing an actual exponent of -127), and the largest is 255 (representing an actual exponent of +128).
Thus, single-precision floating-point numbers can represent values ranging from approximately $2^{-127} \approx 6 \times 10^{-39}$  to $2^{128} \approx 3 \times 10^{38}$.
These correspond to the `float` type in C or `np.single` in Python’s NumPy library.

+++

### Double-Precision Floating Point

```{figure} figures/double.svg
---
height: 120px
---
Memeory layout of a double precision float.
```

For double precision (64-bit), the exponent is stored with a bias of 1023 and the significand with an implicit leading 1 is stored in 52 bits.
The exponent can represent values from -1023 to +1024.
Double-precision numbers can represent values ranging from approximately $10^{-308}$ to  $10^{308}$, corresponding to the `double` type in C or `np.double` in Python/NumPy.

+++

### Other Floating Point

* "half percision" `float16`
* `bfloat16`, used for neural network
* `long double`, could be 80-bit or 128-bit, dependent on the system.

### Encoding of Special Values

```{code}
val    s_exponent_signcnd
+inf = 0_11111111_0000000
-inf = 1_11111111_0000000
```

```{code}
val    s_exponent_signcnd
+NaN = 0_11111111_klmnopq
-NaN = 1_11111111_klmnopq
```
where at least one of `k`, `l`, `m`, `n`, `o`, `p`, or `q` is 1.

### NaN Comparison Rules

```c
/* To run me, save into "demo.c" and then compile with `gcc demo.c -lm -o demo` */

#include <stdio.h>

int
main(int argc, char *argv[])
{
	float x = 0.0 / 0.0;

	if (x == 0.0 / 0.0)
		printf("NaN\n");
	else
		printf("Total confused...\n");

	return 0;
}
```

In C, `NaN` (Not a Number) has some special comparison rules:
* Comparing `NaN` with anything always returns false: this means that `x == NaN`, `x != NaN`, `x < NaN`, `x > NaN`, `x <= NaN`, and `x >= NaN` are all false, regardless of the value of `x` (even if `x` is also `NaN`).
* Use `isnan()` to check for `NaN`: the standard C library provides the `isnan()` function in `math.h` to check if a floating-point value is `NaN`.

+++

```{note}

### Data type conversion errors can be very costly

Ariane 5 is the primary launch vehicle of the European Space Agency (ESA) that operates from the Guiana Space Center near Kourou in the French Guiana.
Its first successful operational flight took place in December 1999, when it carried to space the European X-ray Multi Mirror (XMM) satellite.
Its first test launch, however, on June 4, 1996 resulted in failure, with the rocket exploding 40 seconds into the flight.

A study of the accident by an inquiry board found that the failure was caused by a software problem in the Inertial Reference System that was guiding the rocket.
In particular, the computer program, which was written in the Ada programming language and was inherited from the previous launch vehicle Ariane 4, required at some point in the calculation a conversion of a 64-bit floating point number to a 16-bit integer.
The initial trajectory of the Ariane 5 launch vehicle, however, is significantly different than that of Ariane 4 and the 16 bits are not enough to store some of the required information.
The result was an error in the calculation, which the inertial system misinterpreted and caused the rocket to veer off its flight path and explode.
The cost of the failed launch was upwards of 100 million dollars!
```

+++

## Machine Accuracy

In order to quantify truncation errors, we define:
\begin{align}
\mbox{(relative error)} \equiv \frac{x - \bar{x}}{x}.
\end{align}

If we use a numeral system of base b and keep p significant digits, the machine accuracy is
\begin{align}
\epsilon = \left(\frac{b}{2}\right) b^{-p}.
\end{align}

A single-precision floating-point number, which stores 23 significant digits in binary (the mantissa), provides a machine accuracy of approximately $\epsilon_\mathrm{single} = 2^{-23} \approx 10^{-7}$ in decimal.
In contrast, a double-precision floating-point number, with 52 significant binary digits, corresponds to a much finer machine accuracy of about $\epsilon_\mathrm{double} = 2^{-52} \approx 2\times10^{-16}$ in decimal.

+++

## Another Moment of ZEN

We all learned in high school that the solutions (roots) to the qudratic equation $a x^2 + b x + c = 0$ is
\begin{align}
x = \frac{-b \pm \sqrt{b^2 - 4 a c}}{2a}.
\end{align}

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
\begin{align}
x \approx -\frac{b}{a} \mbox{ or } -\frac{c}{b}
\end{align}
Is it possible to recover the small root $-c/b$?

+++

When $b > 0$, a catastrophic cancellation (see below) happens only in the "+" equation.
We may replace the first qudratic equation by its "conjugate" form
\begin{align}
x = \frac{2c}{-b \mp \sqrt{b^2 - 4 a c}}
\end{align}

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

## Catastrophic Cancellation

Catastrophic cancellation occurs in numerical computing when subtracting two nearly equal numbers, leading to significant loss of precision.
This issue arises because when two close numbers are subtracted, the leading digits cancel out, and the result retains only the less significant digits.
In floating-point arithmetic, these less significant digits may already have been affected by rounding errors, leading to a final result with much lower accuracy than expected.

For example, consider subtracting $x = 1.00000001$ and $y = 1.00000000$ in floating-point arithmetic.
Although the true result is $0.00000001$, if both numbers are rounded to eighth significant digits during storage (e.g., in single precision float), they might both be stored as $1.000000$.
Subtracting these will result in $0$, entirely losing the small difference.

This effect is particularly problematic in computations involving functions where nearly equal terms naturally occur (such as in the calculation of certain derivatives or in algorithms for solving linear systems).
Techniques like reformulating equations to avoid such subtractions or using higher precision arithmetic can help mitigate catastrophic cancellation, as we seen in the qudratic equation solver example above.

+++

## Fast Inverse Square Root

One fascinating example of taking advantage of floating-point representation is the fast inverse square root algorithm, famously used in computer graphics, particularly in the development of 3D engines like those in video games.
This algorithm provides an efficient way to compute $1/\sqrt{x}$, a common operation in tasks like normalizing vectors.
Instead of relying on the standard floating-point math libraries, the algorithm "hacks" the floating-point representation by manipulating the bits directly to achieve a remarkably fast approximation of the result.
This method takes advantage of how floating-point numbers are stored, using a clever combination of bitwise operations and mathematical magic to drastically speed up computation.

The following code is the fast inverse square root implementation from Quake III Arena, stripped of C preprocessor directives, but including [modified comment text](https://web.archive.org/web/20170729072505/https://github.com/id-Software/Quake-III-Arena/blob/master/code/game/q_math.c#L552):
```c
float Q_rsqrt( float number )
{
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;                     // evil floating point bit level hacking
	i  = 0x5f3759df - ( i >> 1 );             // what the f**k?
	y  = * ( float * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) ); // 1st iteration
//	y  = y * ( threehalfs - ( x2 * y * y ) ); // 2nd iteration, this can be removed

	return y;
}
```

+++

## `Arepo`

```{figure} figures/arepo1.png
---
height: 240px
---
The moving mesh code Arepo.
```

```{figure} figures/arepo2.png
---
height: 480px
align: right
---
The fast maping between 53-bit integer and double precision float.
```

The Arepo cosmology code, widely used in astrophysical simulations, employs a moving mesh method to solve fluid dynamics problems with high accuracy and efficiency.
To avoid numerical issues commonly encountered in floating-point arithmetic, Arepo implements an innovative approach to ensure robustness in geometric computations.
Specifically, it maps the significant digits of double-precision floating-point numbers to a 53-bit integer, which allows for exact arithmetic when evaluating geometric predicates.
By doing this, Arepo avoids round-off errors that could otherwise compromise the accuracy of simulations.
This precise handling of numerical data is critical for the accurate modeling of complex cosmological structures.

+++

## References

* [Preliminary discussion of the logical design of an electronic computing instrument](https://www.ias.edu/sites/default/files/library/Prelim_Disc_Logical_Design.pdf) by Burks, Goldstine, and von Neumann (1964)
* [What Every Computer Scientist Should Know About Floating-Point Arithmetic](https://docs.oracle.com/cd/E19957-01/800-7895/800-7895.pdf) by D. Goldberg, in the March 1991 of Computing Surveys
