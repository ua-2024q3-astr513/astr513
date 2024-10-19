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
people think

-- Blaise Pascal (1623--1662)
```

+++

## Your Moment of ZEN

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

```{code-cell} ipython3

```
