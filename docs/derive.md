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

# Numerical and Automatic Derivatives

+++

## Introduction to Derivatives

### Brief recap of derivatives in calculus

### Importance of derivatives in science and engineering

### Applications

+++

## Symbolic Differentiation

+++

## Numerical Differentiation

### Finite Difference Methods

* Low-Order Finite Difference Formulas
  * Forward difference formula.
  * Backward difference formula.
  * Central difference formula.
  * Error analysis and truncation errors.

* High-Order Finite Difference Methods
  * Derivation of higher-order formulas.
  * Error reduction and convergence rates.

### Spectral Methods and Fourier Transform

* Connection between finite differences and Fourier Transform.
* Introduction to spectral differentiation.
* Advantages for smooth periodic functions.

### Complex Step Differentiation

* Introduction to Complex Step Method
* Concept and mathematical foundation.
* Elimination of subtractive cancellation errors.

+++

## Automatic Differentiation

### Introduction to Automatic Differentiation

* Motivation and need for AD.
* Differences between AD, symbolic differentiation, and numerical approximation.

### Dual Numbers and Forward Mode AD

* Understanding Dual Numbers
  * Definition and algebra of dual numbers.
  * How dual numbers enable forward mode AD.
* Connection to Complex Step Method
  * Mathematical similarities and practical differences.

### Reverse Mode AD and Backpropagation

* Concept of Reverse Mode AD
  * Computational graphs and the chain rule.
* Backpropagation in Neural Networks
  * Application in machine learning.

### Higher-Order Derivatives: Jacobians and Hessians

* Importance of higher-order derivatives
  * Applications in optimization algorithms.
  * Computing Jacobians and Hessians with AD
* Techniques and computational considerations.

### Limitations and Challenges of AD

* Handling non-differentiable functions and discontinuities.
* Computational overhead and memory usage.
* Control flow and dynamic graphs.
* Best practices for efficient AD implementation.

+++

## Comparison of Differentiation Methods

### Numerical vs. Automatic vs. Symbolic Differentiation

* Strengths and weaknesses of each method.

### Choosing the right method for different applications.

+++

## Summary

### Open questions and discussion

### Suggestions for further reading and exploration.

```{code-cell} ipython3

```
