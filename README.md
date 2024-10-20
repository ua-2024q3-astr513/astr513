# ASTR 513: Computation and Statistical Methods


This is a
[Jupyter Book](https://jupyterbook.org/) for
[University of Arizona](https://www.arizona.edu/)'s
[ASTR 513: Computational and Statistical Methods for Astrophysics](https://ua-2024q3-astr513.github.io/astr513/),
taught in Fall 2024.
It covers the "computation" part of the course.


## Contributing

All materials in this Jupyter Book are stored a special flavour of
Markdown called
[MyST (or Markedly Structured Text)](https://myst-parser.readthedocs.io/).

To edit the materials locally as Jupyter Notebooks, simply use
Jupytext to create a pair notebook:
```
jupytext --sync [chapter].md
```
All changed made the notebook will be automatically synced to the MyST
markdown file.

To build the Jupyter Book locally, run
```
git clone git@github.com:ua-2024q3-astr513/astr513.git
cd astr513
jupyter-book build docs
```
