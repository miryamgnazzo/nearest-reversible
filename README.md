# nearest-reversible
This repository contains the code for computing the numerical approximation of the nearest reversible Markov Chain, to a given one. The approach relies on a Riemannian optimization problem, over the manifold

${S \in \mathbb{R}^{n\times n} : S > 0,\; S = S^\top,\; S \hat{{\pi}} = \hat{{\pi}}}$

## Main functions:
* **multinomialsymmetricfixedfactory.m**: contains an implementation of the manifold;
* **riemannian_nearest_reversible.m**: computes the nearest reversible matrix to a given one via Riemannian optimization.

## Dependencies:
* The use of the Riemannian optimization-based approach requires the [manopt package](https://www.manopt.org/index.html);
* The function **markov_generator.m** employs Python functions.
