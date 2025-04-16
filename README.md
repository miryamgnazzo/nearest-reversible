# nearest-reversible
This repository contains the code for computing the numerical approximation of the nearest reversible Markov Chain, to a given one. The approach relies on a Riemannian optimization problem, over the manifold

$$ \mathcal{M}_{\pi} = \lbrace S \in \mathbb{R}^{n\times n} : S > 0, S = S^\top, S \hat{{\pi}} = \hat{{\pi}}, \rbrace$$

where $\pi$ is the stationary distribution associated with the original chain, and $\hat{\pi}$ is the vector containing the square roots of the entries of $\pi$, i.e., $\hat{\pi}_{i} = \sqrt{\pi_i}$ for $i=1,\ldots,n$.

> [!NOTE]
> If you use the code contained here and based on Riemannian optimization, please cite the relevant paper:
> - F. Durastante, M. Gnazzo, B. Meini. A Riemannian Optimization Approach for Finding the Nearest
>   Reversible Markov Chain

## Collaborators
* Miryam Gnazzo [📧](mailto:miryam.gnazzo@dm.unipi.it),
* Fabio Durastante [📧](mailto:fabio.durastante@unipi.it) [💻](https://github.com/Cirdans-Home),
* Beatrice Meini [📧](mailto:beatrice.meini@unipi.it).

## Main functions:
* `multinomialsymmetricfixedfactory.m`: contains an implementation of the manifold;
* `riemannian_nearest_reversible.m`: computes the nearest reversible matrix to a given one via Riemannian optimization.

## Dependencies and third party codes:
* The use of the Riemannian optimization-based approach requires the [manopt package](https://www.manopt.org/index.html);
* The function `markov_generator.m` employs Python functions which use the [NetworkX library](https://networkx.org/);
* The function `getClosestSparse.m` is the implementation from [iwhasherefirst2](https://github.com/iwasherefirst2/closest-reversible-markov-chain) implemeting the algorithm from
  * Nielsen, A. J. N., and Marcus Weber. "Computing the nearest reversible Markov chain." Numerical Linear Algebra with Applications 22.3 (2015): 483-499.
> [!IMPORTANT]
> If you end up using the `getClosestSparse.m` routine please cite the relevant paper:
> ```bibtex
> @article {MR3338930,
>    AUTHOR = {Nielsen, A. J. N. and Weber, M.},
>     TITLE = {Computing the nearest reversible {M}arkov chain},
>   JOURNAL = {Numer. Linear Algebra Appl.},
>  FJOURNAL = {Numerical Linear Algebra with Applications},
>    VOLUME = {22},
>      YEAR = {2015},
>    NUMBER = {3},
>     PAGES = {483--499},
>      ISSN = {1070-5325,1099-1506},
>   MRCLASS = {65C40 (65F35)},
>  MRNUMBER = {3338930},
>MRREVIEWER = {Myron\ Hlynka},
>       DOI = {10.1002/nla.1967},
>       URL = {https://doi.org/10.1002/nla.1967},
>}
> ```
* The function `getClosestSparse_gurobi.m` is a modified version of the original code in which Matlab's own quadratic programming solver `quadprog` has been subsistuted by the commercial [Gurobi](https://www.gurobi.com/) solver. This can be obtained with an academic license for research purposes for free. 
