# On Affordable High-Order Entropy-Conservative/Stable and Well-Balanced Methods for Nonconservative Hyperbolic Systems

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17106781.svg)](https://zenodo.org/doi/TODO)
DOI TODO

This repository contains information and code to reproduce the results presented in the article

```bibtex
@online{artiano2026nonconservative,
  title={{O}n Affordable High-Order Entropy-Conservative/Stable
         and Well-Balanced Methods for Nonconservative Hyperbolic Systems},
  author={Artiano, Marco and Ranocha, Hendrik},
  year={2026},
  month={TODO},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={TODO}
}
```

If you find these results useful, please cite the article mentioned above. If you use the implementations provided here, please also cite this repository as

```bibtex
@misc{artiano2026nonconservativeRepo,
  title={Reproducibility repository for
         "{O}n Affordable High-Order Entropy-Conservative/Stable and
           Well-Balanced Methods for Nonconservative Hyperbolic Systems"},
  author={Artiano, Marco and Ranocha, Hendrik},
  year={2026},
  howpublished={\url{TODO}},
  doi={TODO}
}
```

## Abstract
Many entropy-conservative and entropy-stable (summarized as entropy-preserving) methods for hyperbolic conservation laws rely on Tadmor's theory for two-point entropy-preserving numerical fluxes and its higher-order extension via flux differencing using summation-by-parts (SBP) operators, e.g., in discontinuous Galerkin spectral element methods (DGSEMs).
The underlying two-point formulations have been extended to nonconservative systems using fluctuations by Castro et al.\ (2013, \href{https://doi.org/doi/10.1137/110845379}{doi:10.1137/110845379}) with follow-up generalizations to SBP methods.
We propose specific forms of entropy-preserving fluctuations for nonconservative hyperbolic systems that are simple to interpret and allow an algorithmic construction of entropy-preserving methods.
We analyze necessary and sufficient conditions, and obtain a full characterization of entropy-preserving three-point methods within the finite volume framework.
This formulation is extended to SBP methods in multiple space dimensions on Cartesian and curvilinear meshes.
Additional properties such as well-balancedness extend naturally from the underlying finite volume method to the SBP framework.
We use the algorithmic construction enabled by the chosen formulation to derive several new entropy-preserving schemes for nonconservative hyperbolic systems, e.g., the compressible Euler equations of an ideal gas using the internal energy equation and a dispersive shallow-water model.
Numerical experiments show the robustness and accuracy of the proposed schemes.

## Numerical experiments
To reproduce the numerical experiments presented in this article, you need to install Julia. The numerical experiments presented in this article were performed using Julia v1.10.6.

First, you need to download this repository, e.g., by cloning it with git or by downloading an archive via the GitHub interface. Then, you need to start Julia in the code directory of this repository and follow the instructions described in the README.md file therein.

## Authors
- Marco Artiano
- [Hendrik Ranocha](https://ranocha.de/) (Johannes Gutenberg University Mainz, Germany)

## License

The code in this repository is published under the MIT license, see the `LICENSE` file.

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
