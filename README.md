# Density Gradient Example

## References

This repository implements the Density Gradient method in a way similar to this paper:


[https://doi.org/10.1080/1065514021000012363](https://doi.org/10.1080/1065514021000012363)
```
@ARTICLE{WettsteinVLSI2002,
author={Andreas Wettstein and Oleg Penzin and Eugeny Lyumkis},
title={Integration of the Density Gradient Model into a General Purpose Device Simulator},
journal={VLSI Design},
volume={15},
number={4},
pages={751--759},
year={2002},
doi={10.1080/1065514021000012363},
}
```

With the insulator boundary condition described in:

[https://doi.org/10.1109/TCAD.2011.2107990](https://doi.org/10.1109/TCAD.2011.2107990)
```
@ARTICLE{GarciaAsenov2011,
author={Garcia-Loureiro, A.J. and Seoane, N. and Aldegunde, M. and Valin, R. and Asenov, A. and Martinez, A. and Kalna, K.},
journal=ieeetcad,
title={Implementation of the Density Gradient Quantum Corrections for 3-D Simulations of Multigate Nanoscaled Transistors},
year={2011},
month=jun,
volume={30},
number={6},
pages={841--851},
doi={10.1109/TCAD.2011.2107990},
ISSN={0278-0070},}
```

## Notes

While there is code for full drift diffusion, it has only been really tested for MOSCAP simulation with coupling of the Potential equation with either the electron or hole density equation.

A description of the derivation of this model is in this document:
[https://github.com/devsim/devsim_misc/blob/main/devsim_docs/TCADdocs.pdf](https://github.com/devsim/devsim_misc/blob/main/devsim_docs/TCADdocs.pdf)

* `test_1D.py` - 1D simulation plots
* `runs.sh` - Runs 1D plots for different doping and oxide thicknesses
* `moscap2d.geo` - 3nm tox from gmsh
* `test_2d.py` - 2d example

