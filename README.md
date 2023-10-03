# Higher Order Curved Element Quadrature

## Notice: This Repository is No Longer Maintained

**Important:** This repository is no longer actively maintained. It has been archived due to project completion and has been replaced by our new repository.

For the latest version and active development, please visit our new repository: [PuncturedFEM](https://github.com/samreynoldsmath/PuncturedFEM)


## Description
This repository contains the MATLAB code used for the numerical examples in the paper
Jeffrey S. Ovall, and Samuel E. Reynolds, "Quadrature for Implicitly-defined Finite Element Functions on Curvilinear Polygons," Computers & Mathematics with Applications (2022), Vol. 107 (1), pp. 1–16.
[https://doi.org/10.1016/j.camwa.2021.12.003](https://doi.org/10.1016/j.camwa.2021.12.003)

The quadrature approximates the integrals $\int_K v\,w\;dx$ and $\int_K \nabla v\cdot\nabla w\;dx$, where $K \subset \mathbb{R}^2$ is a curvilinear polygon and $v,w$ have polynomial Laplacians and have continuous Dirichlet traces on the boundary $\partial K$, such that the restriction to one edge is the trace of a polynomial. In our paper, we show that these volumetric integrals can be reduced to boundary integrals, and provide a corresponding quadrature, which is implemented here.

A demonstration of the quadrature is provided in [main.m](main.m). Please see the comments in main.m on how to use the code.

We ultimately plan to use this quadrature method to implement an efficient, high-order FEM on curvilinear meshes in 2-d domains. Our previous work demonstrated that the lowest order case is feasible, see "Trefftz Finite Elements on Curvilinear Polygons" by A. Anand, J. Ovall, S. Reynolds, S. Weisser, SIAM J. of Sci. Comp. (2018).


## License
Copyright (C) 2021 Jeffrey S. Ovall and Samuel E. Reynolds.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
