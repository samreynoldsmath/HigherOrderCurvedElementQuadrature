# HigherOrderCurvedElementQuadrature

This repository contains the MATLAB code used for the numerical examples in the paper "Quadrature for Implicitly-defined Finite Element Functions on Curvilinear Polygons" by Jeffrey S. Ovall and Samuel E. Reynolds, submitted to SIAM Journal of Scientific Computing.

The quadrature approximates the integrals \int_K vw dx and \int_K (\nabla v)(\nabla w) dx, where K \subset R^2 is a curvilinear polygon and v,w have polynomial Laplcians and continuous.Dirichlet traces on the boundary \partial K, such that the restriction to one edge is the trace of a polynomial. In our paper, we show that these volumetric integrals can be reduced to boundary integrals, and provide a corresponding quadrature, which is implemented here.

A demonstration of the quadrature is provided in main.m. Please see the comments in main.m on how to use the code.

Any future updates and additional demos will be posted to our GitHub repository: 
https://github.com/samreynoldsmath/HigherOrderCurvedElementQuadrature
