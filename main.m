%{
    Demonstration of IntegrateOverCurvedElement

    ***********************************************************************

    Copyright (C) 2021 Jeffrey S. Ovall, Samuel E. Reynolds

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Contact:
        Jeff Ovall:     jovall@pdx.edu
        Sam Reynolds:   ser6@pdx.edu 

    ***********************************************************************
    
    This code is a demonstration of the numerical method described in 
    "Quadrature for Implicitly-defined Finite Element Functions on 
    Curvilinear Polygons" by J. Ovall, S. Reynolds, submitted to SIAM J. 
    of Sci. Comp. (2021).

    The goal is to compute the L^2(K) inner product and H^1(K) semi-inner 
    product of functions v and w, where K \subset R^2 is a curvilinear
    polygon, and v,w \in H^1(K) have polynomial Laplacians and piecewise
    polynomial Dirichlet traces. 

    The code is split among four files:
        IntegrateOverCurvedElement.m
        class_element.m
        class_poly.m
        class_quad.m
    The first file is the primary function used to execute the quadratures
    considered in our paper, and the other three are supporting code that
    define classes used in the primary file. A fifth file, BuildElement.m, 
    contains code to construct the domains considered in our paper, but is
    not necessary for the primary algorithm to work.

    This file provides a simple demonstration of the algorithm. The default
    settings consider the unit square K = (0,1) x (0,1) with 

        \Delta v = x^2*y in K ,     v = 0 on \partial K , 
        \Delta w = y^2 in K ,       w = 0 on \partial K .

    In the notation of the paper, we have the multi-indices 

        \alpha = (2,1) , \beta = (0,2) .

    This script approximates the quantities 

        L2 = \int_K v*w dx = 8.101386165180633e-05 ,
        H1 = \int_K (\nabla v)*(\nabla w) dx = 1.905102279276017e-03 ,

    by using the quadrature method described in the paper.

    The Kress quadrature parameter n may be increased to improve the
    accuracy of the calculation, at a greater computational cost.

    Different Laplacians for v and w can be obtained by changing the
    polynomial objects p and q, respectively. For instance,     

        p_idx = [3,5];
        p_coef = [2,-7];
        p=class_poly;
        p=p.init(p_idx,p_coef);

    will define the Laplacian of v to be p(x,y) = 2*x^2 - 7*y^2 . To 
    verify the definition of p, one may use p.print() to get

        >> p.print()
         + 2(x^2) - 7(y^2)

    The Laplacian of w is given by q. See eqn. (5) and Tables 1 and 7 in 
    our paper to see how we enumerate the multi-indices. See also the
    commments in class_poly.m for more info on how polynomial objects are
    used in this code.
    
    The Dirichlet traces of v and w are given by f and g, respectively. In
    the example demonstrated here, we choose f=g=0 for the sake of
    simplicity. The user is free to change them, but be aware that the
    tangential derivatives will need to change accordingly. In particular,
    if the boundary \partial K is parameterized by x(t) and F(t) = f(x(t))
    is the provided Dirichlet data, then it holds that the tangential 
    derivative is given by F'(t)/|x'(t)|. 

    Different domains can be selected by changing the string element_name.
    Supported options are:
        'square'    (the unit square (0,1) x (0,1))
        'circle'    (the unit circle, a bigon)
        'pacman'    (the unit circle with pizza sliced removed)
        'shuriken'  (a perturbation of the unit square with sinusoidal edges)
        'jigsaw'    (a puzzle piece)
    Alternatively, the user may constuct their own domain. See the comments
    in class_element.m and BuildElement.m on how to do this.

    Please see the comments made in the supporting files for more info on
    the specifics of the code. For a derivation of the quadrature method,
    please refer to our paper.

    Any future updates and addtional demos will be posted to our GitHub
    repository found at 

    https://github.com/samreynoldsmath/HigherOrderCurvedElementQuadrature

%}

clear
close all
format long e

%% Choose element type
element_name='square';
options=[];

%% Kress quadrature
n=16;       % 2*n points sampled per edge
sigma=7;    % controls distribution of sampled points

QUAD=class_quad;
QUAD=QUAD.KressWeights(n,sigma);

%% Construct element boundary
ELEMENT=BuildElement(element_name,QUAD,options);

%% Define z, as used in the shifted monomial basis {(x-z)^alpha}
z=[0,0];

%% Total number of sampled boundary points 
N=2*ELEMENT.num_edges*QUAD.n;

%% Implicitly define v

% Laplacian of v
p=class_poly;

p_idx=[7];      % Integer index of multi-index \alpha
p_coef=[1];     % Coefficient on (x-z)^\alpha

p=p.init(p_idx,p_coef);

% Boundary data for v
f=zeros(N,1);       % Dirichlet trace
v_td=zeros(N,1);    % Tangential derivative

%% Implicitly define w

% Laplacian of w
q=class_poly;

q_idx=[5];
q_coef=[1];

q=q.init(q_idx,q_coef);

% Boundary data for w
g=zeros(N,1);
w_td=zeros(N,1);

%% Use quadrature method to obtain L2, H1 inner products 
[L2,H1]=IntegrateOverCurvedElement(ELEMENT,QUAD,z,p,f,v_td,q,g,w_td)
