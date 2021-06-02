%{
    IntegrateOverCurvedElement

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

    This function computes the integrals

        \int_K v*w dx  ~and~  \int_K (\nabla v)*(\nabla w) dx

    only in terms of boundary integrals, where K is a curvilinear polygon,
    and v,w are functions in H^1(K) satisfying

        \Delta v = p in K , v = f on \partial K , 
        \Delta w = q in K , w = g on \partial K .

    We assume p,q : K -> R are polynomials (of degree at most m-2) and 
    f,g : \partial K -> R are continuous, and when restricted to any edge
    of \partial K, it holds that f,g are traces of polynomials of degree 
    at most m.
    
    Method is adapted from "Quadrature for Implicitly-defined Finite
    Element Functions on Curvilinear Polygons" by J. Ovall, S. Reynolds,
    submitted to SIAM J. of Sci. Comp. (2021).

    The Dirichlet-to-Neumann map is adapted from "A High-order Method for
    Evaluating Derivatives of Harmonic Functions in Planar Domains" by J.
    Ovall, S. Reynolds, SIAM J. of Sci. Comp. (2018).

    DEPENDENCIES:
        class_element.m
        class_poly.m
        class_quad.m

    INPUT:
        ELEMENT is an object of type class_element, a custom data structure
            that contains relevant data about the curvilinear polygon K.
                * See comments in class_element.m
        QUAD is an object of type class_quad, which contains the quadture
            weights and parameter samplings on [0,2*pi] for the 
            parameterization of the edges of K. 
                * The boundary is traversed counterclockwise
                * Each edge is sampled at 2n+1 points, including the
                endpoints. To avoid redundancy, only the first endpoint is
                included in calculations. Since the Kress quadrature
                weights are zero at the endpoints, the quadrature used in
                solving the Nystrom uses only 2n-1 points per edge. No
                interpolation is used with FFT, which also uses 2n points
                per edge.
                * See comments in class_quad.m
        z is a point in R^2, as used in the shifted monomial basis,
            { (x-z)^alpha : alpha is a multi-index }
        p,q are sparse polynomial objects (type class_poly) that specify
            the coefficients of the Laplacians of v and w, repectively, in
            the shifted monomial basis
                * See comments in class_poly.m 
        f,g are vectors containing the values of the traces of v and w, 
            respectively, on the boundary \partial K. The length of each
            vector is equal to the number of sampled boundary points:
            2*QUAD.n*ELEMENT.num_edges
        v_td,w_td are vectors containing the values of the tangential 
            derivatives of v and w, respectively, with the same length as 
            f and g.

    OUTPUT:
        L2 is the inner product \int_K v*w dx 
        H1 is the semi-inner product \int_K (\nabla v)*(\nabla w) dx

    POTENTIAL IMPROVEMENTS:
        * Compute coefficients of anti-Laplacians of polynomials offline
            (use a look-up table, for example)
        * Set up and solve Nystrom systems in parallel
        * Make option to compute only L2 or H1
        * Write simplified routines for v,w lying in one of V_m^K(K) or
            V_m^{\partial K}(K)
        * Automate computation of (weighted) tangential derivatives of 
            provided Dirichlet data
        * Singularity subtraction to improve performance on domains with 
            reentrant corners

    Last updated 6/2/2021

%}

function [L2,H1]=IntegrateOverCurvedElement(ELEMENT,QUAD,z,p,f,v_td,q,g,w_td)
    
    %% Coefficients of polynomial anti-Laplacians
    P = p.anti_laplacian_poly();
    Q = q.anti_laplacian_poly();
    
    Pstar = P.anti_laplacian_poly();
    Qstar = Q.anti_laplacian_poly();
    
    %% Coefficients of polynomial products
    pQ = prod_poly(p,Q);
    PQ = prod_poly(P,Q);
    
    %% Values of polynomials on the boundary
    Q_trace = PolynomialTrace(Q,z,ELEMENT);
    P_trace = PolynomialTrace(P,z,ELEMENT);
    
    Pstar_trace = PolynomialTrace(Pstar,z,ELEMENT);
    Qstar_trace = PolynomialTrace(Qstar,z,ELEMENT);
    
    fP_trace = f(:)-P_trace;
    gQ_trace = g(:)-Q_trace;
        
    %% Values of normal derivatives of polynomials
    P_nd = PolynomialNormalDerivative(P,z,ELEMENT);
    Pstar_nd = PolynomialNormalDerivative(Pstar,z,ELEMENT);
    Qstar_nd = PolynomialNormalDerivative(Qstar,z,ELEMENT);
    
    %% Values of tangential derivatives of P and Q
    P_td = PolynomialTangentialDerivative(P,z,ELEMENT);
    Q_td = PolynomialTangentialDerivative(Q,z,ELEMENT);
    
    %% Dirichlet-to-Neumann map for harmonic functions
    %{
        These normal derivatives, by construction, are "pre-weighted."
        See eqns (32) and (33) in "A High-order Method for Evaluating 
        Derivatives of Harmonic Functions in Planar Domains" by J. Ovall, 
        S. Reynolds, SIAM J. of Sci. Comp. (2018).
    %}
    
    % Martensen quadrature weights
    MART = class_quad;
    MART = MART.MartensenWeights(QUAD.n);
    
    % (pre-weighted) normal derivatives
    wQ_nd_wgt = Dirichlet2Neumann(w_td-Q_td,ELEMENT,QUAD,MART);
    v1_nd_wgt = Dirichlet2Neumann(v_td,ELEMENT,QUAD,MART);
    v2P_nd_wgt = Dirichlet2Neumann(-P_td,ELEMENT,QUAD,MART);
    
    %% Anti-Laplacian of phi = v-P 
    phi_td = v_td-P_td;
    [PHI,PHI_nd,phi_nd_wgt] = BigPHI(fP_trace,phi_td,ELEMENT,QUAD,MART);
    
    %% Integrate
    
    % Integration weights 
    h = pi/(QUAD.n*ELEMENT.num_edges); % for pre-weighted terms
    WGT = ELEMENT.dx_norm.*kron(ones(ELEMENT.num_edges,1),QUAD.wgt(:));
        
    % Integrate v*w
    L2 = IntegratePolynomial(PQ,z,ELEMENT,WGT)...
        -h*sum( (PHI+Pstar_trace).*wQ_nd_wgt + Qstar_trace.*phi_nd_wgt )...
        +sum(WGT.*( (PHI_nd+Pstar_nd).*gQ_trace + Qstar_nd.*fP_trace ));
        
    % Integrate \nabla v * \nabla w
    H1 = - IntegratePolynomial(pQ,z,ELEMENT,WGT)...
        +h*sum( v1_nd_wgt.*g(:) + v2P_nd_wgt.*Q_trace )...
        +sum(WGT.*P_nd.*Q_trace);
    
end

%% Integate a polynomial p over the element volume
function val=IntegratePolynomial(p,z,ELEMENT,WGT)
    %{
        Computes the volumetric integral \int_K p dx, where p is a
        polynomial. 

        This volumetric integral is reduced to a boundary integral via
            \int_K (x-z)^\alpha dx 
            = 1/(2+|\alpha|) \int_{\partial K} (x-z)^\alpha (x-z)*n ds.
        This is eqn (9) in "Quadrature for Implicitly-defined Finite
        Element Functions on Curvilinear Polygons"
    %}

    % precompute shift 
    xz=ELEMENT.x(:,1)-z(1);
    yz=ELEMENT.x(:,2)-z(2);
    xn=xz.*ELEMENT.unit_normal(:,1)+yz.*ELEMENT.unit_normal(:,2);
    
    % integrate each monomial term
    val=0;
    for i=1:p.nz
        alpha=p.multi_index(p.idx(i));
        xa=xz.^alpha(1);
        ya=yz.^alpha(2);
        val=val+p.coef(i)*sum(xa.*ya.*xn.*WGT)/(2+alpha(1)+alpha(2));
    end
    
end

%% Obtain boundary trace of a polynomial
function trace=PolynomialTrace(p,z,ELEMENT)
    %{
        Obtains the Dirichlet trace of the polynomial p on the boundary
        \partial K, returns vector of boundary values.
    %}

    % precompute shift
    xz=ELEMENT.x(:,1)-z(1);
    yz=ELEMENT.x(:,2)-z(2);
    
    %
    trace=zeros(size(xz));
    for i=1:p.nz
        alpha=p.multi_index(p.idx(i));
        xa=xz.^alpha(1);
        ya=yz.^alpha(2);
        trace=trace+p.coef(i)*xa.*ya;
    end
end

%% Obtain tangential derivative of a polynomial
function td=PolynomialTangentialDerivative(p,z,ELEMENT)
    %{
        Obtains the tangential derivative of the polynomial p on the 
        boundary \partial K, returns vector of these values.
    %}
    
    % precompute shift
    xz=ELEMENT.x(:,1)-z(1);
    yz=ELEMENT.x(:,2)-z(2);
    
    td=zeros(size(xz));
    
    % obtain coefficients of gradient
    [px,py]=p.grad_poly();
    
    % sum over coefficients of x component
    for i=1:px.nz
        alpha=px.multi_index(px.idx(i));
        xa=xz.^alpha(1);
        ya=yz.^alpha(2);
        td=td+px.coef(i)*xa.*ya.*ELEMENT.unit_tangent(:,1);
    end
    
    for i=1:py.nz
        alpha=px.multi_index(py.idx(i));
        xa=xz.^alpha(1);
        ya=yz.^alpha(2);
        td=td+py.coef(i)*xa.*ya.*ELEMENT.unit_tangent(:,2);
    end
    
end

%% Obtain normal derivative of a polynomial
function nd=PolynomialNormalDerivative(p,z,ELEMENT)
    %{
        Obtains the normal derivative of the polynomial p on the 
        boundary \partial K, returns vector of these values.
    %}

    % precompute
    xz=ELEMENT.x(:,1)-z(1);
    yz=ELEMENT.x(:,2)-z(2);
    
    nd=zeros(size(xz));
    
    % obtain coefficients of gradient
    [px,py]=p.grad_poly();
    
    % sum over coefficients of x component
    for i=1:px.nz
        alpha=px.multi_index(px.idx(i));
        xa=xz.^alpha(1);
        ya=yz.^alpha(2);
        nd=nd+px.coef(i)*xa.*ya.*ELEMENT.unit_normal(:,1);
    end
    
    % and the y component 
    for i=1:py.nz
        alpha=py.multi_index(py.idx(i));
        xa=xz.^alpha(1);
        ya=yz.^alpha(2);
        nd=nd+py.coef(i)*xa.*ya.*ELEMENT.unit_normal(:,2);
    end
    
end

%% Neumann-to-Dirichlet map
function u=Neumann2Dirichlet(u_nd,ELEMENT,QUAD,MART)
    %{
        Soltion u to the Neumann problem where u is harmonic and its normal
        derivative du/dn = u_nd is given on \partial K.
        We use a Nystrom method to solve the integral equation via GMRES to
        obtain the trace of u on the boundary \partial K.
    %}

    % GMRES parameters
    restart=20;
    tol=1e-12;
    maxit=10;
    
    % Copy Kress weights into longer vector for all edges of boundary
    n=QUAD.n;
    wgt=kron(ones(ELEMENT.num_edges,1),QUAD.wgt(:));
    
    % Set up and solve system
    B=NystromRHS(u_nd,ELEMENT,MART,QUAD);
    AFUN=@(X)NystromLHS(X,ELEMENT,n,wgt);
    u=gmres(AFUN,B,restart,tol,maxit);
    
end

%% Dirichlet-to-Neumann map
function u_nd=Dirichlet2Neumann(u_td,ELEMENT,QUAD,MART)
    %{
        Given the tangential derivative u_td of a harmonic function u, 
        returns the normal derivative u_nd. This is accomplished by 
        solving an integral equation via a Nystrom method to obtain the 
        trace of a harmonic conjuate v on the boundary, and then 
        differentiating with respect to the boundary parameter t to obtain 
        the tangential derivative of v. Since v is a harmonic conjugate of 
        u, we have that v_td = u_nd.
    %}
    
    % half the number of sampled bounary points
    N=ELEMENT.num_edges*QUAD.n;
    
    % Solve complementary Neumann problem for the harmonic conjugate
    v=Neumann2Dirichlet(-u_td,ELEMENT,QUAD,MART);
    
    % Differentiate the harmonic conjugate to get the normal derivative
    u_nd=FFT_derivative(v,N);
    
end

%% Right-hand side for the Nystrom system
function B=NystromRHS(u_nd,ELEMENT,MART,QUAD)    
    %{
        Computes right-hand side of the Nystrom system.
    %}

    % Tolerance for square norm
    TOL=1.0e-14;

    % For readability
    m=ELEMENT.num_edges;
    n=QUAD.n;
    
    % Allocate space for output
    B=zeros(2*m*n,1);
    
    % Precompute scaling, common product
    scale = -0.25/pi;
    npi=n/pi;
    
    % outer loop over edges
    for edgei=1:m
        % loop over points on edge i
        for i=1:2*n
            ii=i+(edgei-1)*(2*n);
            % inner loop over edges
            for edgej=1:m
                
                if edgei==edgej
                    
                    % Kress and Martensen quadratures
                    for j=2:2*n
                        
                        jj=j+(edgej-1)*(2*n);
                        ij=abs(i-j)+1;

                        xx=ELEMENT.x(ii,1)-ELEMENT.x(jj,1);
                        yy=ELEMENT.x(ii,2)-ELEMENT.x(jj,2);
                        xy2=xx*xx+yy*yy;

                        if xy2<TOL
                            L1=2*scale*log( npi*QUAD.wgt(j)...
                                *ELEMENT.dx_norm(jj) );
                        else
                            L1=scale*log( xy2/MART.t(ij) );
                        end

                        B(ii)=B(ii)+( L1+npi*MART.wgt(ij) )...
                            *ELEMENT.dx_norm(jj)*QUAD.wgt(j)*u_nd(jj);

                    end
                    
                else
                    
                    % Kress quadrature only
                    for j=2:2*n
                        
                        jj=j+(edgej-1)*(2*n);

                        xx=ELEMENT.x(ii,1)-ELEMENT.x(jj,1);
                        yy=ELEMENT.x(ii,2)-ELEMENT.x(jj,2);
                        xy2=xx*xx+yy*yy;

                        B(ii)=B(ii)+scale*log(xy2)...
                            *ELEMENT.dx_norm(jj)*QUAD.wgt(j)*u_nd(jj);

                    end
                    
                end
            end
        end    
    end
    
end

%% Operator for the Nystrom system
function Y=NystromLHS(V,ELEMENT,n,wgt)
    %{
        Left-hand side operator of the Nystrom system.
    %}
    
    % Tolerance for square norm
    TOL=1.0e-14;

    % For readability
    m=ELEMENT.num_edges;

    % Allocate space for output
    Y=zeros(2*m*n,1);
    
    % Precompute scaling, corner values
    scale = 0.5/pi;
    
    corner_idx=kron( 1:2*n:2*m*n, ones(1,2*n) );
    corner_idx=circshift(corner_idx,-n);
    vc=V(corner_idx);
    
    % Indices for summation: first and last terms on each edge neglected
    J=1:2*m*n;
    J(1:2*n:2*m*n)=[];
    
    % loop over points on boundary
    for i=1:2*m*n
        
        Y(i)=ELEMENT.arc_length*vc(i)+0.5*(V(i)-vc(i));
        
        % summation
        for j=J
            
            xx=ELEMENT.x(i,1)-ELEMENT.x(j,1);
            yy=ELEMENT.x(i,2)-ELEMENT.x(j,2);
            xy2=xx*xx+yy*yy;
            
            % Taylor limit for diagonal entries
            if xy2<TOL
                xx=ELEMENT.ddx(j,1);
                yy=ELEMENT.ddx(j,2);
                xy2=2*ELEMENT.dx_norm(j)*ELEMENT.dx_norm(j);
            end
            
            K=scale*(xx*ELEMENT.dx(j,2)-yy*ELEMENT.dx(j,1))/xy2;
            
            Y(i)=Y(i)+(K+ELEMENT.dx_norm(j))*wgt(j)*(V(j)-vc(i));
            
        end
    end
    
end

%% Differentiation via Fast Fourier Transform
function DF=FFT_derivative(F,N)
    %{
        Given time-series data F (periodic, equispaced in time), 
        return the derivative DF. F is sampled at 2*N points.

        Differentiation is performed by differentiation of a 
        truncated Fourier series, obtained with FFT.
    %}

    % Fourier coefficients
    alpha=fft(F); 
    
    % Differentiate with IFFT
    k=1i*[0:(N-1), 0, (-N+1):(-1)]; % scale coeficients by chain rule
    DF=real(ifft(k(:).*alpha(:)));
    
end

%% Anti-Laplacian of phi
function [PHI,PHI_nd,phi_nd_wgt]=BigPHI(phi,phi_td,ELEMENT,QUAD,MART)

    %{
        phi = v-P on the boundary, harmonic in the interior
        phi_td is the tangential derivative of phi
        phi_nd_wgt is the weighted normal derivative of phi
        PHI = (x_1 * rho + x_2 * rhohat)/4 where
            rho is harmonic
            rhohat is a harmonic conjugate of rho
            phihat is a harmonic conjugate of phi
            rho_nd = (phi,-phihat)*n
            rhohat_nd = (phihat,phi)*n
        By construction, it holds that
            \Delta PHI = phi
        See Prop. 2.2 in "Quadrature for Implicitly-Defined Finite Element
        Functions on Curvilinear Polygons"
    %}
    
    % obtain trace of the harmonic conjugate of phi 
    phihat=Neumann2Dirichlet(-phi_td,ELEMENT,QUAD,MART);
    
    % weighted normal derivative of phi
    N=ELEMENT.num_edges*QUAD.n;
    phi_nd_wgt=FFT_derivative(phihat,N);
    
    % compute normal derivatives of rho and its harmonic conjugate
    rho_nd = phi.*ELEMENT.unit_normal(:,1)...
        - phihat.*ELEMENT.unit_normal(:,2);
    rhohat_nd = phihat.*ELEMENT.unit_normal(:,1)...
        + phi.*ELEMENT.unit_normal(:,2);
        
    % solve Neumann problem to get trace of rho
    rho=Neumann2Dirichlet(rho_nd,ELEMENT,QUAD,MART);
    
    % solve Neumann problem to get trace of harmonic conjugate of rho
    rhohat=Neumann2Dirichlet(rhohat_nd,ELEMENT,QUAD,MART);
    
    % compute PHI and its normal derivative
    PHI=0.25*( ELEMENT.x(:,1).*rho + ELEMENT.x(:,2).*rhohat );
    PHI_nd=0.25*(...
        ELEMENT.x(:,1).*rho_nd...
        +ELEMENT.x(:,2).*rhohat_nd...
        +ELEMENT.unit_normal(:,1).*rho...
        +ELEMENT.unit_normal(:,2).*rhohat);
    
end