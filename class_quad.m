%{
    class_quad

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

    Defines a quadrature object, used for Kress and Martensen quadratures.
    
    Each type of quadrature samples the boundary at 2n+1 points per edge,
    including the endpoints. Since each endpoint is shared by exactly two
    edges, to avoid redundancy we associate the first endpoint (proceeding
    counterclockwise around the boundary) with the current edge, and the
    second endpoint with the next edge, thus we associate only 2n points
    with each edge. 

    For a description of these quadratures, see "A High-order Method for 
    Evaluating Derivatives of Harmonic Functions in Planar Domains" by 
    Ovall and Reynolds, SISC 2018. Eqns (28) and (29) define the Kress
    quadrature, and eqns (17) and (18) for the Martensen quadrature.

%}

classdef class_quad
    
    properties
        n % 2*n sampled points per edge, including first point 
        t % sampled points
        wgt % weights 
    end
    
    methods
        %% Kress quadrature weights 
        function obj=KressWeights(obj,n,sigma)

            s=-1:(1/n):1; 

            c=(0.5-1/sigma)*s.^3+s/sigma+0.5; 
            d=fliplr(c); 
            denom=c.^sigma+d.^sigma;
            
            obj.n=n;
            
            % parameter sampling
            obj.t=2*pi*c.^sigma./denom;
            
            % weights
            obj.wgt=pi/n*(3*(sigma-2)*s.^2+2).*((c.*d).^(sigma-1))./(denom.^2);  
        
            % Delete repeated vertex
            obj.t(2*n+1)=[];
            obj.wgt(2*n+1)=[];
            
        end
        
        %% Martensen quadrature weights
        function obj=MartensenWeights(obj,n)
            
            h=pi/n;
            T=(0:h:2*pi-h)';
            
            % term used to split the right-hand side, see (13) and (14)
            obj.t = 4*sin(T/2).^2;
            
            % weights 
            obj.wgt=zeros(2*n,1);
            for m=1:n
                obj.wgt = obj.wgt + cos(m*T)/m;
            end
            obj.wgt = obj.wgt/(2*n); 
            
        end
        
    end
end