%{
    class_element

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

    Defines object that contains information about the parameterization of
    the boundary of an element K \subset R^2. Below, we use N to denote the
    number of sampled points, which is equal to
    N=2*QUAD.n*ELEMENT.num_edges.
        
        num_edges (positive integer):  number of edges
        x (N by 2 float array): physical points sampled on the boundary, 
            x(:,1) are x_1 coordinates
            x(:,2) are x_2 coordinates
        dx (same size as x): the derivative of x with respect to the
            parameter t, same coordinate slicing as x
        ddx (same size and slicing as x): the second derivative
        dx_norm (N by 1 float): the Euclidean norm of the derivative
        unit_tangent (same size as x): the unit tangent vector
        unit_normal (same size as x): the unit normal vector
        arc_length (float): the arc length of the boundary

    The boundary must be traversed counterclockwise. The user is 
    responsible for this, and there is no check of the orientation.

    An element may be created by the user by specifying the number of 
    edges, the x,y coordinates of the parameterized boundary, as well as 
    the first and second derivatives, and by providing the associated 
    quadrature weights. This may be accomplished by using the init method, 
    which automatically computes dx_norm, unit_tangent, unit_normal, and 
    arc_length (which is computed via the trapezoid rule).

    See BuildElement.m for examples of how to use the init method to
    initialize an element object.
    
%}

classdef class_element
    
    properties
        num_edges
        x
        dx
        ddx
        dx_norm
        unit_tangent
        unit_normal
        arc_length
    end
    
    methods
        function obj=init(obj,num_edges,x,y,dx,dy,ddx,ddy,QUAD)
            
            obj.num_edges=num_edges;
            
            obj.x=[x(:),y(:)];
            obj.dx=[dx(:),dy(:)];
            obj.ddx=[ddx(:),ddy(:)];
            
            % derivative magnitude and outward unit normal
            obj.dx_norm=sqrt(sum(obj.dx.*obj.dx,2));
            obj.unit_tangent=obj.dx./obj.dx_norm;
            
            obj.unit_normal=[...
                obj.unit_tangent(:,2),...
                -obj.unit_tangent(:,1)];
            
            % boundary length by trapezoid rule 
            obj.arc_length=0;
            N=2*QUAD.n;
            
            for k=1:obj.num_edges
                idx=(k-1)*N+1:k*N;
                obj.arc_length=obj.arc_length...
                    +sum(QUAD.wgt(:).*obj.dx_norm(idx));
            end
            
        end
    end
    
end