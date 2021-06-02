%{ 
    BuildElement

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

    Defines the element K from predefined list. The boundary is traversed
    counterclockwise, with each edge parameterized on [0,2*pi]. 

    ELEMENT is an object of type class_element containing info about the
        parameterization of the boundary of K. See comments in 
        class_element.m for more info.

    element_name is a string specifying which type of domain to construct. 
    Supported options:
        'square'    (the unit square (0,1) x (0,1))
        'circle'    (the unit circle, a bigon)
        'pacman'    (the unit circle with pizza sliced removed)
        'shuriken'  (a perturbation of the unit square with sinusoidal edges)
        'jigsaw'    (a puzzle piece)

    QUAD is an object of type class_quad containing the quadrature weights
        and sampled values of the parameter t \in [0,2*pi].
 
    options is a struct containing any additional parameters needed to
        build a specific type of domain.

%}

%% construct element object
function ELEMENT=BuildElement(element_name,QUAD,options) 
    switch element_name
        case 'square'
            ELEMENT=Square(QUAD);
        case 'circle'
            ELEMENT=Circle(QUAD);
        case 'pacman'
            ELEMENT=PacMan(QUAD,options.mu);
        case 'shuriken'
            ELEMENT=Shuriken(QUAD,options.amp);
        case 'jigsaw'
            r=0.22; % tab radius
            b=0.17; % tab offset
            ELEMENT=Jigsaw(QUAD,r,b);
        otherwise
            error('Element type "%s" not supported (yet).',element_name)
    end
end

%% 
function ELEMENT=Square(QUAD)
    
    % define number of sampled points on boundary
    num_edges=4;
    N=2*QUAD.n;
    
    % precompute edges
    T=QUAD.t(:)/(2*pi);
    DT=ones(N,1)/(2*pi);
    Z=zeros(N,1);
    E=ones(N,1);
    
    % parameterization of boundary
    x=[ T; E; 1-T; Z ];
    y=[ Z; T; E; 1-T ];
    
    % first derivative
    dx=[ DT; Z; -DT; Z ];
    dy=[ Z; DT; Z; -DT ];
     
    % second derivative
    ddx=zeros(num_edges*N,1);
    ddy=zeros(num_edges*N,1);
    
    % construct element object
    ELEMENT = class_element;
    ELEMENT=ELEMENT.init(num_edges,x,y,dx,dy,ddx,ddy,QUAD);
    
    % overwrite boundary length to exact value
    ELEMENT.arc_length=4;
    
end

%%
function ELEMENT=Circle(QUAD)
    
    % define number of edges
    num_edges=2;
    
    % precompute edges
    C=cos(QUAD.t(:)/2);
    C2=cos(QUAD.t(:)/2+pi);
    S=sin(QUAD.t(:)/2);
    S2=sin(QUAD.t(:)/2+pi);
    
    % parameterization of boundary
    x=[ C; C2 ];
    y=[ S; -S ];
    
    % first derivative
    dx=0.5*[ -S; -S2 ];
    dy=0.5*[ C; -C ];
     
    % second derivative
    ddx=0.25*[ -C; -C2 ];
    ddy=0.25*[ -S; S ];
    
    % construct element object
    ELEMENT = class_element;
    ELEMENT=ELEMENT.init(num_edges,x,y,dx,dy,ddx,ddy,QUAD);
    
    % overwrite boundary length to exact value
    ELEMENT.arc_length=2*pi;
    
end

%%
function ELEMENT=PacMan(QUAD,mu)
    
    num_edges=3;
    N=2*QUAD.n;
    
    % precompute
    scale=1/(2*mu);
    
    T=QUAD.t(:);
    S=sin(scale*T);
    C=cos(scale*T);
    Z=zeros(N,1);
    E=ones(N,1);
    
    x0=cos(pi/mu);
    y0=sin(pi/mu);
    
    % parameterization of boundary
    x=[ C; (1-T/(2*pi))*x0; T/(2*pi) ];
    y=[ S; (1-T/(2*pi))*y0; Z ];
    
    % first derivative
    dx=[ -scale*S; -E*x0/(2*pi); E/(2*pi) ];
    dy=[ scale*C; -E*y0/(2*pi); Z];
    
    % second derivative
    ddx=[-scale*scale*C; Z; Z ];
    ddy=[ -scale*scale*S; Z; Z ];
     
    % construct element object
    ELEMENT = class_element;
    ELEMENT=ELEMENT.init(num_edges,x,y,dx,dy,ddx,ddy,QUAD);
    
    % overwrite boundary length to exact value
    ELEMENT.arc_length=pi/mu+2;

end

%% 
function ELEMENT=Shuriken(QUAD,amp)
    
    % define number of sampled points on boundary
    num_edges=4;
    N=2*QUAD.n;
    
    % precompute edges
    T=QUAD.t(:)/(2*pi);
    DT=ones(N,1)/(2*pi);
    Z=zeros(N,1);
    
    S=amp*sin(QUAD.t(:));
    C=amp*cos(QUAD.t(:));
    
    % parameterization of boundary
    x=[ T; 1-S; 1-T; S];
    y=[ S; T; 1-S; 1-T];
    
    % first derivative
    dx=[ DT; -C; -DT; C ];
    dy=[ C; DT; -C; -DT ];
     
    % second derivative
    ddx=[ Z; S; Z; -S ];
    ddy=[ -S; Z; S; Z ];
     
    % construct element object
    ELEMENT = class_element;
    ELEMENT=ELEMENT.init(num_edges,x,y,dx,dy,ddx,ddy,QUAD);
    
end

%% 
function ELEMENT=Jigsaw(QUAD,r,b)

    if b>=r
        error('Construction of Jigsaw Element requires b < r')
    end
    if r<=0
        error('Construction of Jigsaw Element requires r > 0')
    end
    if abs(b)+r>=0.5
        error('Construction of Jigsaw Element requires b+r < 0.5')
    end
    if sqrt(r*r-b*b)>=0.5
        error('Construction of Jigsaw Element requires sqrt(r^2-b^2) < 0.5')
    end

     % define number of sampled points on boundary
    num_edges=12;
    N=2*QUAD.n;
    
    tab_angle=pi+2*asin(b/r);
    alpha=tab_angle/(2*pi);
    theta=alpha*(QUAD.t(:)-pi);
    
    straight_length=0.5-sqrt(r*r-b*b);
    
    % precompute edges
    T=straight_length*QUAD.t(:)/(2*pi);
    DT=straight_length*ones(N,1)/(2*pi);
    Z=zeros(N,1);
    
    S=r*sin(theta);
    C=r*cos(theta);
    
    % first construct single three-edged 'side' with tab,
    % then rotate, reflect, and translate to construct remaining sides
    
    side_x = [ Z, C+b, Z];
    side_dx = [ Z, -alpha*S, Z ];
    side_ddx = [ Z, -alpha*alpha*C, Z ];
    
    side_y = [ T-0.5, S, T+0.5-straight_length];
    side_dy = [ DT, alpha*C, DT ];
    side_ddy = [ Z, -alpha*alpha*S, Z ];
    
    x = [ side_y, 0.5+side_x, -side_y, -0.5-side_x ];
    dx = [ side_dy, side_dx, -side_dy, -side_dx ];
    ddx = [ side_ddy, side_ddx, -side_ddy, -side_ddx ];
    
    y = [ -0.5+side_x, side_y, 0.5-side_x, -side_y ];
    dy = [ side_dx, side_dy, -side_dx, -side_dy ];
    ddy = [ side_ddx, side_ddy, -side_ddx, -side_ddy ];
    
    % construct element object
    ELEMENT = class_element;
    ELEMENT=ELEMENT.init(num_edges,x,y,dx,dy,ddx,ddy,QUAD);
    
    % overwrite boundary length to exact value
    ELEMENT.arc_length=8*straight_length+4*r*tab_angle;
    
end
