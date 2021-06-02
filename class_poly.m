%{
    class_poly

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

    Defines a sparse polynomial object for polynomials of two variables.
    The terms are sorted according to the bijection defined by (5) in 
    "Quadrature for Implicitly-defined Finite Element Functions on
    Curvilinear Polygons" by J. Ovall, S. Reynolds.

    A polynomial object can be initialized with
        p=class_poly;
        p=p.init(idx,coef);
    where idx is a vector containing the indicies of the nonzero 
    coefficients and coef is a vector containing the corresponding real
    coefficients. The init() method automatically consolidates repeated
    terms and sorts the index list.

    The zero polynomial is defined by 
        p=p.zero_poly();
    in which case p.idx=0 and p.coef=0 are chosen to represent the
    polynomial. Note that p.nz is the length of the arrays, so p.nz=1 for
    the zero polynomial (despite p having no "nonzero terms"
    mathematically). This is done to avoid empty arrays.

    To display the polynomial in a readable format, one my use the command
        p.print();

    To retrieve a multi-index alpha corresponding to an integer n:
        alpha = p.multi_index(n);
        n = p.inv_multi_index(alpha); % inverse map

    Algebraic methods:
        Sum of two polynomials: new_p = sum_poly(p1,p2);
        Product of two polynomials: new_p = prod_poly(p1,p2);

    Calculus methods:
        Gradient: [px,py] = p.grad_poly();
        Laplacian: L = p.laplacian_poly();
        Anti-Laplacian: P = p.anti_laplacian_poly();

    EXAMPLE:
        >> idx=[3,5]; coef=[1,-2];
        >> p=class_poly;
        >> p=p.init(idx,coef)
        p = 
          class_poly with properties:

             deg: 2
              nz: 2
             idx: [3 5]
            coef: [1 -2]
        >> p.print()
         + 1(x^2) - 2(y^2)

%}

classdef class_poly
    
    properties
        deg % polynomial degree
        nz % number of nonzero coefficients (nz=1 for zero polynomial)
        idx % indices of nonzero coefficients
        coef % nonzero coefficients
    end
    
    methods
        
        %% initialization method
        function obj=init(obj,idx,coef)
            
            % check that arrays are same length
            if(numel(coef)~=numel(idx))
                error('idx and coef must have same number of elements');
            end
            
            % assign to element fields
            obj.idx=idx(:).';
            obj.coef=coef(:).';
            
            % sort terms in ascending order
            [obj.idx,d]=sort(obj.idx);
            obj.coef=obj.coef(d);
            
            % consolidate coefficients
            N=numel(idx);
            for k=1:N
                for j=k+1:N
                    if obj.idx(k)==obj.idx(j)
                        obj.coef(k)=obj.coef(k)+obj.coef(j);
                        obj.coef(j)=0;
                    end
                end
            end
            
            % delete zero coefficients
            zero_entries = (abs(obj.coef)<1e-14);
            obj.idx(zero_entries)=[];
            obj.coef(zero_entries)=[];
            
            % when all coefficients are zero, return zero polynomial
            if isempty(obj.idx)
                obj=zero_poly(obj);
            else
                % number of nonzero coefficients
                obj.nz=numel(obj.idx);
                % polynomial degree
                n=obj.idx(obj.nz);
                obj.deg=floor(0.5*(sqrt(8*n+1)-1));
            end

        end
        
        %% zero polynomial
        function obj=zero_poly(obj)
            obj.deg=0;
            obj.nz=1;
            obj.idx=0;
            obj.coef=0;
        end
        
        %% print polynomial in a form a human can read
        function print(obj)
            
            for i=1:obj.nz

                alpha=obj.multi_index(obj.idx(i));

                if obj.coef(i)<0
                    fprintf(' - ');
                else
                    fprintf(' + ');
                end

                fprintf('%d',abs(obj.coef(i)));

                if alpha(1)==1
                    fprintf('(x)');
                elseif alpha(1)>1
                    fprintf('(x^%d)',alpha(1));
                end

                if alpha(2)==1
                    fprintf('(y)');
                elseif alpha(2)>1
                    fprintf('(y^%d)',alpha(2));
                end

            end
            
            fprintf('\n');
            
        end

        %% sum of two polynomials
        function obj=sum_poly(obj1,obj2)

            obj=class_poly;
            obj=obj.init([obj1.idx,obj2.idx],[obj1.coef,obj2.coef]);

        end

        %% product of two polynomials
        function obj=prod_poly(obj1,obj2)
            
            obj=class_poly;

            t=obj1.deg+obj2.deg;
            sz=(t+1)*(t+2)/2;

            new_idx=0:sz-1;
            new_coef=zeros(sz,1);

            for i=1:obj1.nz
                alpha=obj1.multi_index(obj1.idx(i));
                for j=1:obj2.nz
                    beta=obj2.multi_index(obj2.idx(j));
                    n=obj.inv_multi_index(alpha+beta);
                    new_coef(n+1)=new_coef(n+1)+obj1.coef(i)*obj2.coef(j);
                end
            end

            obj=obj.init(new_idx,new_coef);

        end

        %% gradient of a polynomial
        function [obj1,obj2]=grad_poly(obj)

            obj1=class_poly; % x derivative
            obj2=class_poly; % y derivative

            % initialize number of nonzeros to one
            sz_x=1;
            sz_y=1;

            idx_x=zeros(obj.nz,1); % indices for px
            idx_y=zeros(obj.nz,1); % indices for py

            coef_x=zeros(obj.nz,1); % coefficients for px
            coef_y=zeros(obj.nz,1); % coefficients for py

            for k=1:obj.nz
                alpha=obj.multi_index(obj.idx(k));
                % x derivative
                if alpha(1)>0
                    beta=alpha;
                    beta(1)=beta(1)-1;
                    idx_x(sz_x)=obj.inv_multi_index(beta);
                    coef_x(sz_x)=obj.coef(k)*alpha(1); % power rule
                    sz_x=sz_x+1;
                end
                % y derivative
                if alpha(2)>0
                    beta=alpha;
                    beta(2)=beta(2)-1;
                    idx_y(sz_y)=obj.inv_multi_index(beta);
                    coef_y(sz_y)=obj.coef(k)*alpha(2); % power rule
                    sz_y=sz_y+1;
                end
            end

            obj1=obj1.init(idx_x,coef_x);
            obj2=obj2.init(idx_y,coef_y);

        end

        %% Laplacian of a polynomial
        function new_obj=laplacian_poly(obj)

            pxx=class_poly; % two x derivatives
            pyy=class_poly; % two y derivatives

            % initialize number of nonzeros to one
            sz_x=1;
            sz_y=1;

            idx_x=zeros(obj.nz,1); % indices for pxx
            idx_y=zeros(obj.nz,1); % indices for pyy

            coef_x=zeros(obj.nz,1); % coefficients for pxy
            coef_y=zeros(obj.nz,1); % coefficients for pyy

            for k=1:obj.nz
                alpha=obj.multi_index(obj.idx(k));
                % x derivative
                if alpha(1)>1
                    beta=alpha;
                    beta(1)=beta(1)-2;
                    idx_x(sz_x)=obj.inv_multi_index(beta);
                    coef_x(sz_x)=obj.coef(k)*alpha(1)*(alpha(1)-1); % power rule twice
                    sz_x=sz_x+1;
                end
                % y derivative
                if alpha(2)>1
                    beta=alpha;
                    beta(2)=beta(2)-2;
                    idx_y(sz_y)=obj.inv_multi_index(beta);
                    coef_y(sz_y)=obj.coef(k)*alpha(2)*(alpha(2)-1); % power rule twice
                    sz_y=sz_y+1;
                end
            end

            pxx=pxx.init(idx_x,coef_x);
            pyy=pyy.init(idx_y,coef_y);

            new_obj=sum_poly(pxx,pyy);

        end

        %% anti-Laplacian of a polynomial
        function new_obj=anti_laplacian_poly(obj)
            new_obj=class_poly;
            new_obj=new_obj.zero_poly;
            for n=1:obj.nz
                P_alpha = obj.anti_laplacian_mono(obj.idx(n));
                P_alpha.coef = P_alpha.coef*obj.coef(n);
                new_obj = sum_poly(new_obj,P_alpha);
            end
        end

    end
    
    methods (Static)
        
        %% generate multi-index from master index
        function alpha=multi_index(n)
            % order of multi-index with master index n
            t=floor(0.5*(sqrt(8*n+1)-1));
            N=t*(t+1)/2;
            alpha(1)=t-n+N;
            alpha(2)=n-N;
        end
        
        %% retrieve master index from multi-index
        function n=inv_multi_index(alpha)
            t=alpha(1)+alpha(2);
            n=alpha(2)+t*(t+1)/2;
        end

        %% contruct anti-Laplacian of a monomial x^alpha, with n |-> alpha
        function obj=anti_laplacian_mono(n)

            m=floor(0.5*(sqrt(8*n+1)-1));
            if mod(m,2)==0
                N=m/2;
            else
                N=(m-1)/2;
            end

            % polynomial representation of |x-z|^2
            new_idx=[3,5];
            new_coef=[1,1];

            p1=class_poly;
            p1=p1.init(new_idx,new_coef);

            % (|x-z|^2)^{k+1}
            pk=p1;

            % (\Delta)^k |x-z|^\alpha
            new_idx=n;
            new_coef=1;

            Lk=class_poly;
            Lk=Lk.init(new_idx,new_coef);

            % first term, k=0
            obj=prod_poly(pk,Lk);

            % rescale coefficients
            denom=0.25/(m+1);
            obj.coef=obj.coef*denom;

            % sum over 1<= k <= N
            for k=1:N

                pk=prod_poly(pk,p1);
                Lk=laplacian_poly(Lk);

                temp=prod_poly(pk,Lk);
                denom=denom*(-0.25)/((k+1)*(m+1-k));
                temp.coef=temp.coef*denom;

                obj=sum_poly(obj,temp);

            end

        end
        
    end
end