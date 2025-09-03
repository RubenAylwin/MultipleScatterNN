classdef fem
%+========================================================================+
%|                                                                        |
%|              OPENFEM - LIBRARY FOR FINITE ELEMENT METHOD               |
%|           openFem is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : fem.m                                         |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Finite element class definition               |
%|  `---'  |                                                              |
%+========================================================================+

properties
    typ = [];         % FINITE ELEMENT TYPE (P0, P1, P2, RWG, NED)
    opr = [];         % OPERATOR APPLIED TO FINITE ELEMENT
    msh = [];         % MESH FOR FINITE ELEMENT SPACE
    dir = [];         % MESH FOR DIRICHLET CONDITION
    jct = [];         % MESHES AND LINEAR COEFFESCIENT FOR JUNCTIONS
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function fe = fem(mesh,str)
        fe.typ = str;
        fe.opr = '[psi]';
        fe.msh = mesh;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot(varargin)
        fe  = varargin{1};
        spc = 'ob';
        if (nargin == 2)
            spc = varargin{2};
        end
        X = fe.unk;
        plot3(X(:,1),X(:,2),X(:,3),spc)
    end

    function surf(fe,V)
        V = feval(fe,V,fe.msh);
        if iscell(V)
            V = sqrt( V{1}.^2 + V{2}.^2 + V{3}.^2 );
        end
        plot(fe.msh,V);
    end

    function surf0(fe,V)
        plot3(fe.msh.ctr,V);
    end

    function graph(fe,V)
        V = feval(fe,V,fe.msh);
        if iscell(V)
            V = sqrt( V{1}.^2 + V{2}.^2 + V{3}.^2 );
        end
        nrm        = feval(fem(fe.msh,'P0'),fe.msh.nrm,fe.msh);
        fe.msh.vtx = fe.msh.vtx + (V*ones(1,3)).*nrm;
        plot(fe.msh,V)
    end

    function PlotVals(fe,vals)

       trisurf(fe.msh.elt,fe.msh.vtx(:,1),fe.msh.vtx(:,2),fe.msh.vtx(:,3),vals,'facecolor','interp','linestyle','none');


    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LENGTH
    function s = length(fe)
        s = size(fe.unk,1);
    end

    % SIZE
    function s = size(varargin)
        s = size(varargin{1}.unk);
        if (nargin == 2)
            s = s(varargin{2});
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESTRICTION OF THE BASIS FUNCTION
    function fe = rest(fe,i)
        fe.opr = ['[psi]',num2str(i)];
    end

    % CURL OF THE BASIS FUNCTION
    function fe = curl(fe)
        fe.opr = 'curl[psi]';
    end

    % DIVERGENCE OF THE BASIS FUNCTION
    function fe = div(fe)
        fe.opr = 'div[psi]';
    end

    % GRADIENT OF THE BASIS FUNCTION
    function fe = grad(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'grad[psi]';
        else
            fe.opr = ['grad[psi]',num2str(varargin{2})];
        end
    end

    % NORMAL TIMES BASIS FUNCTION
    function fe = ntimes(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'n*[psi]';
        else
            fe.opr = ['n*[psi]',num2str(varargin{2})];
        end
    end

    function fe = ntimes23(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'n23*[psi]';
        else
            fe.opr = ['n23*[psi]',num2str(varargin{2})];
        end
    end

    % QUADRATURE TIMES BASIS FUNCTION
    function fe = qtimes(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'q*[psi]';
        else
            fe.opr = ['q*[psi]',num2str(varargin{2})];
        end
    end

    % QUADRATURE DOT NORMAL OF THE BASIS FUNCTION
    function fe = qdotn(fe)
        fe.opr = 'qdotn*[psi]';
    end

    % NORMAL WEDGE BASIS FUNCTION
    function fe = nx(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'nx[psi]';
        else
            fe.opr = ['nx[psi]',num2str(varargin{2})];
        end
    end

    % NORMAL CROSS GRADIENT OF THE BASIS FUNCTION
    function fe = nxgrad(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'nxgrad[psi]';
        else
            fe.opr = ['nxgrad[psi]',num2str(varargin{2})];
        end
    end

    % DIVERGENCE NORMAL CROSS OF THE BASIS FUNCTION
    function fe = divnx(fe)
        fe.opr = 'curl[psi]';
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%% UNKNOWNS AND DOF %%%%%%%%%%%%%%%%%%%%%%%%%
    % DIRICHLET
    function fe = dirichlet(fe,mesh)
       fe.dir = mesh;
    end

    % JUNCTION
    function fe = junction(varargin)
        fe     = varargin{1};
        fe.jct = varargin(2:end);
    end

    % DEGREES OF FREEDOM
    function [X,elt2dof] = dof(fe)
        [X,elt2dof] = femDof(fe);
    end

    % UNKNOWNS AND REDUCTION MATRIX
    function [X,P] = unk(fe)
        [X,P] = femUnk(fe);
    end

    % UNKNOWM TO QUADRATURE MATRIX -> Mqud2dof x Mdof2unk
    function M = uqm(fe,domain)
        M = femUnk2Qud(fe,domain);
    end

    % UNKNOWNS DATA TO VERTEX
    function I = feval(v,f,mesh)
        bool   = (size(mesh.elt,2) == 4);
        gss    = 4*bool + 3*~bool;
        domain = dom(mesh,gss);
        u      = fem(mesh,'P1');
        M      = integral(domain,u,u);
        Fv     = integral(domain,u,v);
        if iscell(Fv)
            I{1} = M \ (Fv{1} * f);
            I{2} = M \ (Fv{2} * f);
            I{3} = M \ (Fv{3} * f);
        else
            I = M \ (Fv * f);
        end
    end

    function f = point2dof(fe,X)

      [elts,cords] = fe.msh.point2elts(X);

      Npoints = size(X,1);

      f = zeros(Npoints,fe.length());

      for ii=1:Npoints

##        if(isnan(elts(ii)))
##
##          ii
##
##        end

        inds   =  fe.msh.elt(elts(ii),:);

        f(ii,inds(1)) = 1-cords(ii,1)-cords(ii,2);
        f(ii,inds(2)) = cords(ii,1);
        if(size(inds,2) ==3)
          f(ii,inds(3)) = cords(ii,2);
        end

      endfor


  end

  function f = point2dof0(fe,X)

      [elts,cords] = fe.msh.point2elts(X);

      Npoints = size(X,1);

      f = zeros(Npoints,fe.length());

      for ii=1:Npoints

##        if(isnan(elts(ii)))
##
##          ii
##
##        end

        ind   =  elts(ii);

        f(ii,ind) = 1;


      endfor


    end



    function [f,ns] = point2dofQper(fe,X)

      ns = IntCZero(fe,0.5*IntFZero(fe,X(:,1)/(pi)));

      x = X;
      x (:,1) = x(:,1) - 2 * ns * pi;

      f = point2dof(fe,x);

   end

   ##extension of a function defined in [-pi, -pi+a] U [pi-a,pi]
   function f = point2dofExt(fe,X,a)

     ind1 = (abs(X(:,1)) < pi-a+1e-6);

     x = X;

     x(ind1,1) = (pi-a+1e-6) *( x(ind1,1)>=0)+ (-pi+a-1e-6) *( x(ind1,1)<0);

     f = point2dof(fe,x);

     f = f.*(1.0-W(X(:,1),(pi-a)*0.8,pi-a));

 end

    function [f,ns] = point2dofExtQper(fe,X,a)

     ns = IntCZero(fe,0.5*IntFZero(fe,X(:,1)/(pi)));

     x2 = X;

     x2 (:,1) = X(:,1) - 2 * ns * pi;

     ind1 = (abs(x2(:,1)) < pi-a+1e-6);

     x = x2;

     x(ind1,1) = (pi-a+1e-6) *( x(ind1,1)>=0)+ (-pi+a-1e-6) *( x(ind1,1)<0);

     f = point2dof(fe,x);

     f = f.*(1.0-W(x2(:,1),(pi-a)*0.8,pi-a));

end



    function [f,ns] = point2dofExtQper0(fe,X,a)

     ns = IntCZero(fe,0.5*IntFZero(fe,X(:,1)/(pi)));

     x2 = X;

     x2 (:,1) = X(:,1) - 2 * ns * pi;

     ind1 = (abs(x2(:,1)) < pi-a+1e-6);

     x = x2;

     x(ind1,1) = (pi-a+1e-6) *( x(ind1,1)>=0)+ (-pi+a-1e-6) *( x(ind1,1)<0);

     f = point2dof0(fe,x);

     f = f.*(1.0-W(x2(:,1),(pi-a)*0.8,pi-a));

 end


     function [f,ns] = point2dofExtQper2(fe,X,a)

     ns = IntCZero(fe,0.5*IntFZero(fe,X(:,1)/(pi)));

     x2 = X;

     x2 (:,1) = X(:,1) - 2 * ns * pi;

     ind1 = ((x2(:,1)) > pi-a-1e-6);
     ind2 = ((x2(:,1)) < -pi+a+1e-6);

     x = x2;

     x(ind1,1) = (pi-a+1e-6);
     x(ind2,1) = -pi+a+1e-6;

     f = point2dof(fe,x);

     f = f.*(W(x2(:,1),(pi-a)*0.8,pi-a));

end


##    function f = point2dof(fe,X)
##
##      [elts,cords] = fe.msh.point2elts(X);
##
##      Npoints = size(X,1);
##
##      f = zeros(Npoints,fe.length());
##
##      for ii=1:Npoints
##
####        if(isnan(elts(ii)))
####
####          ii
####
####        end
##
##        inds   =  fe.msh.elt(elts(ii),:);
##
##        f(ii,inds(1)) = 1-cords(ii,1)-cords(ii,2);
##        f(ii,inds(2)) = cords(ii,1);
##        if(size(inds,2) ==3)
##          f(ii,inds(3)) = cords(ii,2);
##        end
##
##      endfor
##
##
##  end

  function f = point2dof2dP0(fe,X)

    [elts,cords] = fe.msh.point2elts2d(X);

    Npoints = size(X,1);

    f = zeros(Npoints,fe.length());

    for ii=1:Npoints

      f(ii,elts(ii)) = 1;


    end



  endfunction


 function y =IntCZero(fe,x)
   y =x;
   y(x>0) =floor(x(x>0));
   y(x<0) = ceil(x(x<0));
 end

 function y = IntFZero(fe,x)
   y =  x;
   y(x>0) = ceil(x(x>0));
   y(x<0) = floor(x(x<0));
 end

    %%%%%%%%%%%%%%%%%%%%%%%%%%% USEFULL MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESTRICTION MATRIX
    function M = restriction(u,mesh)
        v         = fem(mesh,u.typ);
        [~,I1,I2] = intersect(v.dof,u.dof,'rows');
        M         = sparse(I1,I2,1,size(v.dof,1),size(u.dof,1));
    end

    % ELIMINATION MATRIX
    function M = elimination(u,mesh)
        v         = fem(mesh,u.typ);
        [~,nodir] = setdiff(u.dof,v.dof,'rows');
        M         = sparse(nodir,1:length(nodir),1,size(u.dof,1),length(nodir));
    end
end
end
