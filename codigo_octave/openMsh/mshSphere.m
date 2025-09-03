function mesh = mshSphere(N,cx,cy,cz,r)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : mshSphere.m                                   |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Build uniform mesh for a sphere               |
%|  `---'  |                                                              |
%+========================================================================+

% Fibonnacci rules
##or     = (1+sqrt(5))/2;
##theta  = (mod((2*pi/or) .* (0:N-1),2*pi))';
##phi    = asin( -1 + 2/(N-1) * (0:N-1))';
##
##% Carthesian coordinates
##[x,y,z] = sph2cart(theta,phi,rad);
##X       = [0 0 0 ; x y z];
##
##% Delaunay triangulation
##DT        = delaunayTriangulation(X);
##[elt,vtx] = freeBoundary(DT);
##
##% Mesh
##mesh = msh(vtx,elt);

fid = fopen('sphere.geo','w');

fprintf(fid, 'SetFactory("OpenCASCADE"); \n');
fprintf(fid, 'Cx = DefineNumber[ %f, Name "Center X" ]; \n',cx);
fprintf(fid, 'Cy = DefineNumber[ %f, Name "Center Y" ]; \n',cy);
fprintf(fid, 'Cz = DefineNumber[ %f, Name "Center Z" ]; \n',cz);
fprintf(fid, 'r = DefineNumber[ %f, Name "radious." ]; \n',r);
fprintf(fid, 'Sphere(1) = {Cx, Cy, Cz, r, -Pi/2, Pi/2, 2*Pi}; \n');
fprintf(fid,'Mesh 2; \n');
fprintf(fid,'For k In{1:%d:1.0} \n',N(1));
fprintf(fid,'RefineMesh; \n');
fprintf(fid,'EndFor \n');
fprintf(fid,'Save "sphere.msh"; ');

fclose(fid);
 #linux...windows need some directory...
  system( "flatpak-spawn --host gmsh sphere.geo - -v 0");
## system( "gmsh geoaux.geo - -v 0");
 ##system( "gmsh-4.11.1-Windows64\\gmsh geoaux.geo - -v 0");

[vtx,elt,col] = mshReadMsh2('sphere.msh');

mesh =  msh(vtx,elt,col);
end
