function [vtx,elt,col] = mshReadMsh(filename)
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
%|    #    |   FILE       : mshReadMsh.m                                  |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.12.2018                                    |
%| ( === ) |   SYNOPSIS   : Read .msh files (particle, edge, triangular   |
%|  `---'  |                and tetrahedral                               |
%+========================================================================+

% Open
fid = fopen(filename,'r');
if (fid==-1)
    error('mshReadMsh.m : cant open the file');
end

% Read file using keywords
str = fgets(fid);
while ~(str==-1)
    if strncmp(str,'$Nodes',6)
        vtx = readNodes(fid);
    elseif strncmp (str,'$Elements',9)
        [elt,col] = readElements(fid);
    end
    str = fgets(fid);
end

% Close file
fclose(fid);

% Only keep higher elements (tetra > triangle > edge > particles)
ord = sum(elt>0,2);
dim = max(ord);
elt = elt(ord==dim,1:dim);
col = col(ord==dim);
end


function vtx = readNodes(fid)
% Nodes
a = fgets(fid);
anum = str2num(a);
Nvtx = anum(end);
vtx  = zeros(Nvtx,3);

str = fgets(fid);

ii=1;
while  ~strncmp (str,'$EndNodes',9)


    strn = str2num(str);
    if(length(strn) > 1)

        toread = strn(end);
        for jj=1:toread
           str = fgets(fid);
        end

        for jj=1:toread
           tmp      = str2num(fgets(fid));
           vtx(ii,:) = tmp(1:3);
           ii=ii+1;
        end


    end


    str = fgets(fid);

end
end


function [elt,col] = readElements(fid)
% Initialize
a = str2num(fgets(fid));
Nelt = a(end);
elt  = zeros(Nelt,4);
col  = zeros(Nelt,1);

str = fgets(fid);
ii=1;
cc=1;

while  ~strncmp (str,'$EndElements',12)


    strn = str2num(str);
    if(strn(1)  !=2)
      str = fgets(fid);
      continue;
    end
    toread = strn(end);

    col(ii:(ii+toread))= cc;
    cc= cc+1;

    for jj=1:toread

          tmp      = str2num(fgets(fid));
##          while(true)
##            if( length(tmp) >2)
##             break;
##           end
##           tmp      = str2num(fgets(fid));
##         end
##         tmp
          elt(ii,1:3) = tmp(2:end);
          ii=ii+1;


    end


    str = fgets(fid);

end

end


