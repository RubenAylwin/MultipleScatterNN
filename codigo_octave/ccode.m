clear all
close all
clc
addpathGypsilab

N = 180;
k = 0;
L = 1;

rx = [linspace(0,L,N/4), L*ones(1,N/4), linspace(L,0,N/4), zeros(1,N/4)];
ry = [zeros(1,N/4), linspace(0,L,N/4), L*ones(1,N/4), linspace(L,0,N/4)];
rz = zeros(1, N);  % z = 0

vertices = [rx(:), ry(:), rz(:)];
elements = [(1:N)', [2:N 1]'];
meshb = msh(vertices, elements);

gamma = dom(meshb,3);
phi = fem(meshb, 'P0');

g = @(X) X(:,2);

Gxy = @(X,Y) femGreenKernel(X,Y,'[log(r)]',k);
Kxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]1',k);
Kxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]2',k);
Kxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]3',k);

% Matrices BEM
V = -0.5/pi * integral(gamma,gamma,phi,Gxy,phi,1e-6) ...
    - 0.5/pi * regularize(gamma,gamma,phi,'[log(r)]',phi);

K = -0.5/pi * integral(gamma,gamma,phi,Kxy,ntimes(phi),1e-6) ...
    - 0.5/pi * regularize(gamma,gamma,phi,'grady[log(r)]',ntimes(phi));

Id = integral(gamma,phi,phi);
gint = integral(gamma,phi,g);

rhs = (gint/2 + K*(Id\gint));
u = V \ rhs;

[x, y] = meshgrid(linspace(0,L,128), linspace(0,L,128));
pts = [x(:), y(:), zeros(numel(x),1)];

U = -0.5/pi * integral(pts, gamma, Gxy, phi, 1e-6) * u ...
    + 0.5/pi * integral(pts, gamma, Kxy, ntimes(phi), 1e-6) * (Id\gint);

contourf(x, y, reshape(real(U), size(x)), 56, 'LineColor', 'none');
colorbar;
title('Solución de Laplace con condiciones de contorno dadas');
xlabel('x'); ylabel('y');
