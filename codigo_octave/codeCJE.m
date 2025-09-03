clc;
clear all;
close all;
addpathGypsilab

L = 1;
N = 50;

rx1 = linspace(0, L, N);
ry1 = zeros(1, N);
vertices1 = [rx1(:), ry1(:), zeros(N,1)];
elements1 = [(1:N-1)', (2:N)'];
meshD1 = msh(vertices1, elements1);

rx2 = linspace(L, 0, N);
ry2 = L * ones(1, N);
vertices2 = [rx2(:), ry2(:), zeros(N,1)];
elements2 = [(1:N-1)', (2:N)'];
meshD2 = msh(vertices2, elements2);

meshD = meshD1.union(meshD2);
meshD = meshD.swap();

rx3 = L * ones(1, N);
ry3 = linspace(0, L, N);
vertices3 = [rx3(:), ry3(:), zeros(N,1)];
elements3 = [(1:N-1)', (2:N)'];
meshN1 = msh(vertices3, elements3);

rx4 = zeros(1, N);
ry4 = linspace(L, 0, N);
vertices4 = [rx4(:), ry4(:), zeros(N,1)];
elements4 = [(1:N-1)', (2:N)'];
meshN2 = msh(vertices4, elements4);

meshN = meshN1.union(meshN2);
meshN = meshN.swap();

gammaD = dom(meshD, 3);
gammaN = dom(meshN, 3);
phiD = fem(meshD, 'P1');
phiN = fem(meshN, 'P1');

k = 0;
Gxy = @(X,Y) femGreenKernel(X,Y,'[log(r)]',k);
Kxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]1',k);
Kxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]2',k);
Kxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]3',k);

g = @(X) X(:,1) + X(:,2);

A = -1/(2*pi) * integral(gammaN, gammaN, phiN, Gxy, phiN, 1e-6) ...
    - 1/(2*pi) * regularize(gammaN, gammaN, phiN, '[log(r)]', phiN);  %VgammaN,gammaN

C = -1/(2*pi) * integral(gammaN, gammaD, phiN, Kxy, ntimes(phiD), 1e-6) ...
    - 1/(2*pi) * regularize(gammaN, gammaD, phiN, 'grady[log(r)]', ntimes(phiD));     %KgammaN,gammD

M = integral(gammaD, phiD, phiD);    %IgammaD
gint = integral(gammaD, phiD, g);  %g en gammaD

rhs = 0.5 * gint + C * (M \ gint);
lambda = A \ rhs;    %% V u  = (I/2 +Kgamma_N,gammaD) g

[x, y] = meshgrid(linspace(0, L, 100), linspace(0, L, 100));
pts = [x(:), y(:), zeros(numel(x), 1)];

U = -1/(2*pi) * integral(pts, gammaD, Gxy, phiD, 1e-6) * (M \ gint) ...
  + 1/(2*pi) * integral(pts, gammaN, Kxy, ntimes(phiN), 1e-6) * lambda;

figure;
contourf(x, y, reshape(real(U), size(x)), 56, 'LineColor', 'none');
colorbar;
title('Solución de Laplace con condiciones mixtas');
xlabel('x'); ylabel('y');
