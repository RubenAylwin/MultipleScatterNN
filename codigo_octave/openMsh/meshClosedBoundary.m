function Z = meshClosedBoundary(N,rx,ry);

  h     = (2*pi)/N;
  theta = (0:h:(2*pi-h))';

  % Element
  vtx  =  [rx(theta) ry(theta) zeros(N,1)];
  elt  = [[(2:N)';1] (1:N)'];
  Z = msh(vtx,elt);


end
