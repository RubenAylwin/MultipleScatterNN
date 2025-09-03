clear all
close all
clc
addpathGypsilab



N= 128;

Nper = 5;



Nsnapshotsself = 128;

k0=5;

pw = @(X,thetainc)  exp(i*k0*(cos(thetainc)*X(:,1)+sin(thetainc)*X(:,2)));

Gxy = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k0);

mesh0 = meshClosedBoundary(N,@(t) cos(t),@(t) sin(t));

gamma0 = dom(mesh0,3);
phi0 = fem(mesh0, 'P0');

invarg =@(X) atan2(X(:,2),X(:,1));

#initial configuration (define centers and radious)

Nscatters = 16;   %%
LBox  = 1; %sccaters contained in [-Lbox, Lbox]x [-Lbox, Lbox]

scaterxdim = ceil(sqrt(Nscatters));
grid = linspace(-LBox,LBox, scaterxdim+1);

hgrid = grid(2)-grid(1);

usablespacegrid = 0.9*hgrid;

rmin = 0.2*usablespacegrid/2;
rmax = 0.8*usablespacegrid/2;

dmin = 2*(hgrid-usablespacegrid);
dmax =  2*LBox;

cn = rmin*0.5*(1./(1:Nper).^2)';

rads  = rand(scaterxdim,scaterxdim)*(rmax-rmin)+rmin;

cmax = usablespacegrid/2-rads;

centers0 = (grid(2:end)+grid(1:end-1))*0.5;
centers0m = centers0+zeros(scaterxdim,scaterxdim);



centersX  = centers0m + (rand(scaterxdim,scaterxdim)-0.5).*cmax;
centersY  = centers0m' + (rand(scaterxdim,scaterxdim)-0.5).*cmax';



centers0X  = reshape(centers0m, size(centers0m,1)*size(centers0m,2),1);
centers0Y  = reshape(centers0m', size(centers0m,1)*size(centers0m,2),1);
centersX  = reshape(centersX, size(centersX,1)*size(centersX,2),1);
centersY  = reshape(centersY, size(centersY,1)*size(centersY,2),1);

rads  = reshape(rads, size(rads,1)*size(rads,2),1);

rXperN = cell(Nscatters,1);
rYperN = cell(Nscatters,1);

rXperNP = cell(Nscatters,1);
rYperNP = cell(Nscatters,1);

extraspaceX = usablespacegrid/2-(abs(centersX-centers0X)+rads);
extraspaceY = usablespacegrid/2-(abs(centersY-centers0Y)+rads);
extrarad = min(extraspaceX,extraspaceY);

for ii=1:Nscatters

    rXperN{ii} = @(t,ycenter,yrad,ys,yc) (centersX(ii)+ycenter*0.5*extraspaceX(ii))+...
        (rads(ii)+(yrad+0.5)*0.5*extrarad(ii))*cos(t)+...
        (sin(t*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
       (cos(t*(1:Nper))*(reshape(yc,Nper,1).*cn));

    rYperN{ii} = @(t,ycenter,yrad,ys,yc) (centersY(ii)+ycenter*0.5*extraspaceY(ii))+...
        (rads(ii)+(yrad+0.5)*0.5*extrarad(ii))*sin(t)+...
        (sin(t*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
       (cos(t*(1:Nper))*(reshape(yc,Nper,1).*cn));

##   rXperNP{ii} = @(t,ycenter,yrad,ys,yc) ...
##        -(rads(ii)+(yrad+0.5)*0.5*extrarad(ii))*sin(t)+...
##        ((cos(t*(1:Nper)).*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
##       ((-sin(t*(1:Nper)).*(1:Nper))*(reshape(yc,Nper,1).*cn));
##
##    rYperNP{ii} = @(t,ycenter,yrad,ys,yc) ...
##        (rads(ii)+(yrad+0.5)*0.5*extrarad(ii))*cos(t)+...
##        ((cos(t*(1:Nper)).*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
##       ((-sin(t*(1:Nper)).*(1:Nper))*(reshape(yc,Nper,1).*cn));

end

##ts = linspace(0,2*pi,64)';
##
##hold on;
##
##Pert = rand((3+4*Nper)*Nscatters,1)-0.5;
##
##for ii=1:Nscatters
##
##  Pertii = Pert( (1:(3+4*Nper))+(ii-1)*(3+4*Nper));
##
##  rx = rXperN{ii}(ts,Pertii(1),Pertii(2),Pertii(4:Nper+3),Pertii(Nper+4:2*Nper+3));
##  ry = rYperN{ii}(ts,Pertii(3),Pertii(2),Pertii(2*Nper+4:3*Nper+3),Pertii(3*Nper+4:4*Nper+3));
##
##  plot(rx,ry);
##
##end
##
##
##hold off;


Nsnapshots = 2048;
tolsvd = 1e-5;
toleim = 1e-4;

Us = zeros(N,Nsnapshots);
SelfMats = cell(Nsnapshots,1)


rX0 = @(t,yrad,ys,yc,yp)    (2*yp)*LBox+((rmax-rmin)*(yrad+0.5)+rmin)*cos(t)+...
    (sin(t*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
   (cos(t*(1:Nper))*(reshape(yc,Nper,1).*cn));

rY0 = @(t,yrad,ys,yc,yp)    (2*yp)*LBox+((rmax-rmin)*(yrad+0.5)+rmin)*sin(t)+...
    (sin(t*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
   (cos(t*(1:Nper))*(reshape(yc,Nper,1).*cn));

##rX0p = @(t,yrad,ys,yc)   -((rmax-rmin)*(yrad+0.5)+rmin)*sin(t)+...
##    ((cos(t*(1:Nper)).*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
##   ((-sin(t*(1:Nper)).*(1:Nper))*(reshape(yc,Nper,1).*cn));
##
##rY0p = @(t,yrad,ys,yc)   ((rmax-rmin)*(yrad+0.5)+rmin)*cos(t)+...
##    ((cos(t*(1:Nper)).*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
##   ((-sin(t*(1:Nper)).*(1:Nper))*(reshape(yc,Nper,1).*cn));

reg = 0.5/pi*regularize(gamma0,gamma0,phi0,'[log(r)]',phi0);

for nn=1:Nsnapshots

    Yh = halton_value(nn-1,4*Nper+4)-0.5;

    rx =@(t) rX0(t,Yh(4),Yh(5:Nper+4),Yh(Nper+5:2*Nper+4),Yh(2));
    ry =@(t) rY0(t,Yh(4),Yh(2*Nper+5:3*Nper+4),Yh(3*Nper+5:4*Nper+4),Yh(3));

    r =@(t) [rx(t),ry(t),zeros(size(t))];

##    rxp =@(t) rX0p(t,Yh(4),Yh(5:Nper+4),Yh(Nper+5:2*Nper+4),Yh(2));;
##    ryp =@(t) rY0p(t,Yh(4),Yh(2*Nper+5:3*Nper+4),Yh(3*Nper+5:4*Nper+4),Yh(3));
##    j =@(t) sqrt(rxp(t).^2+ryp(t).^2);

    ##discretization with included jacobian, it only affect if we want to plot the slution!
    Gxynn = @(X,Y) femGreenKernel(r(invarg(X)),r(invarg(Y)),'[H0(kr)]',k0);

    %%need to check the regularization under a pullback!
    v = -0.5/pi*integral(gamma0,gamma0,phi0,Gxynn,phi0,1e-6) -...
    reg;

    gnn =@(X) pw(r(invarg(X)),2*pi*(Yh(1)+0.5));

    gint = integral(gamma0,phi0,gnn);

    rhs = gint;

    Us(:,nn) = v\rhs;
    SelfMats{nn} = v;

    clear gnn;
    clear Gxynn;
    clear rx;
    clear ry;
    clear g;

end

Ndof = length(Us(:,1));

[U, S, V] = svd (Us,"econ");

singvals = diag(S)

s2sum = sum(singvals.^2);

sum2R = 0

RR = length(singvals);

for rr=1:length(singvals)

  sum2R = sum2R+singvals(rr)^2;


  if( sum2R/s2sum >= 1- tolsvd^2)
    RR = rr;
    break;
  end

end

UR = U(:,1:RR);






##################################RB TEST
##
##
##
##Ntest = 8;
##errors = zeros(Ntest,1);
##
##for nn=1:Ntest
##
##  Ytest = rand(1+Nscatters+Nscatters*4*Nper,1)-0.5;
##
##  mat = cell(Nscatters*(Nscatters+1)/2,1);
##  matc = cell(Nscatters,Nscatters);
##  rhs = zeros(Ndof*Nscatters,1);
##  rhsc = zeros(RR*Nscatters,1);
##
##  for ii=1:Nscatters
##
##    rxii =@(t) rXperN{ii}(t,Ytest(1+(4*Nper+1)*(ii-1)),Ytest(2+(4*Nper+1)*(ii-1)),...
##                          Ytest(3+(4*Nper+1)*(ii-1):(Nper+2)+(4*Nper+1)*(ii-1)),...
##                          Ytest(Nper+3+(4*Nper+1)*(ii-1):(2*Nper+2)+(4*Nper+1)*(ii-1)));
##
##
##
##    ryii =@(t) rYperN{ii}(t,Ytest(1+(4*Nper+1)*(ii-1)),Ytest(2+(4*Nper+1)*(ii-1)),...
##                          Ytest((2*Nper+3)+(4*Nper+1)*(ii-1):(3*Nper+2)+(4*Nper+1)*(ii-1)),...
##                          Ytest((3*Nper+3)+(4*Nper+1)*(ii-1):(4*Nper+2)+(4*Nper+1)*(ii-1)));
##
##
##    rii =@(t) [rxii(t),ryii(t),zeros(size(t))];
##
##    for jj=1:ii
##
##      rxjj =@(t) rXperN{jj}(t,Ytest(1+(4*Nper+1)*(jj-1)),Ytest(2+(4*Nper+1)*(jj-1)),...
##                            Ytest(3+(4*Nper+1)*(jj-1):(Nper+2)+(4*Nper+1)*(jj-1)),...
##                            Ytest(Nper+3+(4*Nper+1)*(jj-1):(2*Nper+2)+(4*Nper+1)*(jj-1)));
##
##      ryjj =@(t) rYperN{jj}(t,Ytest(1+(4*Nper+1)*(jj-1)),Ytest(2+(4*Nper+1)*(jj-1)),...
##                            Ytest((2*Nper+3)+(4*Nper+1)*(jj-1):(3*Nper+2)+(4*Nper+1)*(jj-1)),...
##                            Ytest((3*Nper+3)+(4*Nper+1)*(jj-1):(4*Nper+2)+(4*Nper+1)*(jj-1)));
##
##
##      rjj =@(t) [rxjj(t),ryjj(t),zeros(size(t))];
##
##
##
####        rxjj =@(t) rXperN{jj}(t,Ytest(1+(4*Nper+1)*(jj-1)),...
####                          Ytest(2+(4*Nper+1)*(jj-1):(Nper+1)+(4*Nper+1)*(jj-1)),...
####                          Ytest(Nper+2+(4*Nper+1)*(jj-1):(2*Nper+1)+(4*Nper+1)*(jj-1)));
####
####        ryjj =@(t) rYperN{jj}(t,Ytest(1+(4*Nper+1)*(jj-1)),...
####                              Ytest((2*Nper+2)+(4*Nper+1)*(jj-1):(3*Nper+1)+(4*Nper+1)*(jj-1)),...
####                              Ytest((3*Nper+2)+(4*Nper+1)*(jj-1):(4*Nper+1)+(4*Nper+1)*(jj-1)));
##
##
##        Gxynn = @(X,Y) femGreenKernel(rii(invarg(X)),rjj(invarg(Y)),'[H0(kr)]',k0);
##
##        ind = jj+ii*(ii-1)/2;
##
##        mat{ind} = -0.5/pi*integral(gamma0,gamma0,phi0,Gxynn,phi0,1e-6);
##
##        if (ii==jj)
##          mat{ind} = mat{ind}-reg;
##        end
##
##        matc{ii,jj} = UR'*mat{ind}*UR;
##        if(jj<ii)
##          matc{jj,ii} = UR'*mat{ind}.'*UR;
##        end
##
##        clear Gxynn;
##        clear rxjj;
##        clear ryjj;
##        clear rjj;
##
##
##
##    end
##
##    gii =@(X) pw(rii(invarg(X)),2*pi*(Ytest(end)+0.5));
##
##    rhs((1:Ndof)+(ii-1)*Ndof) = integral(gamma0,phi0,gii);
##
##    rhsc((1:RR)+(ii-1)*RR) = UR'* rhs((1:Ndof)+(ii-1)*Ndof);
##
##    clear gii;
##    clear rxii;
##    clear ryii;
##    clear rii;
##
##
##
##
##  end
##
##
##   prodfull = @(v) blockproductsym(mat,v,Ndof,Nscatters);
##   prodcomp = @(v) blockproduct(matc,v,RR,Nscatters);
##
##   ufull = gmres(prodfull, rhs, Ndof*Nscatters+1, 1e-10);
##
##   ucomp = gmres(prodcomp,rhsc,RR*Nscatters+1,1e-10);
##
##   ucomp2full = zeros(size(ufull));
##   ufull2comp = zeros(size(ucomp));
##   ufull2full = zeros(size(ufull));
##
##   for ii=1:Nscatters
##
##     ucomp2full((1:Ndof)+(ii-1)*Ndof) = UR*ucomp((1:RR)+(ii-1)*RR);
##
##     ufull2comp((1:RR)+(ii-1)*RR) = UR'*ufull((1:Ndof)+(ii-1)*Ndof);
##
##     ufull2full((1:Ndof)+(ii-1)*Ndof) = UR*ufull2comp((1:RR)+(ii-1)*RR);
##
##   end
##
##    errors(nn) = norm(ufull-ucomp2full)/norm(ufull);
##    errors2(nn) = norm(ufull-ufull2full)/norm(ufull);
##
##    clear prodfull;
##    clear prodcomp;
##
##end
##
##
##########################EIM Direct:

for nn=1:Nsnapshots

  A= SelfMats{nn}.full();

  Vmatss(:,nn) = reshape(A,size(A,1)*size(A,2),1);


end


##svd for interpolate matrix
[Uvs, Svs, Vvs] = svd (Vmatss(:,:),'econ');

singvalsvs = diag(Svs);

s2sum = sum(singvalsvs.^2);

sum2Q = 0

Q = length(singvalsvs);

for rr=1:length(singvalsvs)

  sum2Q = sum2Q+singvalsvs(rr)^2;

  if( sum2Q/s2sum >= 1- toleim^2)
    Q = rr;
    break;
  end

end

UQ = Uvs(:,1:Q);
##
##
############################store self comppressed:
##
for nn=1:Nsnapshots

  SelfMatsc{nn} = UR'*SelfMats{nn}*UR;


end
##
##
##################################################CROSS MATIX
##
rX1 = @(t,ypos,yrad,ys,yc) (dmin+(dmax-dmin)*(ypos+0.5)) +...
     ((rmax-rmin)*(yrad+0.5)+rmin)*cos(t)+...
    (sin(t*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
   (cos(t*(1:Nper))*(reshape(yc,Nper,1).*cn));

rY1 = @(t,yrad,ys,yc)   ((rmax-rmin)*(yrad+0.5)+rmin)*sin(t)+...
    (sin(t*(1:Nper))*(reshape(ys,Nper,1).*cn))+...
   (cos(t*(1:Nper))*(reshape(yc,Nper,1).*cn));

CrossMats = cell(2*Nsnapshots,1);

for nn=1:2*Nsnapshots

  Yh = halton_value(nn-1,8*Nper+3)-0.5;

  rx1 =@(t) rX1(t,Yh(1),Yh(2),Yh(4:Nper+3),Yh(Nper+4:2*Nper+3));
  ry1 =@(t) rY1(t,Yh(2),Yh(2*Nper+4:3*Nper+3),Yh(3*Nper+4:4*Nper+3));

  r1 = @(t) [rx1(t),ry1(t),zeros(size(t))];

  rx0 =@(t) rX0(t,Yh(3),Yh(4*Nper+4:5*Nper+3),Yh(5*Nper+4:6*Nper+3),0);
  ry0 =@(t) rY0(t,Yh(3),Yh(6*Nper+4:7*Nper+3),Yh(7*Nper+4:end),0);

  r0 = @(t) [rx0(t),ry0(t),zeros(size(t))];

  Gxynn = @(X,Y) femGreenKernel(r1(invarg(X)),r0(invarg(Y)),'[H0(kr)]',k0);


##  meshb0 = meshClosedBoundary(N,rx0,ry0);
##  meshb1 = meshClosedBoundary(N,rx1,ry1);
##
##  gamma0 = dom(meshb0,3);
##  phi0 = fem(meshb0, 'P0');
##
##  gamma1 = dom(meshb1,3);
##  phi1 = fem(meshb1, 'P0');

##  CrossMats{nn} = -0.5/pi*integral(gamma0,gamma1,phi0,Gxy,phi1,1e-6) -...
##    0.5/pi*regularize(gamma0,gamma1,phi0,'[log(r)]',phi1);

  CrossMats{nn} = -0.5/pi*integral(gamma0,gamma0,phi0,Gxynn,phi0,1e-6);

  clear r1;
  clear rx1;
  clear ry1;
  clear r0;
  clear rx0;
  clear ry0;



end


for nn=1:2*Nsnapshots

  A= CrossMats{nn}.full();

  Vmatss(:,nn) = reshape(A,size(A,1)*size(A,2),1);


end


##svd for interpolate matrix
[Uvsc, Svsc, Vvsc] = svd (Vmatss(:,:),'econ');

singvalsvsc = diag(Svs);

s2sumc = sum(singvalsvsc.^2);

sum2Qc = 0

Qc = length(singvalsvsc);

for rr=1:length(singvalsvsc)

  sum2Qc = sum2Qc+singvalsvsc(rr)^2;

  if( sum2Qc/s2sumc >= 1- toleim^2)
    Qc = rr;
    break;
  end

end

UQc = Uvsc(:,1:Q);


##########################store self comppressed:

for nn=1:2*Nsnapshots

  CrossfMatsc{nn} = UR'*CrossMats{nn}*UR;


end


##old!
##[Uvsc, Svsc, Vvsc] = svd (Vcomps(:,:),'econ');
##
##singvalsvsc = diag(Svsc);
##
##s2sumc = sum(singvalsvsc.^2);
##
##sum2Qc = 0
##
##Qc = length(singvalsvsc);
##
##for rr=1:length(singvalsvsc)
##
##  sum2Qc = sum2Qc+singvalsvsc(rr)^2;
##
##  if( sum2Qc/s2sumc >= 1- toleim^2)
##    Qc = rr;
##    break;
##  end
##
##end
##
##UQc = Uvsc(:,1:Qc);
##
##rowsN = 64;
##colsN = 13;
##
##dcomp = zeros(rowsN*colsN,Nsnapshots);
##
##jj=1;
##
##for nn=1:Nsnapshots
##
####  A =  reshape(Vmatss(:,nn) ,Ndof,Ndof);
####
####  A = A(78:128,14:64);
##
##  if(size(Vs{nn}.chd{1,3}.dat{1},2)!= 13)
##    continue;
##  end
##
##  A = Vs{nn}.chd{1,3}.dat{1};
##
##  A = reshape(A,size(A,1)*size(A,2),1);
##
##  dcomp(:,jj) = A;
##
##  jj=jj+1;
##
####  dcomp(:,nn) = diag(A);
##
##
##end
##
##dcomp = dcomp(:,1:jj-1);
##
##[Uvsd, Svsd, Vvsd] = svd (dcomp(:,:),'econ');
##
##singvalsvsd = diag(Svsd);
##
##s2sumd = sum(singvalsvsd.^2);
##
##sum2Qd = 0
##
##Qd = length(singvalsvsd);
##
##for rr=1:length(singvalsvsd)
##
##  sum2Qd = sum2Qd+singvalsvsd(rr)^2;
##
##  if( sum2Qd/s2sumd >= 1- toleim^2)
##    Qd = rr;
##    break;
##  end
##
##end
##
##UQd = Uvsd(:,1:Qd);
##
##
####rs = linspace(0,0.5-1e-2,128);
####
####theta = linspace(0,2*pi,129);
######theta = theta(1:end-1);
####
####[rs,theta] = meshgrid(rs,theta);
####
####x = rs.*cos(theta);
####y = rs.*sin(theta);
####
####
####z = zeros(size(x,1)*size(x,2),1);
####
####U = -0.5/pi*integral([ reshape(x,size(x,1)*size(x,2),1), reshape(y,size(y,1)*size(y,2),1),z],gamma,Gxy,phi,1e-6)*(Id\gint)...
####     +0.5/pi*integral([reshape(x,size(x,1)*size(x,2),1), reshape(y,size(y,1)*size(y,2),1),z],gamma,Kxy,ntimes(phi),1e-6)*u;
####
####contourf(x,y,reshape(real(U),size(x,1),size(x,2)),56,'LineColor','none');



