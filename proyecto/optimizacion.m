%{
=================== TRABAJO DE OPTIMIZACIÓN ESTRUCTURAL ===================
-------------------| Universidad Nacional de Colombia |--------------------

Aplicación del método BESO al diseño de una estructura metálica

    Nelson Esteban Hernández Soto
    Heidy Paola Rodríguez Quevedo
===========================================================================
%}
clear
clc
nelx = 100;
nely = 40;
volfrac = 0.5;
er = 0.02;
rmin = 3;
nu = 0.3;
E = 1;
[U]=softbeso(nelx,nely,volfrac,er,rmin);

d11 = E/(1-nu^2);
d22 = d11;
d12 = E*nu/(1-nu^2);
d21 = d12;
d33 = E/(2*(1+nu));

D = [   d11 d12 0
        d21 d22 0
        0   0   d33];

B = [   -1/2    0       1/2     0       1/2     0       -1/2    0
        0       1/2     0       1/2     0       -1/2    0       -1/2
        1/2     -1/2    1/2     1/2     -1/2    1/2     -1/2    -1/2];

epsx = zeros(nely,nelx);
epsy = zeros(nely,nelx);
gammaxy = zeros(nely,nelx);
sx = zeros(nely,nelx);
sy = zeros(nely,nelx);
txy = zeros(nely,nelx);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    ae = U(edof);
    epsilon = B*ae;
    epsx(ely,elx) = epsilon(1);
    epsy(ely,elx) = epsilon(2);
    gammaxy(ely,elx) = epsilon(3);
    sigma = D*epsilon;
    sx(ely,elx) = sigma(1);
    sy(ely,elx) = sigma(2);
    txy(ely,elx) = sigma(3);
  end
end
epsz = -nu/E*(sx+sy);
tmax = 1/2*sqrt((sx-sy).^2+4*txy.^2); % esfuerzo cortante máximo
s1 = (sx+sy)/2+tmax;
s2 = (sx+sy)/2-tmax;
s3 = zeros(size(s1));
ang = atan2(2*txy, sx-sy)/2; % ángulo de inclinación de s1
sv   = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2); % von Mises
%% 

figure
colormap jet
subplot(311)
hold on;
colorbar;
title('\sigma_x')
imagesc(sx); axis equal tight;
subplot(312)
hold on;
colorbar;
title('\sigma_y')
imagesc(sy); axis equal tight;
subplot(313)
hold on;
colorbar;
title('\tau_{xy}')
imagesc(txy); axis equal tight;

figure
colormap jet
subplot(311)
hold on;
colorbar;
title('\sigma_{I}')
esc = 1;
imagesc(s1); axis equal tight;
quiver(cos(ang), sin(ang),... 
                esc, ...                  % con una escala esc
                'k',...                   % de color negro
                'ShowArrowHead','off') % una flecha sin cabeza
subplot(312)
hold on;
colorbar;
title('\sigma_{II}')
imagesc(s2); axis equal tight;
quiver(cos(ang+pi/2), sin(ang+pi/2),... 
                esc, ...                  % con una escala esc
                'k',...                   % de color negro
                'ShowArrowHead','off') % una flecha sin cabeza

subplot(313)
hold on;
colorbar;
title('\sigma_{v}')
imagesc(sv); axis equal tight;

