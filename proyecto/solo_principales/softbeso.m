%%%%% A SOFT-KILL BESO CODE BY X. HUANG and Y.M. Xie  %%%%%
function [U, s1, s2, tmax, x] = softbeso(nelx,nely,volfrac,er,rmin)
%
% softbeso(100,40,0.5,0.02,3);
%
% INITIALIZE
x(1:nely,1:nelx) = 1.; vol=1.; i = 0; change = 1.; penal = 3.;
%fig = figure('Units','normalized','Position',[0 0 1 1]);
fig = figure;
fig.WindowState = 'maximized';
tiledlayout('flow', 'TileSpacing','compact', 'Padding', 'compact');
%fig.WindowState = 'maximized';
% START iTH ITERATION
while change > 0.001   
  i = i + 1;  vol = max(vol*(1-er),volfrac);
  if i >1; olddc = dc; end
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal);         
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS AND STRESSES
  [KE, D, B, E, nu] = lk;
  e_x = zeros(nely,nelx); e_y = zeros(nely,nelx); g_xy = zeros(nely,nelx);
  s_x = zeros(nely,nelx); s_y = zeros(nely,nelx); t_xy = zeros(nely,nelx);
  c(i) = 0.;
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2]);
      c(i) = c(i) + 0.5*x(ely,elx)^penal*Ue'*KE*Ue;
      dc(ely,elx) = 0.5*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
      epsilon = B*Ue;
      e_x(ely,elx) = epsilon(1);
      e_y(ely,elx) = epsilon(2);
      g_xy(ely,elx) = epsilon(3);
      sigma = D*epsilon;
      s_x(ely,elx) = sigma(1);
      s_y(ely,elx) = sigma(2);
      t_xy(ely,elx) = sigma(3);
    end
  end
  e_z = -nu/E*(s_x+s_y);
  tmax = 1/2*sqrt((s_x-s_y).^2+4*t_xy.^2); % esfuerzo cortante máximo
  s1 = (s_x+s_y)/2+tmax;
  s2 = (s_x+s_y)/2-tmax;
  s3 = zeros(size(s1));
  ang = atan2(2*t_xy, s_x-s_y)/2; % ángulo de inclinación de s1
  sv   = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2); % von Mises
% FILTERING OF SENSITIVITIES
[dc]   = check(nelx,nely,rmin,x,dc);
% STABLIZATION OF EVOLUTIONARY PROCESS
if i > 1; dc = (dc+olddc)/2.; end  
% BESO DESIGN UPDATE
[x]    = ADDDEL(nelx,nely,vol,dc,x);
% PRINT RESULTS
if i>10
 change=abs(sum(c(i-9:i-5))-sum(c(i-4:i)))/sum(c(i-4:i));
end
disp([' It.: ' sprintf('%4i',i) ' Obj.: ' sprintf('%10.4f',c(i)) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES
ax(1) = nexttile(1);
imagesc(-x); colormap(ax(1), gray);axis equal tight off;
title(sprintf('Iteración %d', i))
pause(1e-6);

% ax(2) = nexttile(2);
% imagesc(sv); axis equal tight; colormap(ax(2), jet);colorbar;
% title('\sigma_v');

% ax(3) = nexttile(3); 
% imagesc(e_z); axis equal tight; colormap(ax(3), jet); colorbar;
% title('$\varepsilon_{z}$ ','interpreter','latex');
% 
% ax(4) = nexttile(4); 
% imagesc(e_x); axis equal tight; colormap(ax(4), jet); colorbar;
% title('$\varepsilon_{x}$ ','interpreter','latex');
% 
% ax(5) = nexttile(5);
% imagesc(s_x); axis equal tight; colormap(ax(5), jet);colorbar;
% title('\sigma_x');

ax(6) = nexttile(2); hold off;
imagesc(s1); axis equal tight; colormap(ax(6), jet);colorbar;
hold on; 
quiver(s1.*cos(ang), s1.*sin(ang), 4, 'k','ShowArrowHead','off');
title('\sigma_{I}')

% ax(7) = nexttile(7); 
% imagesc(e_y); axis equal tight; colormap(ax(7), jet); colorbar;
% title('$\varepsilon_{y}$ ','interpreter','latex');
% 
% ax(8) = nexttile(8);
% imagesc(s_y); axis equal tight; colormap(ax(8), jet);colorbar;
% title('\sigma_y');

ax(9) = nexttile(3); hold off;
imagesc(s2); axis equal tight; colormap(ax(9), jet);colorbar;
hold on; 
quiver(s2.*cos(ang+pi/2), s2.*sin(ang+pi/2), 4, 'k','ShowArrowHead','off');
title('\sigma_{II}')

% ax(10) = nexttile(10); 
% imagesc(g_xy); axis equal tight; colormap(ax(10), jet); colorbar;
% title('\gamma_{xy}');
% 
% ax(11) = nexttile(11);
% imagesc(t_xy); axis equal tight; colormap(ax(11), jet);colorbar;
% title('\tau_{xy}');

ax(12) = nexttile(4); hold off;
imagesc(tmax); axis equal tight; colormap(ax(12), jet);colorbar;
hold on; 
quiver(tmax.*cos(ang+pi/4), tmax.*sin(ang+pi/4), 2, 'k','ShowArrowHead','off');
title('\tau_{max}')
end
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=ADDDEL(nelx,nely,volfra,dc,x)
l1 = min(min(dc)); l2 = max(max(dc));
while ((l2-l1)/l2 > 1.0e-5)
   th = (l1+l2)/2.0;
   x = max(0.001,sign(dc-th));
   if sum(sum(x))-volfra*(nelx*nely) > 0;
      l1 = th;
   else
      l2 = th;
   end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcf]=check(nelx,nely,rmin,x,dc)
dcf=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
       fac = rmin-sqrt((i-k)^2+(j-l)^2);
       sum = sum+max(0,fac);
       dcf(j,i) = dcf(j,i) + max(0,fac)*dc(l,k);
      end
    end
    dcf(j,i) = dcf(j,i)/sum;
  end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk; 
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (Cantilever)
F(2*(nelx+1)*(nely+1)-nely,1)=-10000; % 10 kN -> 1 tonf
fixeddofs=[1:2*(nely+1)];
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE, D, B, E, nu]=lk
E = 200000;  % MPa
nu = 0.3;
k=[ 1/2-nu/6  -1/8-nu/8 -1/4-nu/12  1/8-3*nu/8 ... 
   -1/4+nu/12  1/8+nu/8  nu/6      -1/8+3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
d11 = E/(1-nu^2);
d22 = d11;
d12 = E*nu/(1-nu^2);
d21 = d12;
d33 = E/(2*(1+nu));

D = [   d11 d12 0
        d21 d22 0
        0   0   d33];
% Esta matriz es de esta forma, para que 1 lado del elemento finito
% corresponda a 1 cm, teniendo unidades consistentes de mm, N y MPa
B = [   -1/20    0       1/20     0       1/20     0       -1/20    0
        0       1/20     0       1/20     0       -1/20    0       -1/20
        1/20     -1/20    1/20     1/20     -1/20    1/20     -1/20    -1/20];