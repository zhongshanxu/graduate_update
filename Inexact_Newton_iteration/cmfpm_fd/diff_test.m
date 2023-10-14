function diff_test(n,N)

% Computes the first and second derivatives of a function
% using an n-point centered formula and a grid of N nodes
%
% MIT 1.723 Computational methods for flow in porous media - Luis Cueto-Felgueroso April 2008
%


syms x
u= 3*(5-4*cos(2*x)^2)^-1;
du = diff(u,1);
ddu= diff(u,2);
Fu  = inline(vectorize(simplify(u)));
Fdu = inline(vectorize(simplify(du)));
Fddu= inline(vectorize(simplify(ddu)));

%Grid and differentiation matrices
h= 2*pi/N;
x= (0:h:2*pi)';
[D1,D2]=diffmatrix(x,n);

%Finite difference derivatives
u= Fu(x);
duFD = D1*u;
dduFD= D2*u;

%Exact derivatives
duex = Fdu(x);
dduex= Fddu(x);

%Errors and plots
error1= (1/max(abs(duex )))*sqrt(sum( (duFD-duex).^2 )/(N+1));
error2= (1/max(abs(dduex)))*sqrt(sum( (dduFD-dduex).^2 )/(N+1));

figure(1);plot(x,duFD,'k','marker','o','markersize',4);hold on;
figure(2);plot(x,dduFD,'k','marker','o','markersize',4);hold on;

xex= 0:2*pi/200:2*pi;
figure(1);plot(xex,Fdu(xex),'r');
title(['error=  ' num2str(error1)],'fontsize',14);axis tight;axis square
xlabel('X','fontsize',14);ylabel(['du/dx'],'fontsize',14);

figure(2);plot(xex,Fddu(xex),'r');
title(['error=  ' num2str(error2)],'fontsize',14);axis tight;axis square
xlabel('X','fontsize',14);ylabel(['d^2u/dx^2'],'fontsize',14);