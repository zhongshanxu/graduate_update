function poisson_test(n,N)

% Solves the 1D Poisson equation using finite differences
% n= n-point formula (n odd, n= 3,5,7,9,11,...)
% N+1= Number of grid points
%
% MIT 1.723 Computational methods for flow in porous media - Luis Cueto-Felgueroso April 2008
%

syms x
u= cos(3*x/2)-1;
ddu= diff(u,x,2);
Fu= inline(vectorize(simplify(u  )));
FS= inline(vectorize(simplify(ddu)));

%The system of equations for the 1D Poisson equation is L*u= S(x)
%The discrete laplacian is L= D2

%Grid and differentiation matrix
h= 2*pi/N;
x= (0:h:2*pi)';
[D1,D2]= diffmatrix(x,n);
L= D2;

%Source term
S= FS(x);

%Impose boundary conditions (Dirichlet at x= 0, Neumann at x= 2*pi)
%Dirichlet
L(1,:)= 0;L(1,1)= 1;
S(1)= 0;
%Neumann
L(end,:)= D1(end,:);
S(end)= 0;

%Finite difference solution
uFD= L\S;
%Exact solution
uex= Fu(x);

%Errors and plots
error= (1/max(abs(uex)))*sqrt(sum( (uFD-uex).^2 )/N);
figure(1);plot(x,uFD,'k','marker','o','markersize',4);hold on;
xex= 0:2*pi/200:2*pi;plot(xex,Fu(xex),'r');
title(['error=  ' num2str(error)],'fontsize',14);
xlabel('X','fontsize',14);ylabel(['u'],'fontsize',14);axis square;axis tight