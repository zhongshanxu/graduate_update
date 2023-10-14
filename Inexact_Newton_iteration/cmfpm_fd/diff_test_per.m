function diff_test_per(m,n,N)

% Computes the m-th derivative of a periodic function
% using an n-point centered formula and a grid of N nodes
%
% MIT 1.723 Computational methods for flow in porous media - Luis Cueto-Felgueroso April 2008
%

syms x
u= 3*(5-4*cos(2*x)^2)^-1;
du = diff(u,m);
Fu = inline(vectorize(simplify(u )));
Fdu= inline(vectorize(simplify(du)));

%Coefficients of the difference scheme
x= -(n-1)/2:(n-1)/2;
[FDcoefs]= fdcoefs(m,n-1,x,0);

%Grid and differentiation matrix
h= 2*pi/N;
x= 0:h:2*pi;x= x(1:N)';
FDcoefs= FDcoefs/(h^m);
R= [ FDcoefs((n-1)/2+1:n)   zeros(1,N-n)   FDcoefs(1:(n-1)/2) ];
C= [fliplr(FDcoefs(1:(n-1)/2+1))  zeros(1,N-n)  fliplr(FDcoefs((n-1)/2+2:n))];
D= toeplitz(sparse(C),sparse(R));

%Finite difference derivative
u= Fu(x);
duFD= D*u;

%Exact derivative
duex= Fdu(x);

%Errors and plots
error= (1/max(abs(duex)))*sqrt(sum( (duFD-duex).^2 )/N);
figure(2);plot(x,duFD,'k','marker','o','markersize',4);hold on;
xex= 0:2*pi/200:2*pi;plot(xex,Fdu(xex),'r');
title(['error=  ' num2str(error)],'fontsize',14);axis tight;axis square
xlabel('X','fontsize',14);ylabel(['d^' num2str(m) 'u/dx^' num2str(m)],'fontsize',14);