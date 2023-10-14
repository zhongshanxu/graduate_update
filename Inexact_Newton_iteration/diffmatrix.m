function [D1,D2]=diffmatrix(x,nt)

% Computes the first and second order differentiation matrices 
% associated to a grid with nodes x= [x_0 x_1 ... x_N]
% Constructs formulas of nt points (nt odd, nt= 3,5,7,9,11,...)
% Generalization of a function by Greg von Winckel
%
% MIT 1.723 Computational methods for flow in porous media - Luis Cueto-Felgueroso April 2008
%


q= (nt-1)/2;

N1=length(x); N=N1-1;

D1= zeros(N1);
D2=D1;

% Interior points
for k=(q+1):N-(q-1)
   
    n=(k-q):(k+q);
    
    A=collocD(x(n));
    A2=A*A;
    
    D1(k,n)=A(q+1,:);
    D2(k,n)=A2(q+1,:);
end

% Boundary points
A=collocD(x(1:nt)); A2=A*A;
D1(1:q,1:nt)=A(1:q,:);

D2(1:q,1:nt)=A2(1:q,:);

dex=(N-(nt-2)):N1;
A=collocD(x(dex)); A2=A*A;
D1(N+(-(q-2):1),dex)=A(q+2:nt,:);

D2(N+(-(q-2):1),dex)=A2(q+2:nt,:);

D1= sparse(D1);
D2= sparse(D2);



function D=collocD(x)
x=x(:); N=length(x); N1=N+1; N2=N*N;
X=repmat(x,1,N);                    Xdiff=X-X'+eye(N);
W=repmat(1./prod(Xdiff,2),1,N);    
D=W./(W'.*Xdiff); 
D(1:N1:N2)=1-sum(D);                D=-D';