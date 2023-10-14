function [D,D2]=fd3pt(x)

% fd3pt.m 
%
% Construct the 3-point centered difference approximation 
% on the grid points specified by the vector x. Useful
% for forming preconditioners for pseuspectral methods.
%
%
% Written by: Greg von Winckel - 07/12/05
% Contact: gregvw(at)chtm(dot)unm(dot)edu 

N1=length(x); N=N1-1;

D=zeros(N1);
D2=D;

% Interior points
for k=2:N
   
    n=k-1:k+1;
   
    A=collocD(x(n));
    A2=A*A;
    
    D(k,n)=A(2,:);
    D2(k,n)=A2(2,:);
    
end

% Boundary points
A=collocD(x(1:3)); A2=A*A;
D(1,1:3)=A(1,:);
D2(1,1:3)=A2(1,:);

dex=(N-1):N1;
A=collocD(x(dex)); A2=A*A;
D(N1,dex)=A(3,:);
D2(N1,dex)=A2(3,:);


function D=collocD(x)
x=x(:); N=length(x); N1=N+1; N2=N*N;
X=repmat(x,1,N);                    Xdiff=X-X'+eye(N);
W=repmat(1./prod(Xdiff,2),1,N);    
D=W./(W'.*Xdiff); 
D(1:N1:N2)=1-sum(D);                D=-D';