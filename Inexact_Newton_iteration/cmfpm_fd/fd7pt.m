function [D,D2]=fd7pt(x)

% fd7pt.m 
%
% Construct the 7-point centered difference approximation 
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
for k=4:N-2
   
    n=(k-3):(k+3);
    
    A=collocD(x(n));
    A2=A*A;
    
    D(k,n)=A(4,:);
    D2(k,n)=A2(4,:);
    
end

% Boundary points
A=collocD(x(1:7)); A2=A*A;
D(1:3,1:7)=A(1:3,:);

D2(1:3,1:7)=A2(1:3,:);

dex=(N-5):N1;
A=collocD(x(dex)); A2=A*A;
D(N+(-1:1),dex)=A(5:7,:);

D2(N+(-1:1),dex)=A2(5:7,:);


function D=collocD(x)
x=x(:); N=length(x); N1=N+1; N2=N*N;
X=repmat(x,1,N);                    Xdiff=X-X'+eye(N);
W=repmat(1./prod(Xdiff,2),1,N);    
D=W./(W'.*Xdiff); 
D(1:N1:N2)=1-sum(D);                D=-D';