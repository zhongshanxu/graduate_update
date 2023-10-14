function D = diffper(m,n,N,h)

% Computes the m-th derivative of a periodic function
% using an n-point centered formula and a grid of N nodes
%
% MIT 1.723 Computational methods for flow in porous media - Luis Cueto-Felgueroso April 2008
%

    function [coefs]= fdcoefs(m,n,x,xi)
        
        %Finite difference weights (Fornberg's algorithm)
        % Computes the weights for a FD formula for the m-th derivative at x= xi,
        % using the nodes x= [x_0 x_1 ... x_n]
        %
        % m: Differentiation order
        % Number of points in the formula n+1: formal order n-m+1 (irregular points)
        %
        % MIT 1.723 Computational methods for flow in porous media - Luis Cueto-Felgueroso April 2008
        %
        c1= 1;
        c4= x(1)-xi;
        
        c= zeros(n+1,m+1);
        c(1,1)= 1;
        
        for i=1:n
            mn= min([i,m]);
            c2= 1;
            c5= c4;
            c4= x(i+1)-xi;
            for j= 0:i-1
                c3= x(i+1)-x(j+1);
                c2= c2*c3;
                for k= mn:-1:1
                    c(i+1,k+1)= c1*(k*c(i,k)-c5*c(i,k+1))/c2;
                end
                c(i+1,1)= -c1*c5*c(i,1)/c2;
                for k=mn:-1:1
                    c(j+1,k+1)= (c4*c(j+1,k+1)-k*c(j+1,k))/c3;
                end
                c(j+1,1)= c4*c(j+1,1)/c3;
            end
            c1= c2;
        end
        
        coefs= c(:,m+1)';
    end

%Coefficients of the difference scheme
x= -(n-1)/2:(n-1)/2;
[FDcoefs]= fdcoefs(m,n-1,x,0);
%Grid and differentiation matrix
FDcoefs= FDcoefs/(h^m);
R= [ FDcoefs((n-1)/2+1:n)   zeros(1,N-n)   FDcoefs(1:(n-1)/2) ];
C= [fliplr(FDcoefs(1:(n-1)/2+1))  zeros(1,N-n)  fliplr(FDcoefs((n-1)/2+2:n))];
D= toeplitz(sparse(C),sparse(R));
end