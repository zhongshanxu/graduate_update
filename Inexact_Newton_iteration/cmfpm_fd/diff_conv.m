function diff_conv(m,n,N)

syms x
u= 3*(5-4*cos(2*x)^2)^-1;
du = diff(u,m);
Fu = inline(vectorize(simplify(u )));
Fdu= inline(vectorize(simplify(du)));

if isnan(n)==0;
    x= -(n-1)/2:(n-1)/2;
    [FDcoefs0]= fdcoefs(m,n-1,x,0);
    
    for iN=1:length(N);
        h= 2*pi/N(iN);
        x= 0:h:2*pi;x= x(1:N(iN))';
        FDcoefs= FDcoefs0/(h^m);
        R= [ FDcoefs((n-1)/2+1:n)   zeros(1,N(iN)-n)   FDcoefs(1:(n-1)/2) ];
        C= [fliplr(FDcoefs(1:(n-1)/2+1))  zeros(1,N(iN)-n)  fliplr(FDcoefs((n-1)/2+2:n))];
        D= toeplitz(sparse(C),sparse(R));

        u= Fu(x);
        duFD= D*u;
        duex= Fdu(x);
        error(iN)= (1/max(abs(duex)))*sqrt(sum( (duFD-duex).^2 )/N(iN));
    end;
else;
    for iN=1:length(N);
        k= [0:(N(iN)/2-1) (-N(iN)/2):(-1)]';
        h= 2*pi/N(iN);
        x= 0:h:2*pi;x= x(1:N(iN))';
        
        u= Fu(x);
        i= sqrt(-1);
        duSP= real(ifft((i*k).^m.*fft(u)));
        duex= Fdu(x);
        error(iN)= (1/max(abs(duex)))*sqrt(sum( (duSP-duex).^2 )/N(iN));
    end;
end;


figure(1);loglog(N,error,'k','marker','o','markersize',4);axis square;hold on;
if isnan(n)==0;
    if mod(m,2)==0
        order= 5*error(1)*ones(1,length(error))./(N./N(1)).^(n-m+1);
    else;
        order= 5*error(1)*ones(1,length(error))./(N./N(1)).^(n-m);
    end;
    %loglog(N,order,'k--');xlabel('N');ylabel('error');
end;
xlabel('N','fontsize',14);ylabel('error','fontsize',14);
