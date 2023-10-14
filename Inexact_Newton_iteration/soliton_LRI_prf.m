%parameters
mu=5.6;Lz=1.0;Lx=24.0;
Nz=50;Nx=100;Nzx=Nz*Nx;N=2*Nzx;
[dz,z]=ched(Nz);dz=(2/Lz)*dz;z=(1+z)/2;
[dx,x]=ched(Nx);dx=(2/Lx)*dx;x=(Lx/2)*x;
Iz=speye(Nz);Ix=speye(Nx);I=speye(Nzx);A=zeros(N);
Dz=kron(dz,Ix);Dx=kron(Iz,dx);
Dzz=kron(dz^2,Ix);Dxx=kron(Iz,dx^2);
%有限差分微分矩阵
[fdz1, fdz2] = diffmatrix(z,3);
[fdx1, fdx2] = diffmatrix(x,3);
fDz = kron(fdz1, Ix); fDx = kron(Iz, fdx1);
fDzz = kron(fdz2, Ix);fDxx = kron(Iz,fdx2);
[zz,xx]=meshgrid(z,x);
Z=zz(:);X=xx(:);
z0=find(Z==0);
z1=find(Z==1);
xb=find(abs(X)==Lx/2);
bd=union(union(z0,z1),xb);
%initialization
psi=4*Z.*tanh(X);phi=mu*(1-Z);
seed=[psi;phi];sol=zeros(N,1);
fz=1-Z.^3;fpz=-3*Z.^2;
tol=1e-5;tau=4./7;
k=0;kmax=10;eps=1e-6;
%计算loop内非更新量，以存代算！（以空间换时间）
A11_out=fz.*Dxx+fz.^2.*Dzz+fz.*fpz.*Dz;
A22_out=Dxx+fz.*Dzz;
A11z1_out=fpz(z1).*Dz(z1,:)+Dxx(z1,:)-I(z1,:);
fA11_out=fz.*fDxx+fz.^2.*fDzz+fz.*fpz.*fDz;
fA22_0ut=fDxx+fz.*fDzz;
fA11z1_out=fpz(z1).*fDz(z1,:)+fDxx(z1,:)-I(z1,:);
tic
while k<kmax
    %computing b
    b1=phi.^2.*psi+fz.*(-Z.*psi+(Dxx*psi)+fpz.*(Dz*psi)+fz.*(Dzz*psi));
    b2=-2*phi.*psi.^2+(Dxx*phi)+fz.*(Dzz*phi);
    b1(z1)=fpz(z1).*(Dz(z1,:)*psi)+(Dxx(z1,:)*psi)-psi(z1);
    b1(z0)=psi(z0);
    b1(xb)=Dx(xb,:)*psi;
    b2(z1)=phi(z1);
    b2(z0)=phi(z0)-mu;
    b2(xb)=Dx(xb,:)*phi;
    b=-[b1;b2];
    %computing A
    A11=diag(-Z.*fz+phi.^2)+A11_out;
    A12=diag(2*phi.*psi);
    A21=diag(-4*psi.*phi);
    A22=diag(-2*psi.^2)+A22_out;
    A11(z0,:)=I(z0,:);A11(z1,:)=A11z1_out;
    A11(xb,:)=Dx(xb,:);
    A12(bd,:)=0;
    A21(bd,:)=0;
    A22(z0,:)=I(z0,:);A22(z1,:)=I(z1,:);
    A22(xb,:)=Dx(xb,:);
    %A=[A11,A12;A21,A22];
    A(1:Nzx,1:Nzx)=A11;A(1:Nzx,Nzx+1:N)=A12;
    A(Nzx+1:N,1:Nzx)=A21;A(Nzx+1:N,Nzx+1:N)=A22;
    %预处理因子
    fA11=sparse(1:Nzx,1:Nzx,-Z.*fz+phi.^2)+fA11_out;
    fA12=sparse(1:Nzx,1:Nzx,2*phi.*psi);
    fA21=sparse(1:Nzx,1:Nzx,-4*psi.*phi);
    fA22=sparse(1:Nzx,1:Nzx,-2*psi.^2)+fA22_0ut;
    fA11(z0,:)=I(z0,:);fA11(z1,:)=fA11z1_out;
    fA11(xb,:)=fDx(xb,:);
    fA12(bd,:)=0;
    fA21(bd,:)=0;
    fA22(z0,:)=I(z0,:);fA22(z1,:)=I(z1,:);
    fA22(xb,:)=fDx(xb,:);
    fA=[fA11,fA12;fA21,fA22];
    %Rechardson迭代
    r=(b-A*sol);
    while true
        pr = fA\r;Apr = A*pr;
        soln = sol + tau*pr;
        r = r - tau*Apr;
        tolrms = vecnorm(soln-sol)/sqrt(N)
        if tol < tolrms && tolrms <100
            sol = soln;
        else
            break;
        end
    end
    solrms = vecnorm(sol)/sqrt(N);
    seed=seed+sol;
    reseed=reshape(seed,[],2);
    psi=reseed(:,1);
    phi=reseed(:,2);
    if tolrms > 100
        iteration = 'Failed'
        break;
    elseif solrms < eps
        break;
    else
        k=k+1;
    end
end
toc
mesh(zz,xx,reshape(psi,[Nx,Nz]))
