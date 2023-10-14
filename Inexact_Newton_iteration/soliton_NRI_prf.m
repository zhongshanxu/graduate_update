%sparse原则：除非矩阵稀疏，否则不要稀疏化！强行稀疏化效率反而低！
%element-wise运算较慢
tic
%parameters
mu=5.6;Lz=1.0;Lx=24.0;
Nz=100;Nx=300;N=Nz*Nx;
[dz,z]=ched(Nz);dz=(2/Lz)*dz;z=(1+z)/2;
[dx,x]=ched(Nx);dx=(2/Lx)*dx;x=(Lx/2)*x;
Iz=speye(Nz);Ix=speye(Nx);I=speye(N);%O=sparse(N,N);
Dz=kron(dz,Ix);Dx=kron(Iz,dx);
Dzz=kron(dz^2,Ix);Dxx=kron(Iz,dx^2);
%有限差分微分矩阵
[fdz1, fdz2] = diffmatrix(z,3);
[fdx1, fdx2] = diffmatrix(x,3);
fDz = kron(fdz1,Ix); fDx = kron(Iz,fdx1);
fDzz = kron(fdz2,Ix);fDxx = kron(Iz,fdx2);
[zz,xx]=meshgrid(z,x);
Z=zz(:);X=xx(:);
z0=find(Z==0.);
z1=find(Z==1.);
xb=find(abs(X)==Lx/2);
bd=union(union(z0,z1),xb);
%初始化
Psi=4*Z.*tanh(X);Phi=mu*(1-Z);
seed=[Psi;Phi];scale=length(seed);
Fz=1-Z.^3;Fpz=-3*Z.^2;
k=0;kmax=20;eps=1e-6;tau=4./7;
%计算loop中非更新数据
fA11_out = Fz.*fDxx+Fz.^2.*fDzz+Fz.*Fpz.*fDz;
fA11z1_out=Fpz(z1).*fDz(z1,:)+fDxx(z1,:)-I(z1,:);
fA22_out = fDxx+Fz.*fDzz;
while k<kmax
    %谱方法离散化方程
    R1=Phi.^2.*Psi+Fz.*(-Z.*Psi+(Dxx*Psi)+Fpz.*(Dz*Psi)+Fz.*(Dzz*Psi));
    R2=-2*Phi.*Psi.^2+(Dxx*Phi)+Fz.*(Dzz*Phi);
    R1(z1)=Fpz(z1).*(Dz(z1,:)*Psi)+(Dxx(z1,:)*Psi)-Psi(z1);
    R1(z0)=Psi(z0);
    R1(xb)=Dx(xb,:)*Psi;
    R2(z1)=Phi(z1);
    R2(z0)=Phi(z0)-mu;
    R2(xb)=Dx(xb,:)*Phi;
    R=[R1;R2];
    %有限差分离散化雅可比
    fA11=sparse(1:N,1:N,-Z.*Fz+Phi.^2)+fA11_out;
    fA12=sparse(1:N,1:N,2*Phi.*Psi);
    fA21=sparse(1:N,1:N,-4*Psi.*Phi);
    fA22=sparse(1:N,1:N,-2*Psi.^2)+fA22_out;
    fA11(z0,:)=I(z0,:);fA11(z1,:)=fA11z1_out;
    fA11(xb,:)=fDx(xb,:);
    fA12(bd,:)=0.;
    fA21(bd,:)=0.;
    fA22(z0,:)=I(z0,:);fA22(z1,:)=I(z1,:);
    fA22(xb,:)=fDx(xb,:);
    fA=[fA11,fA12;fA21,fA22];
    sol = fA\(-tau*R);
    solrms = sqrt(sum(sol.^2)/scale)
    seed=seed+sol;
    reseed=reshape(seed,[],2);
    Psi=reseed(:,1);
    Phi=reseed(:,2);
    if solrms < eps
        break;
    else
        k=k+1;
    end
end
toc
mesh(zz,xx,reshape(Psi,[Nx,Nz]))