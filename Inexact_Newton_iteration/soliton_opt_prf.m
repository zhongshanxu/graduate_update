%sparse原则：除非矩阵稀疏，否则不要稀疏化！强行稀疏化效率反而低！
%element-wise运算较慢
%parameters
tic
mu=5.6;Lz=1.0;Lx=24.0;
Nz=300;Nx=400;N=Nz*Nx;
%伪谱法微分矩阵
[dz,z]=ched(Nz);dz=(2/Lz)*dz;dz=dz';dzz=dz^2;z=(1+z)/2;
[dx,x]=ched(Nx);dx=(2/Lx)*dx;dxx=dx^2;x=(Lx/2)*x;
Iz=speye(Nz);Ix=speye(Nx);I=speye(N);%O=sparse(N,N);
%有限差分微分矩阵
[fdz1, fdz2] = diffmatrix(z,3);
[fdx1, fdx2] = diffmatrix(x,3);
fDz = kron(fdz1, Ix); fDx = kron(Iz, fdx1);
fDzz = kron(fdz2, Ix);fDxx = kron(Iz,fdx2);
[zz,xx]=meshgrid(z,x);
Z=zz(:);X=xx(:);
z0=find(Z==0.);
z1=find(Z==1.);
xb=find(abs(X)==Lx/2);
bd=union(union(z0,z1),xb);
%initialization
psi=4*zz.*tanh(xx);phi=mu*(1-zz);
Psi=psi(:);Phi=phi(:);
Seed=[Psi;Phi];scale=length(Seed);
fz=1-zz.^3;fpz=-3*zz.^2;
Fz=fz(:);Fpz=fpz(:);
k=0;kmax=20;eps=1e-6;tau=4./7;
%计算loop中非更新数据
fA11_out = Fz.*fDxx+Fz.^2.*fDzz+Fz.*Fpz.*fDz;
fA11z1_out=Fpz(z1).*fDz(z1,:)+fDxx(z1,:)-I(z1,:);
fA22_out = fDxx+Fz.*fDzz;
while k<kmax
    %谱方法离散化方程
    R1=phi.^2.*psi+fz.*(-zz.*psi+(dxx*psi)+fpz.*(psi*dz)+fz.*(psi*dzz));
    R2=-2*phi.*psi.^2+(dxx*phi)+fz.*(phi*dzz);
    R1(:,1)=fpz(:,1).*(psi*dz(:,1))+(dxx*psi(:,1))-psi(:,1);
    R1(:,end)=psi(:,end);
    R1(1,:)=dx(1,:)*psi;R1(end,:)=dx(end,:)*psi;
    R2(:,1)=phi(:,1);
    R2(:,end)=phi(:,end)-mu;
    R2(1,:)=dx(1,:)*phi;R2(end,:)=dx(end,:)*phi;
    R=[R1(:);R2(:)];
    %有限差分离散化雅可比
    fA11=spdiags(-Z.*Fz+Phi.^2,0,N,N)+fA11_out;
    fA12=spdiags(2*Phi.*Psi,0,N,N);
    fA21=spdiags(-4*Psi.*Phi,0,N,N);
    fA22=spdiags(-2*Psi.^2,0,N,N)+fA22_out;
    fA11(z0,:)=I(z0,:);fA11(z1,:)=fA11z1_out;
    fA11(xb,:)=fDx(xb,:);
    fA12(bd,:)=0;
    fA21(bd,:)=0;
    fA22(z0,:)=I(z0,:);fA22(z1,:)=I(z1,:);
    fA22(xb,:)=fDx(xb,:);
    fA=[fA11,fA12;fA21,fA22];
    sol = fA\(-tau*R);
    solrms = vecnorm(sol)/sqrt(scale)
    Seed=Seed+sol;
    reseed=reshape(Seed,[],2);
    Psi=reseed(:,1);
    Phi=reseed(:,2);
    psi=reshape(Psi,[Nx,Nz]);
    phi=reshape(Phi,[Nx,Nz]);
    if solrms < eps
        break;
    else
        k=k+1;
    end
end
toc
mesh(zz,xx,psi)
