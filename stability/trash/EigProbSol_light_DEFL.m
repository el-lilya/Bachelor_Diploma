%function [E0,res0,Ei,Xi,resEi,min_max,l,sing_1,sing_2,sing_3,sing_4]=EigProbSol_light_DEFL(Uss,coef,vtau,ne,nss,ind,deflation)
function [E0,res0,Ei,Xi,resEi,min_max,l,sing_1,sing_2,sing_3,rank]=EigProbSol_light_DEFL(Uss,coef,vtau,ne,nss,ind,deflation)

 [E0,~,res0]=EigenEstimator(Uss,coef,vtau,1:ne);
 %res0
 
%   E0 = [-0.000026 ;
%    -0.010143 ;
%    -0.011944 ;
%    -0.059947 ;
%    -0.080579 ;
%    -1.545909 ;
%    -3.076608 ];
E0(1:7)
if ind==1
    figure(nss+deflation*10)
    plot(real(E0),imag(E0),'bo')
    title('$\Lambda$','Interpreter','latex'); 
    hold on
end

[Ei,Xi,resEi,min_max,l,sing_1,sing_2,sing_3,sing_4,rank]=EigenSolver(Uss,coef,vtau,E0(1:7),deflation);
resEi
if ind==1
    figure(nss+deflation*10)
    plot(real(Ei),imag(Ei),'b+')
    title('$\Lambda$','Interpreter','latex'); 
    hold off
end

    
function [Lambda,X,resLambda]=EigenEstimator(Uss,coef,vtau,ne)

delta=0.01;  
%delta=0.0025;
scheme=2;

[L,~,M]=Matrices(Uss,coef,delta,vtau,scheme);

[X,Mu]=eig(full(M));
Mu=diag(Mu);
X=X(:,Mu~=0);
Mu=Mu(Mu~=0);
Lambda=log(Mu)/delta;
[~,se]=sort(real(Lambda),'descend');
Lambda=Lambda(se);

Mu=Mu(se);
X=X(:,se);

X=X(:,ne);
Lambda=Lambda(ne);
X = X(1:7,:); 
resLambda=zeros(1,length(Lambda));
for k=1:length(Lambda)
   X(:,k)=X(:,k)/Mu(k)^(vtau(end)-1);
   resLambda(k)=norm(formAz(L,vtau,Lambda(k),coef)*X(:,k),2)/...
       norm(formAz(L,vtau,Lambda(k),coef));
end
%nreal=size(Lambda(imag(Lambda)==0),1);
%res0=max(resLambda(1:7));


function [Lambda,Xi,resLambda,min_max,l,sing_1,sing_2,sing_3,sing_4,rank_A]=EigenSolver(Uss,coef,vtau,z0s,deflation)

lmax=150;
[L,~,~]=Matrices(Uss,coef,0.01,vtau,1);
nz0=size(z0s,1);
k=1;
X=[];
S=[];
resLambda = zeros(1,nz0);
min_max = zeros(1,nz0);
l = zeros(1,nz0);
while(k<=nz0)
    [x,nu,z,l(k)]=SLPM(z0s(k),X,S,L,vtau,coef,lmax,deflation);
    X=[X,x];
    if deflation~=0
        S=[S,nu;zeros(1,size(S,1)),z];
    else
        S=blkdiag(S,z);
    end
    k=k+1;
%end
[XS,ES]=eig(S);
Lambda=diag(ES);
[~,se]=sort(real(Lambda),'descend');
Lambda=Lambda(se);
XS=XS(:,se);
Xi=X*XS; 
%for i=1:nz0
for i=1:size(Xi,2)
    Xi(:,i)=Xi(:,i)/norm(Xi(:,i));
end
[~,Z,~]=svd(Xi);
Z=diag(Z);
min_max(k-1)=Z(end)/Z(1);
end

sing_1=zeros(1,nz0);
sing_2=zeros(1,nz0);
sing_3=zeros(1,nz0);
sing_4=zeros(1,nz0);
rank_A=zeros(1,nz0);
for i=1:nz0
    [~,Z,~]=svd(formAz(L,vtau,Lambda(i),coef));
    Z = diag(Z);
    sing_1(i)=Z(end)/Z(1);
    sing_2(i)=Z(end-1)/Z(1);
    sing_3(i)=Z(end-2)/Z(1);
    sing_4(i)=Z(end-3)/Z(1);
    rank_A(i) = rank(formAz(L,vtau,Lambda(i),coef));
end


for i=1:nz0
    resLambda(i)=norm(formAz(L,vtau,Lambda(i),coef)*Xi(:,i))/norm(formAz(L,vtau,Lambda(i),coef));
    %resLambda(i)=norm(formAz(L,vtau,Lambda(i),coef)*Xi(:,i));
end
%res=max(resLambda);


function [x,nu,z,l] = SLPM(z0,X,S,L,vtau,coef,lmax,deflation)
l=0;
z = z0;
figure(9);
hold off
plot(real(z),imag(z),'b+');
hold on
n = size(L{1,1},1);
res_need=1e-13;
res=Inf;
while(l<=lmax && res>res_need)   
%while res>res_need
    A = DEFL_Tau(z,X,S,L,vtau,coef,deflation);
    B = DEFL_dTaudz(z,X,S,L,vtau,coef,deflation);
    [u,e] = eig(A,B);
    e = diag(e);
    [~,se]=sort(abs(e),'ascend');
    e = e(se(1));
    z = z-e;      
    u = u(:,se(1));
    u = u/norm(u);
    res = norm(DEFL_Tau(z,X,S,L,vtau,coef,deflation)*u)/norm(DEFL_Tau(z,X,S,L,vtau,coef,deflation));
    l = l+1;
    figure(9);
    plot(real(z),imag(z),'r+')
    title('$\Lambda$','Interpreter','latex'); 
    hold on
end
if l==lmax
    res
end
x = u(1:n);
nu = u(n+1:end);

function Tau=DEFL_Tau(z,X,S,L,vtau,coef,deflation)
m=size(S,1);
A=[];
B=[];
Tau=formAz(L,vtau,z,coef);
if(m>0&&deflation~=0)
  T=Tau;   
  C=T*X/(z*eye(m)-S);
  if(deflation==1)
     A=X';
     B=zeros(m);
  elseif(deflation==2)  
     A=X'+z*(S')*(X');
     B=(S')*(X')*X;
  elseif(deflation==3)   
     A=X'+z*(S')*(X')+z^2*(S')^2*(X');
     B=(X*S)'*X+(X*S*S)'*X*S+z*(X*S*S)'*X;
  end 
  Tau=[T,C;A,B];
end


function dTau=DEFL_dTaudz(z,X,S,L,vtau,coef,deflation)
m=size(S,1);
n=size(L{1,1},1);
dA=[];
dB=[];
dTau=formdAdz(L,vtau,z,coef);
if(m>0&&deflation~=0)
  dT=dTau; 
  T=formAz(L,vtau,z,coef);
  dC=-T*X/(z*eye(m)-S)^2+dT*X/(z*eye(m)-S);
  if(deflation==1)
     dA=zeros(m,n);
     dB=zeros(m);
  elseif(deflation==2)  
     dA=(X*S)';
     dB=zeros(m);
  elseif(deflation==3)   
     dA=(X*S)'+2*z*(X*S*S)';
     dB=(X*S*S)'*X;
  end 
  dTau=[dT,dC;dA,dB];
end

function A=formAz(L,vtau,z,coef)

mu_C = coef(6);
w_C = vtau(2);
n = size(L{1,1},1); 
A=eye(n)*z-L{1,1};
for i=2:(size(L,2)-1)
    A=A-L{1,i}*exp(-vtau(i-1)*z);
end
eps = 1;
if (abs(z+mu_C)>=eps)
    A=A-(1-exp(-(z+mu_C)*w_C))/(mu_C+z)*L{1,end};
else
    %A = A-vtau(2)*L{1,end};
    A = A-Teyl_A(z,coef,vtau)*L{1,end};
end

function sum=Teyl_A(z,coef,vtau)
mu_C=coef(6);
w_C=vtau(2);
k=1;
next=w_C;
sum=0;
while (next~=0)
    sum=sum+next;
    k=k+1;
    next=-next*(z+mu_C)*w_C/k;
end


function B=formdAdz(L,vtau,z,coef)

mu_C = coef(6);
n=size(L{1,1},1);
B=eye(n);
for i=2:(size(L,2)-1)
    B=B+L{1,i}*(vtau(i-1)*exp(-vtau(i-1)*z));
end
eps=1;
if (abs(z+mu_C)>=eps)
    B=B-(vtau(2)*exp(-(z+mu_C)*vtau(2))/(mu_C+z)+(1-exp(-(z+mu_C)*vtau(2)))/...
    (mu_C+z)^2)*L{1,end};
else
    %B=B-vtau(2)^2/2*L{1,end};
    B=B-Teyl_B(z,coef,vtau);
end


function sum=Teyl_B(z,coef,vtau)
mu_C=coef(6);
w_C=vtau(2);
k=2;
s_k=w_C^2/2;
sum=0;
while (s_k*(1-k)~=0)
    sum=sum+s_k*(1-k);
    k=k+1;
    s_k=-s_k*(z+mu_C)*w_C/k;
end
