function [E0,Ei,s1,s2,s3,s4]=EigProbSol_light(Uss,coef,vtau,ne,nss,ind)
E0=EigenEstimator(Uss,coef,vtau,1:ne)
%E0=[-0.0000;-0.0101;-0.0119; -0.0599; -0.0806; -1.5459;-3.0766];
%clf(nss)
if ind==1
    figure(nss)
    plot(real(E0),imag(E0),'bo')
    title('$\Lambda$','Interpreter','latex'); 
    hold on
end

[Ei,s1,s2,s3,s4]=EigenSolver(Uss,coef,vtau,E0(1:7));
if ind==1
    figure(nss)
    plot(real(Ei),imag(Ei),'b+')
    title('$\Lambda$','Interpreter','latex'); 
    hold off
end

    
function Lambda=EigenEstimator(Uss,coef,vtau,ne)

delta=0.01;  
scheme=2;

[~,~,M]=Matrices(Uss,coef,delta,vtau,scheme);

[~,Mu]=eig(full(M));
Mu=diag(Mu);
Mu=Mu(Mu~=0);
Lambda=log(Mu)/delta;
[~,se]=sort(real(Lambda),'descend');
Lambda=Lambda(se);
Lambda=Lambda(ne);



function [Lambda,s1,s2,s3,s4]=EigenSolver(Uss,coef,vtau,z0s)

lmax=10000;
[L,~,~]=Matrices(Uss,coef,0.01,vtau,1);
nz0=size(z0s,1);
k=1;
Lambda = zeros(nz0,1);
s1 = zeros(1,nz0);
s2 = zeros(1,nz0);
s3 = zeros(1,nz0);
s4 = zeros(1,nz0);
while(k<=nz0)
    [Lambda(k),s1(k),s2(k),s3(k),s4(k)]=SLPM_damp(z0s(k),L,vtau,coef,lmax);
    k=k+1;
end

function [z,s1,s2,s3,s4] = SLPM(z0,L,vtau,coef,lmax)

l=0;
sing_res_need = 1e-15;
z = z0;
sing_res=Inf;
figure(9);
plot(real(z),imag(z),'b+');
hold on
A = formAz(L,vtau,z,coef);
B = formdAdz(L,vtau,z,coef);
while(sing_res>sing_res_need && l<lmax)    
    [~,e] = eig(A,B);
    e = diag(e);
    [~,se]=sort(abs(e),'ascend');
    e = e(se(1));
    z = z-e;
    A = formAz(L,vtau,z,coef);
    B = formdAdz(L,vtau,z,coef);
    [~,Z,~]=svd(A);
    Z=diag(Z);
    sing_res=Z(end)/Z(1);
    %sing_res = norm(formAz(L,vtau,z,coef)*x)/norm(formAz(L,vtau,z,coef));
    l = l+1;
    figure(9);
    plot(real(z),imag(z),'r+')
    hold on   
end
hold off
s1=Z(end)/Z(1);
s2=Z(end-1)/Z(1);
s3=Z(end-2)/Z(1);
s4=Z(end-3)/Z(1);

function [z,s1,s2,s3,s4] = SLPM_damp(z0,L,vtau,coef,lmax)
% s_i=exp(-alpha*i)
alpha=0.5;
i=3;
s=exp(-alpha*i);
l=0;
res_need = 1e-15;
z = z0;
figure(9);
hold off
plot(real(z),imag(z),'b+');
hold on
A = formAz(L,vtau,z,coef);
B = formdAdz(L,vtau,z,coef);
[~,e] = eig(A,B);
e = diag(e);
[~,se]=sort(abs(e),'ascend');
e = e(se(1));
z=z-e*s;
A = formAz(L,vtau,z,coef);
B = formdAdz(L,vtau,z,coef);
[~,Z,~]=svd(A);
Z=diag(Z);
res0=Z(end)/Z(1)
while(res0>res_need && l<lmax)    
    [~,e] = eig(A,B);
    e = diag(e);
    [~,se]=sort(abs(e),'ascend');
    e = e(se(1));
    A = formAz(L,vtau,z-e*s,coef);
    [~,Z,~]=svd(A);
    Z=diag(Z);
    res=Z(end)/Z(1);
    while res>=res0
        i=i+1;
        s=exp(-alpha*i);
        A = formAz(L,vtau,z-e*s,coef);
        [~,Z,~]=svd(A);
        Z=diag(Z);
        res=Z(end)/Z(1);
        if z==z-e*s
            break 
        end
    end  
    if z==z-e*s
        break 
    end
    res0=res    
    z = z-e*s;
    B = formdAdz(L,vtau,z,coef);
    i=3;
    s=exp(-alpha*i);
    l = l+1;
%         figure(9);
%         z
%         plot(real(z),imag(z),'r+')
%         hold on   
    
end
hold off
s1=Z(end)/Z(1);
s2=Z(end-1)/Z(1);
s3=Z(end-2)/Z(1);
s4=Z(end-3)/Z(1);



function A=formAz(L,vtau,z,coef)

mu_C = coef(6);
n = size(L{1,1},1); 
A=eye(n)*z-L{1,1};
for i=2:(size(L,2)-1)
    A=A-L{1,i}*exp(-vtau(i-1)*z);
end
eps = 1;
if (abs(z+mu_C)>=eps)
    A=A-(1-exp(-(z+mu_C)*vtau(2)))/(mu_C+z)*L{1,end};
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
