function [Ei,Xi,resEi]=EigProbSol(Uss,coef,vtau,aoev,Esolver,deflation,addinf2,nss,ind,ifig)
%
% Latest revision 28.04.2019 
%
% Вычисление собственных значений и векторов
%
% Входные данные: Uss - стационарное состояние
%                 coef - парамеры модели
%                 aoev - число собственных значений, которые необходимо найти
%                 Esolver - решатель частичных линейных проблем
%                 deflation - режим дефляции
%                 addinf2 - информация о потери ведущего собственного значения (0/1)
%                           + пошаговой информации о работе решателей частичных проблем (2)
%                 ind - нужно ли отрисовывать собственные значения на комплексной
%                       плоскости
%                 nss - номер стационарного состояния
%                 ifig - номер рисунка, 
%
% Выходные данные: Ei - собственные значения,
%                  Xi - собственные вектора,
%                  resEi - погрешность вычисления собственных значений
%
% Внешние процедуры: Matrices
%
% Внутренние процедуры: dAdz, DEFL, EigenEstimator, EigenSolver, formAz,
%                       formdAdz, SFMlight, SLPM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[Ei,~]=EigenEstimator(Uss,coef,vtau,[1:aoev]);
[Ei,Xi,resEi]=EigenSolver(Uss,coef,vtau,Ei,0,Esolver,deflation,addinf2);

%Отрисовка собственных значений
if ind==1
    figure(nss*10+ifig)
        plot(real(Ei),imag(Ei),'k+')
        title('$\Lambda$','Interpreter','latex'); 
        hold on
end

    
function [Lambda,X]=EigenEstimator(Uss,coef,vtau,ne)
%
% Latest revision 02.01.2018 
%
%
delta=1.0e-2; % шаг сетки
scheme=2; % схема интегрирования: 1 - неявный метод Эйлера,
          %                       2 - BDF2.

[~,~,M]=Matrices(Uss,coef,delta,vtau,scheme);
%
[X,Mu]=eig(full(M));
Mu=diag(Mu);
X=X(:,Mu~=0);
Mu=Mu(Mu~=0);
Lambda=log(Mu)/delta;
[~,se]=sort(real(Lambda),'descend');
Lambda=Lambda(se);
Mu=Mu(se);
X=X(:,se);
%
% Вычисление невязок для Mu:
%Mu=Mu(ne);
X=X(:,ne);
%for j=1:length(ne)
%   fprintf(1,'ne=%d, Mu=%g+%gi\n',ne(j),real(Mu(j)),imag(Mu(j)));
%   X(:,j)=X(:,j)/norm(X(:,j));
%   fprintf(1,'SpectrRes_Mu=%d\n',norm(full(M*X(:,j))-Mu(j)*X(:,j)));
%end
%
Lambda=Lambda(ne);
X=X(1:4,:);
% Вычисление невязок для Lambda:
%resLambda=zeros(size(Lambda));
for k=1:length(Lambda)
   X(:,k)=X(:,k)/norm(X(:,k));
%   resLambda(k)=norm(formAz(MATR,Lambda(k))*X(:,k),2);
end
%maxresLambda=max(resLambda);


function [Lambda,X,maxresLambda]=EigenSolver(Uss,coef,vtau,z0s,ieig,Esolver,...
                                                          deflation,addinf)
%
% Latest revision 22.01.2018  
%
R=10.0; 
%
% SFM parameters:
lmaxsf=25;   
fmin=1.0e-28; 
ld0=1;  %damper=qdamp*exp(-pdamp*ld)      
ldmax=8;     
pdamp=2.0;          
qdamp=2.0;
%
% SLP parameters:
lmaxlp=10;   
resmin=1.0e-12; 
%
[L,~,~]=Matrices(Uss,coef,1,vtau,1);
n=size(L{1,1},1);
X=[];
S=[];
%
nz0=size(z0s,1);
m=0;
k=1;
while(k<=nz0)
   switch Esolver
      case 1 
      [z,~,x,infrm]=SFMlight(L,vtau,X,S,z0s(k),R,...
                    deflation,lmaxsf,fmin,ld0,ldmax,pdamp,qdamp,addinf);
%  infrm=0   no eigenvalues are found for lmaxsf steps  
%  infrm=1   an eigenvalue is found
%  infrm=10  gradf(z)=0
%  infrm=20  ldamp>ldmax
%  infrm=30  computed z is out of region z: |z-z0|<=R
      case 2 
      [z,~,x,infrm]=SLPM(L,vtau,X,S,z0s(k),R,deflation,lmaxlp,resmin,addinf);
%  infrm=0   no eigenvalues are found for lmaxlp steps  
%  infrm=1   an eigenvalue is found
%  infrm=30  computed z is out of region z: |z-z0|<=R      
   end
   if(addinf>0)
      fprintf(1,'Причина останова infrm=%g\n',infrm);       
   end   
   if(infrm==1||ieig==0)
      nx=norm(x(1:n),2);
      X=[X,x(1:n)/nx];
      if(deflation~=0)
         S=[S,x(n+1:end)/nx;zeros(1,m),z];
      else  
         S=blkdiag(S,z);
      end
      m=m+1;
   end
   k=k+1;
end
[XS,ES]=eig(S);
Lambda=diag(ES);
[~,se]=sort(real(Lambda),'descend');
Lambda=Lambda(se);
XS=XS(:,se);
X=X*XS; 
%
% Вычисление невязок:
resLambda=zeros(size(Lambda));
for k=1:length(Lambda)
   X(:,k)=X(:,k)/norm(X(:,k));
   resLambda(k)=norm(formAz(L,vtau,Lambda(k))*X(:,k),2);
end
maxresLambda=max(resLambda);

function [z,f,x,infrm]=SFMlight(L,vtau,X,S,z0,R,... 
                        deflation,lmaxsf,fmin,ld0,ldmax,pdamp,qdamp,addinf)
%
% Latest revision 22.01.2018
%
ldamp=ld0+1;
l=1;
z=z0;
m=size(S,1);
n=size(L{1,1},1);
fold=0;
while (l<=lmaxsf+1)
% для текущего z вычисляем f и x 
   [A,P,Sz,XSz,XS,XS2]=DEFL(L,vtau,z,X,S,deflation);
   [~,~,VV]=svd(A);
   x=VV(:,end);
   Ax=A*x;
   f=norm(Ax,2)^2;
   am2=max(abs(A(:)))^2;
   if(addinf>0) 
      fprintf(1,'l=%g, ldamp=%g, z=%g, f=%g, am2=%g, m=%g\n',...
                                                        l,ldamp,z,f,am2,m);  
   end
% проверяем достижение заданной точности:    
   if(f<=am2*fmin)
      infrm=1; % заданная точность достигнута
      return
   end
% проверяем успешность шага:  
   flag=f<fold||l==1;
   if(flag)
      zold=z;
      fold=f;
      xold=x;
   end   
% проверяем не достигнуто ли максимальное колличесво шагов:        
   if(l>lmaxsf)
      z=zold;
      f=fold;
      x=xold;
      infrm=0; % достигнуто максимальное колличество шагов, но
      return   % заданная точность не достигнута
   end 
% формируем новый шаг SFM: 
   if(flag)
      v=dAdz(L,vtau,z,x(1:n)); % v=dA(z)/dz*x
      if(m>0&&deflation~=0)
         xx=x(n+1:end);
         if(deflation==-1)
            v=[v;zeros(m,1)]; 
         else  
            if(deflation==1)
               vv=zeros(m,1); 
            elseif(deflation==2)
               vv=XS'*x(1:n); 
            elseif(deflation==3) 
               vv=XS'*x(1:n)+2*z*(XS2'*x(1:n))+XS2'*(X*xx); 
            end 
            v=[v+dAdz(L,vtau,z,XSz*xx)-P*(XSz*(Sz\xx));vv];
         end   
      end   
      gradf=2.0*(v'*Ax);
      gfr=gradf; 
      agfr=abs(gfr);
% проверяем неравенство нулю градиента:         
      if(agfr==0.0) 
         infrm=10; % градиент равен нулю
         return    % заданная точность не достигнута
      else
         step=f/agfr;
         gfrdir=gfr/agfr;
         ldamp=max(ldamp-1,0);
      end 
   else
% будем шагать с меньшим демпфером:     
      ldamp=ldamp+1;
      if(ldamp>ldmax)
         z=zold;
         f=fold;
         x=xold; 
         infrm=20; % возможность уменьшения демпфера исчерпана
         return
      end   
   end
% шагаем:       
   damper=qdamp*exp(-pdamp*ldamp);
   z=zold-damper*step*gfrdir;
% проверяем не выскочили ли за пределы заданной области:   
   if(abs(z-z0)>R)
      z=zold;
      f=fold;
      x=xold; 
      infrm=30; % выскочили за пределы заданной области
      return
   end
   l=l+1;    
end
%
function [z,res,x,infrm]=SLPM(L,vtau,X,S,z0,R,deflation,lmaxlp,resmin,addinf)
%
% Latest revision 24.01.2018
%
l=1;
z=z0;
m=size(S,1);
n=size(L{1,1},1);
res=Inf;
x=[];
if(addinf>0) 
   fprintf(1,'l=%g, z=%g, m=%g\n',l,z,m);  
end   
while(l<=lmaxlp+1)
% формируем A(z):       
   [A,P,Sz,XSz,XS,XS2]=DEFL(L,vtau,z,X,S,deflation);
   if(l>1)
      res=norm(A*x,2); 
      am=max(abs(A(:)));
      if(addinf>0) 
         fprintf(1,'l=%g, z=%g, res=%g, am=%g, m=%g\n',l,z,res,am,m);  
      end   
% проверяем достижение заданной точности:    
      if(res<=am*resmin)
         infrm=1; % заданная точность достигнута
         return
      end
% проверяем не достигнуто ли максимальное колличесво шагов:        
      if(l>lmaxlp)
         infrm=0; % достигнуто максимальное колличество шагов, но
         return   % заданная точность не достигнута
      end       
   end 
% формируем B=dA/dz(z):       
   B=formdAdz(L,vtau,z);
   if(m>0&&deflation~=0)
      C=zeros(m,n);D=zeros(m,m);E=zeros(n,m);
      if(deflation>0)
         E=B*XSz-P*(XSz/Sz);
         if(deflation==2)
            C=XS';
         elseif(deflation==3) 
            C=XS'+2*z*XS2'; D=XS2'*X;
         end        
      end
      B=[B,E;C,D];
   end 
% решаем линейную проблему:   
   [u,e]=eig(A,-B); 
   e=diag(e);
   [~,se]=sort(abs(e),'ascend'); 
   e=e(se(1));
% делаем шаг и проверяем не выскочили ли за пределы заданной области:  
   z=z+e; 
   if(abs(z-z0)>R&&l>1)
      z=z-e;
      infrm=30; % выскочили за пределы заданной области
      return
   end
   x=u(:,se(1)); 
   l=l+1;
end

function [A,P,Sz,XSz,XS,XS2]=DEFL(L,vtau,z,X,S,deflation)
%
% Latest revision 22.01.2018
%
m=size(S,1);
P=[];
Sz=[];
XSz=[];
XS=[];
XS2=[];
A=formAz(L,vtau,z);
if(m>0&&deflation~=0)
   P=A;
   if(deflation==-1)
      A=[P,-X;X',zeros(m,m)]; 
   else 
      Sz=z*eye(m)-S;
      XSz=X/Sz;
      if(deflation==1)
         LBC=[X',zeros(m)];
      elseif(deflation==2)  
         XS=X*S;
         LBC=[X'+z*XS',XS'*X];
      elseif(deflation==3)  
         XS=X*S;
         XS2=XS*S;
         LBC=[X'+z*XS'+z^2*XS2',XS'*X+XS2'*XS+z*XS2'*X];    
      end 
      A=[P,P*XSz;LBC];
   end   
end

function [A]=formAz(L,vtau,z)
%
% Latest revision 23.04.2019
%
A=L{1,1};
N = size(A,1); 
for i=2:size(L,2)
    A=A+L{1,i}*exp(-vtau(i-1)*z);
end
A=A-eye(N)*z; 

function [v]=dAdz(L,vtau,z,x)
%
% Latest revision 10.03.2019
%
v=-L{1,2}*(vtau(1)*exp(-vtau(1)*z)*x);
for i=3:size(L,2)
    v=v-L{1,i}*(vtau(i-1)*exp(-vtau(i-1)*z)*x);
end
v=v-x; 

function [B]=formdAdz(L,vtau,z)
%
% Latest revision 10.03.2019
%
n=size(L{1,1},1);
B=-L{1,2}*(vtau(1)*exp(-vtau(1)*z));
for i=3:size(L,2)
    B=B-L{1,i}*(vtau(i-1)*exp(-vtau(i-1)*z));
end
B=B-eye(n); 
 
% function [resAXS]=resnorm(MATR,X,S)
% %
% % Latest revision 22.01.2018
% %
% if(size(S,1)>0)
%    % MATR=[L_0(:),L_tau(:),L_tau_A(:),tau,tau_A];
%    L_0=MATR(:,1:4);
%    L_tau=MATR(:,5:8);
%    L_tau_A=MATR(:,9:12);
%    tau=MATR(1,13);
%    tau_A=MATR(2,13);
%    resAXS=norm(L_0*X+L_tau*(X*expm(-tau*S))...
%                          +L_tau_A*(X*expm(-tau_A*S))-X*S,2);  
% else
%    resAXS=0.0; 
% end 