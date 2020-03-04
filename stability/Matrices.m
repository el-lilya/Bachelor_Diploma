function [L,CM,M]=Matrices(Uss,coef,delta,vtau,scheme)

vm=ceil(vtau/delta);
mp=vm(end);
mu_C = coef(6);
%

L=cell(1,(size(vtau,1)+2));
CM=cell(1,(vm(2)+2));
[U,Ui,F,Fi]=RHS(coef);
for i=1:size(U,2)
    L{1,i}=jacobian(F,U{1,i});
    for j=2:size(U,2)
        for k=1:size(U{1,1},1)
            L{1,i}=subs(L{1,i},U{1,j}(k),U{1,1}(k));
        end
    end
    for k=1:size(U{1,1},1)
        L{1,i}=subs(L{1,i},Ui(k),U{1,1}(k));
    end
    for l=1:size(U{1,1},1)
        L{1,i}=subs(L{1,i},U{1,1}(l),Uss(l));
    end
   L{1,i}=double(L{1,i});
end
i = i+1;
L{1,i}=jacobian(Fi,Ui);
for j=2:size(U,2)
    for k=1:size(U{1,1},1)
        L{1,i}=subs(L{1,i},U{1,j}(k),U{1,1}(k));
    end
end
for k=1:size(U{1,1},1)
    L{1,i}=subs(L{1,i},Ui(k),U{1,1}(k));
end
for l=1:size(U{1,1},1)
    L{1,i}=subs(L{1,i},U{1,1}(l),Uss(l));
end
L{1,i}=double(L{1,i});



I=eye(size(F,2));
if(scheme==1)
   C=(I-delta*L{1,1})\I;
   C1=C;
   C2=zeros(size(F,2),size(F,2)); 
elseif(scheme==2)    
   C=(1.5*I-delta*L{1,1}-delta^2/2*L{1,end})\I;
   C1=2*C;
   C2=-0.5*C;
end

m_C = vm(2);

C_int = zeros(7,7,m_C);
for i = 1:(m_C-1)
C_int(:,:,i) = C*delta^2*exp(-mu_C*delta*i)*L{1,end};
end  
%C_int(:,:,m_C) = C*delta^2/2*exp(-mu_C*delta*m_C)*L{1,2};%подозрительно
C_int(:,:,m_C) = C*delta^2/2*exp(-mu_C*delta*m_C)*L{1,end};

CM{1,1}=C1;
CM{1,2}=C2;
CM{1,3} = C_int(:,:,1:end-1);
CM{1,4} = C_int(:,:,end);
nov=size(U{1,1},1);
%
M=sparse(nov+1:1:nov*mp,1:1:nov*mp-nov,ones(nov*mp-nov,1),nov*mp,nov*mp,nov*nov*mp);

M(1:nov,1:nov)=sparse(C1+C_int(:,:,1));
M(1:nov,nov+1:2*nov)=sparse(C2+C_int(:,:,2));
M(1:nov,2*nov+1:m_C*nov) = sparse(reshape(C_int(:,:,3:end),nov,[]));

for i=2:size(U,2)
    Ci=C*(delta*L{1,i});
    mi=vm(i-1);
    M(1:nov,nov*mi-nov+1:nov*mi)=M(1:nov,nov*mi-nov+1:nov*mi)+sparse(Ci);
    CM{1,i+4}=[CM{1,i+4},Ci];
end
