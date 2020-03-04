function dGdU=NumJacob_poisk_eps(U, U_int, delays, coef,delta,scheme)

epsilons = [1.0e-7,1.0e-8, 1e-9,1e-10,1e-11,1e-12,1e-14]; 
k=size(U,1);
if(scheme==1)
   %dGdU=(1.0/delta)*eye(k); 
elseif(scheme==2)
   %dGdU=(1.5/delta)*eye(k);
   dGdU=zeros(k);
   dGdU_eps = zeros((k+1)*size(epsilons,2),k);
end
I=eye(k);

for j=1:size(epsilons,2)
   epsilon = epsilons(j);
   dGdU_eps((j-1)*(k+1)+1,:)=log10(epsilon);
    for i=1:k
        dGdU(:,i)= dGdU(:,i)-(Function_with_int([U_int,U+epsilon*I(:,i)],delays,coef,delta)-Function_with_int([U_int,U-epsilon*I(:,i)],delays,coef,delta))/(2*epsilon);
    end
   dGdU_eps((j-1)*(k+1)+2:j*(k+1),:) = dGdU;
end
save('dGdU.mat', 'dGdU_eps');
end

