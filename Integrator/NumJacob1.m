function dGdU=NumJacob1(U, U_int, delays, coef,delta,scheme,eps)
k = size(U,1);
if(scheme==1)
   dGdU=(1.0/delta)*eye(k); 
elseif(scheme==2)
   dGdU=(1.5/delta)*eye(k);
end
I=eye(k);
epsilon = eps;
for i=1:k
    dGdU(:,i)= dGdU(:,i)-(Function_with_int([U_int,U+...
        epsilon*I(:,i)],delays,coef,delta)-...
        Function_with_int([U_int,U-epsilon*I(:,i)],delays,coef,delta))/...
        (2*epsilon);
end

end

