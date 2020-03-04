fileID = fopen('st_st.txt','w');
numbers=[1,3,4];
for i=1:3
    [coef,vtau,~] = COEF(numbers(i));
    syms A L I V E0 Q E 
    
    r_A=coef(1);
    r_E0=coef(2);
    mu_A=coef(3);       
    mu_L=coef(4);        
    mu_I=coef(5);        
    mu_C=coef(6);
    mu_E0=coef(7);
    mu_Q=coef(8);
    mu_E=coef(9);
    mu_U=coef(10);
    mu_V=coef(11);
    nu_I=coef(12);
    sigma_I=coef(13);
    p_L=coef(14);
    a_L=coef(15);
    w_C=coef(16);
    w_U=coef(17);
    n_Q1=coef(18);
    n_Q2=coef(19);
    n_E1=coef(20);
    n_E2=coef(21);
    p_E=coef(22);
    g_LV=coef(23);
    g_CV=coef(24);
    g_AV=coef(25);
    g_LI=coef(26);
    g_AI=coef(27);
    g_EI=coef(28);
    g_E0AV=coef(29);
    g_QAV=coef(30);

    rho_E1 = g_E0AV * E0 * A * V;
    rho_E2 = g_QAV * Q * A * V;
    rho_U = nu_I * I;
    rho_C = (1-p_L)*g_AV*A*V+g_AI*A*I+g_LI*L*I+a_L*L;

    F(1) = r_A-mu_A*A-g_AV*A*V-g_AI*A*I;
    F(2) = -(mu_L+a_L)*L-g_LI*L*I+p_L*g_AV*A*V;
    F(3) = -(mu_I+sigma_I*nu_I)*I-g_EI*E*I+exp(-mu_C*w_C)*rho_C;
    F(4) = -mu_V*V-g_AV*A*V-g_LV*L*V-g_CV*rho_C/mu_C*(1-exp(-mu_C*w_C))*V+...
        exp(-mu_U*w_U)*rho_U;
    F(5) = r_E0-mu_E0*E0-rho_E1;
    F(6) = -mu_Q*Q-rho_E2+n_Q1*rho_E1+n_Q2*rho_E2;
    F(7) = -mu_E*E-p_E*g_EI*E*I+n_E1*rho_E1+n_E2*rho_E2;

    
    
    S = solve([F(1)==0,F(2)==0,F(3)==0,F(4)==0,F(5)==0,F(6)==0,F(7)==0]);
    %S = solve([F(1)==0,F(2)==0,F(3)==0,F(4)==0,F(5)==0,F(6)==0,F(7)==0],'Real',true);        
   
    
    res = [vpa(S.A)'; vpa(S.L)';vpa(S.I)'; vpa(S.V)'; vpa(S.E0)'; vpa(S.Q)';vpa(S.E)'];
    isreal(res)
    indexs = ones(1,size(res,2));
    for j =1:size(res,2)
        a=res(:,j);
        if a(res(:,j) < 0) ~= 0
                indexs(j)=0;
        end
    end
    res=res(:,indexs==1)    
    

    fprintf(fileID,'%d set of coefficients\n',i);
    fprintf(fileID,'%24s %24s %24s %24s %24s %24s %24s\n','A','L','I',...
        'V','E0','Q','E');
    fprintf(fileID,'%24.4g %24.4g %24.4g %24.4g %24.4g %24.4g %24.4g\n',res);
    fprintf(fileID,'\n');
%     res(:,2)
%     subs(F(1),[A,L,I,V,E0,Q,E],res(:,1)')
%     vpa(subs(F(1),[A,L,I,V,E0,Q,E],[4.956*10^6,1.139*10^3,1.650*10^2,5.604*10^3,1.005*10^3,8.900*10^6,...
%      1.482*10^6]))
%     log10(subs(F(1),[A,L,I,V,E0,Q,E],res(:,1)'))
    
    
end
fclose(fileID);


