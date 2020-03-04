function Y=Function_with_int(U_int_full, delays, coef, delta)

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

A=U_int_full(1,end);
L=U_int_full(2,end);
I=U_int_full(3,end); 
V=U_int_full(4,end);
E0=U_int_full(5,end);
Q=U_int_full(6,end);
E=U_int_full(7,end);

U_w_C=delays(:,1); 
U_w_U=delays(:,2); 
U_w_E1=delays(:,3);
U_w_E2=delays(:,4); 

A_w_C = U_w_C(1);
L_w_C = U_w_C(2);
I_w_C = U_w_C(3);
V_w_C = U_w_C(4);
I_w_U = U_w_U(3);
E0_w_E1 = U_w_E1(5);
A_w_E1 = U_w_E1(1);
V_w_E1 = U_w_E1(4);
Q_w_E2 = U_w_E2(6);
A_w_E2 = U_w_E2(1);
V_w_E2 = U_w_E2(4);

rho_C_w_C = (1 - p_L) * g_AV * A_w_C * V_w_C + g_AI *  A_w_C * I_w_C +...
    g_LI * L_w_C * I_w_C + a_L * L_w_C;
rho_U_w_U = nu_I * I_w_U;
rho_E1_w_E1 = g_E0AV * E0_w_E1 * A_w_E1 * V_w_E1;
rho_E2_w_E2 = g_QAV * Q_w_E2 * A_w_E2 * V_w_E2;
rho_E1 = g_E0AV * E0 * A * V;
rho_E2 = g_QAV * Q * A * V;

A_int = U_int_full(1,:);
L_int = U_int_full(2,:);
I_int = U_int_full(3,:);
V_int = U_int_full(4,:);

t_s = w_C:-delta:0;

int = trapz(delta,exp(-mu_C*t_s).*((1-p_L)*g_AV*A_int.*V_int+...
    g_AI*A_int.*I_int+g_LI*L_int.*I_int+a_L*L_int));

Y=zeros(size(U_int_full,1),1);

Y(1) = r_A-mu_A*A-g_AV*A*V-g_AI*A*I;
Y(2) = -(mu_L+a_L)*L-g_LI*L*I+p_L*g_AV*A*V;
Y(3) = -(mu_I+sigma_I*nu_I)*I-g_EI*E*I+exp(-mu_C*w_C)*rho_C_w_C;
Y(4) = -mu_V*V-g_AV*A*V-g_LV*L*V-g_CV*int*V+exp(-mu_U*w_U)*rho_U_w_U;
Y(5) = r_E0-mu_E0*E0-rho_E1;
Y(6) = -mu_Q*Q-rho_E2+n_Q1*rho_E1_w_E1+n_Q2*rho_E2_w_E2;
Y(7) = -mu_E*E-p_E*g_EI*E*I+n_E1*rho_E1_w_E1+n_E2*rho_E2_w_E2;

end

