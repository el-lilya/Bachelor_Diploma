function [U,Ui,F,Fi]=RHS(coef)

syms A L I V E0 Q E ...
    A_w_C L_w_C I_w_C V_w_C E0_w_C Q_w_C E_w_C   ...
    A_w_U L_w_U I_w_U V_w_U E0_w_U Q_w_U E_w_U ...
    A_w_E1 L_w_E1 I_w_E1 V_w_E1 E0_w_E1 Q_w_E1 E_w_E1...
    A_w_E2 L_w_E2 I_w_E2 V_w_E2 E0_w_E2 Q_w_E2 E_w_E2...
    A_s L_s I_s V_s E0_S Q_s E_s
    


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

U=cell(1,5);
U{1,1}=[A; L; I; V; E0; Q; E];
U{1,2}=[A_w_U; L_w_U; I_w_U; V_w_U; E0_w_U; Q_w_U; E_w_U];
U{1,3}=[A_w_C; L_w_C; I_w_C; V_w_C; E0_w_C; Q_w_C; E_w_C];
U{1,4}=[A_w_E2; L_w_E2; I_w_E2; V_w_E2; E0_w_E2; Q_w_E2; E_w_E2];
U{1,5}=[A_w_E1; L_w_E1; I_w_E1; V_w_E1; E0_w_E1; Q_w_E1; E_w_E1];

Ui = [A_s; L_s; I_s; V_s; E0_S; Q_s; E_s];

rho_C_w_C = (1 - p_L) * g_AV * A_w_C * V_w_C + g_AI *  A_w_C * I_w_C +...
    g_LI * L_w_C * I_w_C + a_L * L_w_C;
rho_U_w_U = nu_I * I_w_U;
rho_E1_w_E1 = g_E0AV * E0_w_E1 * A_w_E1 * V_w_E1;
rho_E2_w_E2 = g_QAV * Q_w_E2 * A_w_E2 * V_w_E2;
rho_E1 = g_E0AV * E0 * A * V;
rho_E2 = g_QAV * Q * A * V;

rho_C = (1-p_L)*g_AV*A*V+g_AI*A*I+g_LI*L*I+a_L*L;

% 
F(1) = r_A-mu_A*A-g_AV*A*V-g_AI*A*I;
F(2) = -(mu_L+a_L)*L-g_LI*L*I+p_L*g_AV*A*V;
F(3) = -(mu_I+sigma_I*nu_I)*I-g_EI*E*I+exp(-mu_C*w_C)*rho_C_w_C;
F(4) = -mu_V*V-g_AV*A*V-g_LV*L*V-g_CV*rho_C/mu_C*(1-exp(-mu_C*w_C))*V+...
    exp(-mu_U*w_U)*rho_U_w_U;
F(5) = r_E0-mu_E0*E0-rho_E1;
F(6) = -mu_Q*Q-rho_E2+n_Q1*rho_E1_w_E1+n_Q2*rho_E2_w_E2;
F(7) = -mu_E*E-p_E*g_EI*E*I+n_E1*rho_E1_w_E1+n_E2*rho_E2_w_E2;


rho_C_s = (1-p_L)*g_AV*A_s*V_s+g_AI*A_s*I_s+g_LI*L_s*I_s+a_L*L_s;
Fi = sym(zeros(7,1));
Fi(4) = -g_CV*V*rho_C_s;
