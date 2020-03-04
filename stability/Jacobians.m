function L=Jacobians(Uss, coef)

A=Uss(1);
L=Uss(2);
I=Uss(3); 
V=Uss(4);
E0=Uss(5);
Q=Uss(6);
E=Uss(7);

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
w_E1=coef(31);
w_E2=coef(32);

rho_C = (1-p_L)*g_AV*A*V+g_AI*A*I+g_LI*L*I+a_L*L;

% int_V = rho_C/mu_C*(1-exp(-mu_C*w_C));
int_V = 0;

L0=zeros(7);
L0(1,1)=-mu_A-g_AV*V-g_AI*I;
L0(1,3)=-g_AI*A;
L0(1,4)=-g_AV*V;
L0(2,1)=p_L*g_AV*V;
L0(2,2)=-mu_L-a_L-g_LI*I;
L0(2,3)=-g_LI*L;
L0(2,4)=p_L*g_AV*A;
L0(3,3)=-mu_I-sigma_I*nu_I-g_EI*E;
L0(3,7)=-g_EI*I;
L0(4,1)=-g_AV*V;
L0(4,2)=-g_LV*V;
L0(4,4)=-mu_V-g_AV*A-g_LV*L-g_CV*int_V;
L0(5,1)=-g_E0AV*E0*V;
L0(5,4)=-g_E0AV*A*E0;
L0(5,5)=-mu_E0-g_E0AV*A*V;
L0(6,1)=-g_QAV*Q*V;
L0(6,4)=-g_QAV*Q*A;
L0(6,6)=-mu_Q-g_QAV*A*V;
L0(7,3)=-p_E*g_EI*E;
L0(7,7)=-mu_E-p_E*g_EI*I;

L1 = zeros(7);
% L1(4,1) = (1-p_L)*g_AV*V+g_AI*I;
% L1(4,2) = g_LI*I+a_L;
% L1(4,3) = g_AI*A+g_LI*L;
% L1(4,4) = (1-p_L)*g_AV*A;
% L1 = -g_CV*V*L1;

L_U = zeros(7);
L_U(4,3) = exp(-mu_U*w_U)*nu_I;

L_C = zeros(7);
L_C(3,1) = exp(-mu_C*w_C)*((1-p_L)*g_AV*V+g_AI*I);
L_C(3,2) = exp(-mu_C*w_C)*(g_LI*I+a_L);
L_C(3,3) = exp(-mu_C*w_C)*(g_AI*A+g_LI*L);
L_C(3,4) = exp(-mu_C*w_C)*(1-p_L)*g_AV*A;

L_E2 = zeros(7);
L_E2(6,1) = n_Q2*g_QAV*Q*V;
L_E2(6,4) = n_Q2*g_QAV*Q*A;
L_E2(6,6) = n_Q2*g_QAV*A*V;
L_E2(7,1) = n_E2*g_QAV*Q*V;
L_E2(7,4) = n_E2*g_QAV*Q*A;
L_E2(7,6) = n_E2*g_QAV*A*V;

L_E1 = zeros(7);
L_E1(6,1) = n_Q1*g_E0AV*E0*V;
L_E1(6,4) = n_Q1*g_E0AV*E0*A;
L_E1(6,5) = n_Q1*g_E0AV*A*V;
L_E1(7,1) = n_E1*g_E0AV*E0*V;
L_E1(7,4) = n_E1*g_E0AV*E0*A;
L_E1(7,5) = n_E1*g_E0AV*A*V;

L = cat(3,L0,L1,L_U,L_C,L_E2,L_E1);
L = num2cell(L,[1,2]);
