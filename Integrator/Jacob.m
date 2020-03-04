function dGdU=Jacob(U,U_int, coef,delta,scheme)

U_int = [U_int,U];
A=U(1);
L=U(2);
I=U(3); 
V=U(4);
E0=U(5);
Q=U(6);
E=U(7);

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



A_int = U_int(1,:);
L_int = U_int(2,:);
I_int = U_int(3,:);
V_int = U_int(4,:);

t = w_C:-delta:0;
int_A = trapz(delta,exp(-mu_C*t).*((1-p_L)*g_AV*V_int+g_AI*I_int));
int_L = trapz(delta,exp(-mu_C*t).*(g_LI*I_int+a_L));
int_I = trapz(delta,exp(-mu_C*t).*(g_AI*A_int+g_LI*L_int));
int_V_1 = trapz(delta,exp(-mu_C*t).*((1-p_L)*g_AV*A_int.*V_int+g_AI*A_int.*I_int+g_LI*L_int.*I_int+a_L*L_int));
int_V_2 = trapz(delta,exp(-mu_C*t).*((1-p_L)*g_AV*A_int));

if(scheme==1)
   dGdU=(1.0/delta)*eye(7); 
elseif(scheme==2)
   dGdU=(1.5/delta)*eye(7);
end
d=zeros(7);

d(1,1)=-mu_A-g_AV*V-g_AI*I;
d(1,3)=-g_AI*A;
d(1,4)=-g_AV*V;
d(2,1)=p_L*g_AV*V;
d(2,2)=-mu_L-a_L-g_LI*I;
d(2,3)=-g_LI*L;
d(2,4)=p_L*g_AV*A;
d(3,3)=-mu_I-sigma_I*nu_I-g_EI*E;
d(3,7)=-g_EI*I;
d(4,1)=-g_AV*V-g_CV*V*int_A;
d(4,2)=-g_LV*V-g_CV*V*int_L;
d(4,3)=-g_CV*V*int_I;
d(4,4)=-mu_V-g_AV*A-g_LV*L-g_CV*int_V_1-g_CV*V*int_V_2;
d(5,1)=-g_E0AV*E0*V;
d(5,4)=-g_E0AV*A*E0;
d(5,5)=-mu_E0-g_E0AV*A*V;
d(6,1)=-g_QAV*Q*V;
d(6,4)=-g_QAV*Q*A;
d(6,6)=-mu_Q-g_QAV*A*V;
d(7,3)=-p_E*g_EI*E;
d(7,7)=-mu_E-p_E*g_EI*I;

dGdU=dGdU-d;
