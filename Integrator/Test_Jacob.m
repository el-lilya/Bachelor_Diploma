function Test_Jacob

T=150;  
delta=1e-4;

scheme=2; 
%
% err=1.0e-7; 
%

%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
fprintf(1,'T=%g\n',T);
fprintf(1,'delta=%g\n',delta);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
[~,coef]=Parameters1;

%w_C, w_U, w_E1, w_E2
w=[coef(16),coef(17),coef(31),coef(32)];

% discrete delays
m=ceil(w/delta);
m_max=max(m);

K=ceil(T/delta);

% 
t=(-m_max:1:K)*delta;


%
t_V = 0.01;
m_V = ceil(t_V/delta);
U_full=zeros(7,size(t,2));
U_full(1,1:m_max+1)=coef(1)/coef(3);
U_full(2,1:m_max+1)=100;
U_full(3,1:m_max+1)=100;
U_full(4,1:m_max-m_V)=1000;
U_full(5,1:m_max+1)=coef(2)/coef(7);
U_full(6,1:m_max+1)=100;

U_full(4,m_max-m_V+1:m_max+1)=5000*(t(m_max-m_V+1:m_max+1)+0.01);
U_full(7,m_max+1)=1000;

% U_full(1,1:m_max+1)=exp(12.2058)-1;
% U_full(2,1:m_max+1)=exp(2)-1;
% U_full(3,1:m_max+1)=exp(0.5872)-1;
% U_full(4,1:m_max-m_V)=exp(3.3452)-1;
% U_full(5,1:m_max+1)=exp(9.1195)-1;
% U_full(6,1:m_max+1)=exp(9.5289)-1;
% 
% U_full(4,m_max-m_V+1:m_max+1)=exp(3.3452)-1+15000*(t(m_max-m_V+1:m_max+1)+0.01);
% U_full(7,m_max+1)=exp(9.3488)-1;


% max+2 - min index: t>0 
start_index=m_max+2;
k = start_index;
delays=[U_full(:,k-m(1)),U_full(:,k-m(2)),U_full(:,k-m(3)),U_full(:,k-m(4))];
U_int = U_full(:,k-m(1):k-1); 
%U = [4999998.6397854; 0.709636004673691; 0.709636004673691; 17.9010576865359; 0; 0; 0];
%U = [206880.404079681; 1718602.88568155; 1718602.88568155; ...
%    65012290.4874993; 126128.355941138; 75149.7782729832; 126128.355941138];
%U = [677176.129646076; 2312.84496046766; 2312.84496046766; 79784.5220695454; 1333161.81216673; 2213718.13422443; 1333161.81216673];
%U = U_int(:,end);
U = [199998.72427773;	0.356391207464418;	0.460182233524496;	14.7970233893796;	11913.5183045234;	630.624778251331;	805.019730406017];

J1=Jacob1(U,U_int, coef,delta,scheme);
J = Jacob(U,U_int, coef,delta,scheme);
norm(J1-J)
epsilons = [100000,1.0e-6,1.0e-7,1.0e-8, 1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,1e-16];
%epsilons = 1.0e-7;
norm_dif = zeros(1,size(epsilons,2));
for i = 1:size(epsilons,2)
    eps = epsilons(i);
    J_num=NumJacob1(U,U_int,delays, coef,delta,scheme,eps);
    norm_dif(i) = norm(J1-J_num);
end 
norm_J = norm(J);
norm_J1 = norm(J1);
norm_J_num = norm(J_num);
norm_dif
% J_num-J1;
% (J_num-J1)./J_num;
% J_num
% J1
% J_num=NumJacob(U,U_int,delays, coef,delta,scheme)-(1.5/delta)*eye(7);
% %J=Jacob(U,U_int, coef,delta,scheme)-(1.5/delta)*eye(7)
% J1=Jacob1(U,U_int, coef,delta,scheme)-(1.5/delta)*eye(7)

