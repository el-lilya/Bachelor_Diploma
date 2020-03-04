function TEST
%
% Последнее обновление 09.10.2019 
%
% Численное интегрирование модели ВИЧ-инфекци Бочарова-Перцева из работы
% [Pertsev N. V, Loginov K. К., Bocharov G. A. Nonlinear effects 
% in the dynamics of HIV-1 infection predicted by mathematical model 
% with multiple delays // X. X. Vol. X. No. X, P. X—X].
%
% Авторы: Саржина Лилия
%
% Внешние процедуры: Parameters, Ploter, NumJacob
% Внутренние процедуры: Step, Res
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% УСТАНОВКА ПАРАМЕТРОВ:
%
T=10;  % время интегрирования
delta=1.0e-3; % шаг сетки 

scheme=2; % схема интегрирования: 1 - неявный метод Эйлера,
          %                       2 - BDF2.
err=1.0e-7; % точность метода Ньютона  
%
iclean=1;  % предварительная чистка картинок
%
marks={'r-','g-','b-'}; % маркеры для отрисовки результатов интегрирования
mark=marks{3}; % текущий маркер 
%
inc=100; % инкремент выдачи информации на экран 

nfig=1; % номер фигуры (если nfig<1, то картинка не рисуется).
%
tUfilename='tU'; % файл для результатов расчета t, U
%
iresult=0; % 0 - начало расчета,
           % Продолжение расчета: 
           % 1 - дополнить старое решение новыми данными,    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% вывод значений параметров на экран:
fprintf(1,'T=%g\n',T);
fprintf(1,'delta=%g\n',delta);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ИНИЦИАЛИЗАЦИЯ:
%
% параметры 
[~,coef]=Parameters;

%массив задержек w_C, w_U, w_E1, w_E2
w=[coef(16),coef(17),coef(31),coef(32)];
%массив дискретных аналогов задержек
m=ceil(w/delta);
m_max=max(m);

K=ceil(T/delta);
if(iresult==0) 
% сетка:
t=(-m_max:1:K)*delta;
% Параметры, используемые для задания начальных данных 



% Задание начальных данных
t_V = 0.01;
m_V = ceil(t_V/delta);
U_full=zeros(7,size(t,2));
U_full(1,1:m_max+1)=coef(1)/coef(3);
U_full(2,1:m_max+1)=0;
U_full(3,1:m_max+1)=0;
U_full(4,1:m_max-m_V)=0;
U_full(5,1:m_max+1)=coef(2)/coef(7);
U_full(6,1:m_max+1)=0;

U_full(4,m_max-m_V+1:m_max+1)=5000*(t(m_max-m_V+1:m_max+1)+0.01);
U_full(7,m_max+1)=0;


% max+2 - индекс, отвечающий минимальному узлу сетки такому, что t>0 
start_index=m_max+2;
else
    % ПРОДОЛЖЕНИЕ РАСЧЁТА       
    load(tUfilename); 
    start_index=size(t,2)+1;
    t_new=(1:1:K)*delta+t(end);
    U_full=[U_full,zeros(7,size(t_new, 2))];
    t=[t,t_new];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Нелинейное интегрирование:
tic
nns=0; 
% max+2 - индекс, отвечающий минимальному узлу сетки такому, что t>0 
for k=start_index:size(t,2) 
    delays=[U_full(:,k-m(1)),U_full(:,k-m(2)),U_full(:,k-m(3)),U_full(:,k-m(4))];
    U_int = U_full(:,k-m(1):k-1); %%возможно +1
    [U_full(:,k),nnsk]=Step(U_int, delays,coef,delta,err,scheme);
    nns=nns+nnsk;
    if (mod(k-start_index+1,inc)==0) % информация о работе метода
        fprintf(1,'t=%g, nns=%g\n',t(k),nns);
        
        nns=0;
    end
end
toc
U_full(:,start_index-1:inc:k);
%Отрисовка графиков 
if (nfig>0)
   figure(nfig)
   Ploter(t,U_full,mark,iclean) 
end

% Запись t, U в файл:
save(tUfilename,'t','U');
%

function [U,l]=Step(U_int,delays,coef,delta,err,scheme)
%
% Последнее обновление 24.09.2018
%
% Шаг решения задачи Коши.
%
% Авторы: Саржина Лилия
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ВХОДНЫЕ ПАРАМЕТРЫ:
% Udelta,Udelta2,delays - заданные значения в узлах сетки
% U_int - значения U от t-w_C до t, не включая t (U_int = U(:,k-m(1):k-1))
% coef - коэффициенты модели
% delta - шаг сетки
% err - точность решения нелинейного уравнения, к которому сводиться
%       реализация шага метода 
% scheme=1 - неявный метод Эйлера 
%        2 - BDF2
% ВЫХОДНЫЕ ПАРАМЕТРЫ:
% U - значение в следующем узле сетки
% l - количество шагов, сделанных методом Ньютона  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l_max=10; %максимально возможное количество шагов метода Ньютона 

U=U_int(:,end); %задаем начальное приближение 
l=0;
G=Res(U,U_int,delays,coef,delta,scheme);
normG=norm(G,2);
while(normG>err&&l<l_max)
   l=l+1;
   U=U-NumJacob(U,U_int,delays, coef,delta,scheme)\G;
   %NumJacob(U,U_int,delays, coef,delta,scheme)\G;
   G=Res(U,U_int,delays,coef,delta,scheme);
   normG=norm(G,2);
end

function G=Res(U,U_int,delays,coef,delta,scheme)
%
% Последнее обновление 24.09.2018
%
% Считает невязку
%
% Авторы: Саржина Лилия
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ВХОДНЫЕ ПАРАМЕТРЫ:
% U, Udelta,Udelta2,delays - заданные  значения
% U - значения U от t-w_C до t
% coef - коэффициенты модели
% delta - шаг сетки
% scheme=1 - неявный метод Эйлера 
%        2 - BDF2
% ВЫХОДНЫЕ ПАРАМЕТРЫ:
% G - значение невязки 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Udelta = U_int(:,end);
Udelta2 = U_int(:,end-1);

if(scheme==1)
   G=(U-Udelta)/delta; 
elseif(scheme==2)    
   G=(1.5*(U-Udelta)+0.5*(Udelta2-Udelta))/delta;
end
G=G-Function_with_int([U_int,U],delays,coef,delta);
  
% % function dGdU=Jacob(U,U_int,coef,delta,scheme)
% % 
% % Последнее обновление 24.09.2018
% % 
% % Считает Якобиан системы
% % 
% % Авторы: Саржина Лилия
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ВХОДНЫЕ ПАРАМЕТРЫ:
% % U,delays - заданные  значения
% % U - значения U от t-w_C до t
% % coef - коэффициенты модели
% % delta - шаг сетки
% % scheme=1 - неявный метод Эйлера 
% %        2 - BDF2
% % ВЫХОДНЫЕ ПАРАМЕТРЫ:
% % dGdU - Якобиан системы 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U_int = [U_int,U];
% A=U(1);
% L=U(2);
% I=U(3); 
% V=U(4);
% E0=U(5);
% Q=U(6);
% E=U(7);
% 
% r_A=coef(1);
% r_E0=coef(2);
% mu_A=coef(3);       
% mu_L=coef(4);        
% mu_I=coef(5);        
% mu_C=coef(6);
% mu_E0=coef(7);
% mu_Q=coef(8);
% mu_E=coef(9);
% mu_U=coef(10);
% mu_V=coef(11);
% nu_I=coef(12);
% sigma_I=coef(13);
% p_L=coef(14);
% a_L=coef(15);
% w_C=coef(16);
% w_U=coef(17);
% n_Q1=coef(18);
% n_Q2=coef(19);
% n_E1=coef(20);
% n_E2=coef(21);
% p_E=coef(22);
% g_LV=coef(23);
% g_CV=coef(24);
% g_AV=coef(25);
% g_LI=coef(26);
% g_AI=coef(27);
% g_EI=coef(28);
% g_E0AV=coef(29);
% g_QAV=coef(30);
% w_E1=coef(31);
% w_E2=coef(32);
% 
% 
% 
% A_int = U_int(1,:);
% L_int = U_int(2,:);
% I_int = U_int(3,:);
% V_int = U_int(4,:);
% 
% t = w_C:-delta:0;
% int_A = trapz(delta,exp(-mu_C*t).*((1-p_L)*g_AV*V_int+g_AI*I_int));
% int_L = trapz(delta,exp(-mu_C*t).*(g_LI*I_int+a_L));
% int_I = trapz(delta,exp(-mu_C*t).*(g_AI*A_int+g_LI*L_int));
% int_V_1 = trapz(delta,exp(-mu_C*t).*((1-p_L)*g_AV*A_int.*V_int+g_AI*A_int.*I_int+g_LI*L_int.*I_int+a_L*L_int));
% int_V_2 = trapz(delta,exp(-mu_C*t).*((1-p_L)*g_AV*A_int));
% 
% if(scheme==1)
%    dGdU=(1.0/delta)*eye(7); 
% elseif(scheme==2)
%    dGdU=(1.5/delta)*eye(7);
% end
% d=zeros(7);
% 
% d(1,1)=-mu_A-g_AV*V-g_AI*I;
% d(1,3)=-g_AI*A;
% d(1,4)=-g_AV*V;
% d(2,1)=p_L*g_AV*V;
% d(2,2)=-mu_L-a_L-g_LI*I;
% d(2,3)=-g_LI*L;
% d(2,4)=p_L*g_AV*A;
% d(3,3)=-mu_I-sigma_I*nu_I-g_EI*E;
% d(3,7)=-g_EI*I;
% d(4,1)=-g_AV*V-g_CV*V*int_A;
% d(4,2)=-g_LV*V-g_CV*V*int_L;
% d(4,3)=-g_CV*V*int_I;
% d(4,4)=-mu_V-g_AV*A-g_LV*L-g_CV*int_V_1-g_CV*V*int_V_2;
% d(5,1)=-g_E0AV*E0*V;
% d(5,4)=-g_E0AV*A*E0;
% d(5,5)=-mu_E0-g_E0AV*A*V;
% d(6,1)=-g_QAV*Q*V;
% d(6,4)=-g_QAV*Q*A;
% d(6,6)=-mu_Q-g_QAV*A*V;
% d(7,3)=-p_E*g_EI*E;
% d(7,7)=-mu_E-p_E*g_EI*I;
% 
% dGdU=dGdU-d;









