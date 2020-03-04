function Test_stability
title = 'steady_state3';
T=150;  
delta = 0.5e-3;
%delta = 1e-4;

scheme=2; 
%
err=1e-7; %Newton's error
%
inc=10000;  

tUfilename='tU'; 
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

U_full=zeros(7,size(t,2));

%Uss = [199929., 6.84703, 0.833115, 28.5402, 8574.97, 13952., 11626.2];
%Uss = [200000,0,0,0,12000,0,0];
Uss = [5000000,0,0,0,1200,0,0];
Uss = [4.956*10^6;1.388*10^3;1.650*10^2;5.604*10^3;1.005*10^3;8.897*10^6;...
     1.482*10^6];
U_full(1,1:m_max+1)=Uss(1);
U_full(2,1:m_max+1)=Uss(2);
U_full(3,1:m_max+1)=Uss(3);
U_full(4,1:m_max+1)=Uss(4);
U_full(5,1:m_max+1)=Uss(5);
U_full(6,1:m_max+1)=Uss(6);
U_full(7,1:m_max+1)=Uss(7);


% max+2 - min index: t>0 
start_index=m_max+2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
nns=0; 
for k=start_index:size(t,2) 
    delays=[U_full(:,k-m(1)),U_full(:,k-m(2)),U_full(:,k-m(3)),U_full(:,k-m(4))];
    U_int = U_full(:,k-m(1):k-1); 
    [U_full(:,k),nnsk]=Step(U_int, delays,coef,delta,err,scheme);
    nns=nns+nnsk;
    if (mod(k-start_index+1,inc)==0) 
        fprintf(1,'t=%g, nns=%g\n',t(k),nns);

        nns=0;
    end
end
toc

% Plots

%Pertzev's results from excel
Table = readtable('Результат_теста_3.xls','Range','S6:Y157');
size(Table);
U_Pertz = zeros(size(Table,2),size(Table,1));
U_Pertz(1,:)=Table.Var1;
U_Pertz(2,:)=Table.Var2;
U_Pertz(3,:)=Table.Var3;
U_Pertz(4,:)=Table.Var4;
U_Pertz(5,:)=Table.Var5;
U_Pertz(6,:)=Table.Var6;
U_Pertz(7,:)=Table.Var7;
t_Pertz = 0:1:(size(U_Pertz,2)-1);

figure('Name', title,'Position', [10,10,700,900]);

Ploter(t_Pertz,U_Pertz)
Ploter(t,U_full) 

% 
U = U_full(:,start_index-1:1000:end);
save(tUfilename,'t','U');
%

function [U,l]=Step(U_int,delays,coef,delta,err,scheme)
l_max=5; 
U=U_int(:,end); 
l=0;
G=Res(U,U_int,delays,coef,delta,scheme);
normG=norm(G,2);
while(normG>err&&l<l_max)
   l=l+1;
   %U=U-NumJacob(U,U_int,delays, coef,delta,scheme)\G;
   U=U-Jacob1(U,U_int,coef,delta,scheme)\G;
   G=Res(U,U_int,delays,coef,delta,scheme);
   normG=norm(G,2);
end



function G=Res(U,U_int,delays,coef,delta,scheme)
Udelta = U_int(:,end);
Udelta2 = U_int(:,end-1);

if(scheme==1)
   G=(U-Udelta)/delta; 
elseif(scheme==2)    
   G=(1.5*(U-Udelta)+0.5*(Udelta2-Udelta))/delta;
end
G=G-Function_with_int([U_int,U],delays,coef,delta);
  
