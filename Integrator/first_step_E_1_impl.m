function first_step_E_1_impl

%first step with implicit Euler
T=1;  
delta = 0.5e-3;
u0 = 1;
%
err=1e-7; %Newton's error
%
inc=1000;  

%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
fprintf(1,'T=%g\n',T);
fprintf(1,'delta=%g\n',delta);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 


K=ceil(T/delta);

% 
t_u=(0:K)*delta;

u = zeros(1,size(t_u,2));
u(1) = u0;
u(2) = Step1(t_u(2),u(1),delta,err);
start_index_u = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
nns=0; 
for k=start_index_u:size(t_u,2) 
    [u(k),nnsk]=Step(t_u(k),u(k-1),u(k-2),delta,err);
    nns=nns+nnsk;
    if (mod(k-start_index_u+1,inc)==0) 
        fprintf(1,'t=%g, nns=%g\n',t_u(k),nns);

        nns=0;
    end
end
toc

% Plots

hold on
plot(t_u,u) 
end
% 

%

function [u,l]=Step(t,udelta,udelta2,delta,err)
l_max=5; 
u=udelta; 
l=0;
G=Res(t,u,udelta,udelta2,delta);
normG=norm(G,2);
while(normG>err&&l<l_max)
   l=l+1;
   u=u-Jacob(t,u,delta)\G;
   G=Res(t,u,udelta,udelta2,delta);
   normG=norm(G,2);
end
end


function G=Res(t,u,udelta,udelta2,delta)
G=(1.5*(u-udelta)+0.5*(udelta2-udelta))/delta;
G=G-f(t,u);
end

function f = f(t,u)
f = 5*t;
end  

function Jacob = Jacob(u,t,delta)
Jacob = 1.5/delta;
end  

function u=Step1(t,udelta,delta,err)
l_max=5; 
u=udelta; 
l=0;
G=Res1(t,u,udelta,delta);
normG=norm(G,2);
while(normG>err&&l<l_max)
   l=l+1;
   u=u-Jacob1(t,u,delta)\G;
   G=Res1(t,u,udelta,delta);
   normG=norm(G,2);
end
end


function G=Res1(t,u,udelta,delta)
G=(u-udelta)/delta;
G=G-f(t,u);
end


function Jacob = Jacob1(u,t,delta)
Jacob = 1./delta;
end  
