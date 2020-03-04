function first_step_E_1

%first step with implicit Euler or constant
T=150;  
delta = 0.5e-3;
y0 = 1;
%
err=1e-7; %Newton's error
%
inc=100000;  

%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
fprintf(1,'T=%g\n',T);
fprintf(1,'delta=%g\n',delta);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% y - 1 step: constants
% u - 1 step: Implicit Euler


K=ceil(T/delta);

% 
t=(-1:K)*delta;

y = zeros(1,size(t,2));
y(1) = y0;
y(2) = y0;
start_index = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_u=(0:K)*delta;

u = zeros(1,size(t_u,2));
u(1) = y0;
u(2) = Step1(t_u(2),u(1),delta,err);
start_index_u = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
nns=0; 
for k=start_index:size(t,2) 
    [y(k),nnsk]=Step(t(k),y(k-1),y(k-2),delta,err);
    nns=nns+nnsk;
    if (mod(k-start_index+1,inc)==0) 
        fprintf(1,'t=%g, nns=%g\n',t(k),nns);

        nns=0;
    end
end
toc



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

norm(y(2:end)-u)/norm(u)
% Plots

figure('Position', [10,10,700,900]);

plot(t,y) 
hold on
plot(t_u,u) 
end
% 

%

function [y,l]=Step(t,ydelta,ydelta2,delta,err)
l_max=5; 
y=ydelta; 
l=0;
G=Res(t,y,ydelta,ydelta2,delta);
normG=norm(G,2);
while(normG>err&&l<l_max)
   l=l+1;
   y=y-Jacob(t,y,delta)\G;
   G=Res(t,y,ydelta,ydelta2,delta);
   normG=norm(G,2);
end
end


function G=Res(t,y,ydelta,ydelta2,delta)
G=(1.5*(y-ydelta)+0.5*(ydelta2-ydelta))/delta;
G=G-f(t,y);
end

function f = f(t,y)
f = -y;
end  

function Jacob = Jacob(t,y,delta)
Jacob = 1.5/delta+1;
end  

function u=Step1(t,udelta,delta,err)
l_max=5; 
u=udelta; 
l=0;
G=Res1(t,u,udelta,delta);
normG=norm(G,2);
while(normG>err&&l<l_max)
   l=l+1;
   u=u-Jacob_u(t,u,delta)\G;
   G=Res1(t,u,udelta,delta);
   normG=norm(G,2);
end
end


function G=Res1(t,u,udelta,delta)
G=(u-udelta)/delta;
G=G-f(t,u);
end
function Jacob = Jacob_u(t,u,delta)
Jacob = 1./delta;
end  
