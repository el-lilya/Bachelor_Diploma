function Ploter(t,U) 
subplot(421)
   plot(t,U(1,:))
   xlabel('$t$','Interpreter','latex');           
   title('$A$','Interpreter','latex'); 
   hold on
subplot(422)
   plot(t,U(2,:))
   xlabel('$t$','Interpreter','latex');           
   title('$L$','Interpreter','latex'); 
   hold on
subplot(423) 
   plot(t,U(3,:))
   xlabel('$t$','Interpreter','latex');           
   title('$I$','Interpreter','latex'); 
   hold on
subplot(424)
   plot(t,U(4,:))
   xlabel('$t$','Interpreter','latex');           
   title('$V$','Interpreter','latex'); 
   hold on
subplot(425)
   plot(t,U(5,:))
   xlabel('$t$','Interpreter','latex');           
   title('$E_0$','Interpreter','latex'); 
   hold on
subplot(426)
   plot(t,U(6,:))
   xlabel('$t$','Interpreter','latex');           
   title('$Q$','Interpreter','latex'); 
   hold on
subplot(427)
   plot(t,U(7,:))
   xlabel('$t$','Interpreter','latex');           
   title('$E$','Interpreter','latex'); 
   hold on

   hold on
   
   
