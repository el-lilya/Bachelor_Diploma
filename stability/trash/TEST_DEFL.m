function TEST_DEFL

Uss = cell(1,4);
Uss{1,1} = [5*10^6;0;0;0;1.2*10^3;0;0];
Uss{1,2} = [4.956*10^6;1.139*10^3;1.650*10^2;5.604*10^3;1.005*10^3;8.900*10^6;...
     1.482*10^6];
Uss{1,3} = [200000,0,0,0,12000,0,0];
Uss{1,4} = [200000,0,0,0,12000,0,0];
Uss{1,5} = [199929., 6.84703, 0.833115, 28.5402, 8574.97, 13952., 11626.2];


DEFL = 0;
fileID = fopen('test_DEFL.txt','w');
for deflation=DEFL
fprintf(fileID,'Deflation = %d\n',deflation);
    for i = 2
        [coef,vtau,~] = COEF(i);
        [E0,~,Ei,~,resEi,min_max,l,sing_1,sing_2,sing_3,rank]=EigProbSol_light_DEFL(Uss{1,i},coef,vtau,15,i,0,deflation);
        data = [E0(1:7)';Ei';round(log10(resEi));round(log10(min_max));l;round(log10(sing_1));round(log10(sing_2));round(log10(sing_3));rank];
        fprintf(fileID,'%d steady state\n',i);
        fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s %12s %12s\n' ,'E0','Ei','res','1/cond','l','sing_1','sing_2','sing_3','rank');
        fprintf(fileID,'%12f %12f %12d %12d %12d %12d %12d %12d %12d\n',data);
        fprintf(fileID,'\n');
    end
end
fclose(fileID);


% % All
% % 
% fileID = fopen('test1.txt','w');
% for i=1:5
% [coef,vtau,~] = COEF(i);
% fprintf(fileID,'%d steady state\n',i);
% fprintf(fileID,'\n');
% data = []; res = []; sing = [];
% % data = E0,Ei,resEi,sigma_min/sigma_max
%     for deflation = [0,2]
%         [E0,resE0,Ei,~,resEi,min_max]=EigProbSol_light_DEFL(Uss{1,i},coef,vtau,15,i,0,deflation);
%         data = [data;Ei';round(log10(resEi))];
%         resE0 = round(log10(resE0));
%         %res = [res;round(log10(resEi))];
%         sing = [sing;round(log10(min_max))];
%     end
%     fprintf(fileID,'%12s %6s\n','E0','res0');
%     fprintf(fileID,'%12f %6d\n',[E0(1:7)';resE0(1:7)]);     
%     fprintf(fileID,'\n'); 
%     fprintf(fileID,'\n');
%     fprintf(fileID,'%12s \n','Deflation');
%     fprintf(fileID,'%12d %18d \n',[0,2]); 
%     fprintf(fileID,'\n'); 
%     fprintf(fileID,'%12s %6s\n','Ei','resEi');
%     fprintf(fileID,'%12f %6d %12f %6d \n',data); 
%     fprintf(fileID,'\n'); 
%     fprintf(fileID,'\n');
%     fprintf(fileID,'%12s \n','1/rcond');
%     fprintf(fileID,'%12d %18d \n',sing); 
%     fprintf(fileID,'\n'); 
%     fprintf(fileID,'\n');
% end
% fclose(fileID);