function TEST

Uss = cell(1,4);
Uss{1,1} = [5*10^6;0;0;0;1.2*10^3;0;0];
Uss{1,2} = [4.956*10^6;1.388*10^3;1.650*10^2;5.604*10^3;1.005*10^3;8.897*10^6;...
     1.482*10^6];
Uss{1,3} = [200000,0,0,0,12000,0,0];
Uss{1,4} = [200000,0,0,0,12000,0,0];
Uss{1,5} = [199929., 6.84703, 0.833115, 28.5402, 8574.97, 13952., 11626.2];

fileID = fopen('test.txt','w');

for i = 1:5
    [coef,vtau,~] = COEF(i);
    [E0,Ei,s1,s2,s3,s4]=EigProbSol_light(Uss{1,i},coef,vtau,15,i,1);
    data = [E0(1:7)';Ei';round(log10(s1));round(log10(s2));...
        round(log10(s3));round(log10(s4))];
    fprintf(fileID,'%d steady state\n',i);
    fprintf(fileID,'%12s %12s %12s %12s %12s %12s\n','E0','Ei','s1',...
        's2','s3','s4');
    fprintf(fileID,'%12f %12f %12d %12d %12d %12d\n',data);
    fprintf(fileID,'\n');
end

fclose(fileID);

