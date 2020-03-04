function [ncase,par]=Parameters2

ncase=1;
switch ncase
case 1
    r_A = 2000;   %1
    r_E0 = 120;   %2
    mu_A = 0.01;    %3
    mu_L = 0.01;    %4
    mu_I = 0.01;    %5
    mu_C = 0.01;    %6
    mu_E0 = 0.01;   %7
    mu_Q = 0.01;    %8
    mu_E = 0.08;    %9
    mu_U = 3.0;     %10
    mu_V = 3.0;     %11
    nu_I = 110.0;   %12
    sigma_I = 0.0032;%13
    p_L = 0.3;      %14
    a_L = 0.1;      %15
    w_C = 0.2;      %16
    w_U = 0.02;     %17
    n_Q1 = 12.0;    %18
    n_Q2 = 10.0;    %19
    n_E1 = 16.0;    %20
    n_E2 = 12.0;    %21
    p_E = 0.12;     %22
    g_LV = 0.35e-7; %23
    g_CV = 0.35e-7; %24
    g_AV = 0.35e-7; %25
    g_LI = 5e-7;    %26
    g_AI = 5e-7;    %27
    g_EI = 3.5e-5;  %28
    g_E0AV = 7e-10; %29
    g_QAV = 4e-10;  %30
    w_E1 = 1.5;     %31
    w_E2 = 1.2;     %32
    
    
end
%
par=[r_A;r_E0;mu_A;mu_L;mu_I;mu_C;mu_E0;mu_Q;mu_E;mu_U;mu_V;nu_I;sigma_I;
    p_L;a_L;w_C;w_U;n_Q1;n_Q2;n_E1;n_E2;p_E;g_LV;g_CV;g_AV;g_LI;g_AI;g_EI;
    g_E0AV;g_QAV;w_E1;w_E2];
%   
