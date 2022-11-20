%% Data simulation
clear
clc
rng(1)
X=readmatrix('FinalMatrix_100.txt');
r = 2; % rank of the data matrix X
m = size(X,1); n = size(X,2); % size of the data matrix X (m x n)


Gtheo = round(255.*rand(m,r));
Ftheo = rand(r,n);

% Xtheo = Gtheo*Ftheo; % simulation of theoretical data matrix Xtheo

N = zeros(m,n); % simulating noise matrice N
% X = Xtheo+N; % simulating data matrix X


%% Matrix initialisation
Ginit = round(255.*rand(m,r));
Finit = rand(r,n);


%% NMF parameters
Iter_max = 1.e4;
rho_G = .01;
rho_F = .01;

%% Simulating matrices for weighted methods (missing entry example)
W = ones(m,n);
prop_missing = .1; % proportion of missing entries in X
idx_missing = randperm(m*n);
W(idx_missing(1:round(prop_missing*m*n))) = 0;
X = X.*W;

%% Running EM_WMU_NMF
Iter_max_E_step=2;
Iter_max_M_step=5;
Nesterov_Max_Iterations=100;
Nesterov_Min_Iterations=10;
% WMU_NMF
% fprintf('\nRunning NMF Lee and Seung multiplicative update...');
% [ G_WMU_NMF , F_WMU_NMF ] = WMU_NMF( W , X , Ginit , Finit , Iter_max );

% EM_WMU_NMF
fprintf('\n Running EM_WMU_NMF...');
% [EM_WMU_NMF,G_EM_WMU_NMF,F_EM_WMU_NMF,T_EM_WMU_NMF,I_EM_WMU_NMF,f_EM_WMU_NMF]  = WNMFLibrary.EM_WMU_NMF( W , X , Ginit , Finit , Iter_max,Iter_max_E_step );

% EM_WNE_NMF
fprintf('\nRunning EM_WNE_NMF...');
[EM_WNE_NMF,G_EM_WNE_NMF,F_EM_WNE_NMF,T_EM_WNE_NMF,I_EM_WNE_NMF,f_EM_WNE_NMF]  = WNMFLibrary.EM_WNE_NMF( W , X , Ginit , Finit , Iter_max_M_step,Iter_max_E_step,Nesterov_Max_Iterations,Nesterov_Min_Iterations );



% EM_L1R_WNE_NMF_
fprintf('\nRunning EM_L1R_WNE_NMF_...');
[L1R_WNE_NMF, G_EM_L1R_WNE_NMF,F_EM_L1R_WNE_NMF, T_EM_L1R_WNE_NMF , I_EM_L1R_WNE_NMF , f_EM_L1R_WNE_NMF] = WNMFLibrary.EM_L1R_WNE_NMF( W , X , G_EM_WNE_NMF , F_EM_WNE_NMF ,Iter_max_M_step,Iter_max_E_step,Nesterov_Max_Iterations,Nesterov_Min_Iterations  );

% EM_GWNE_NMF
% fprintf('\nRunning EM_GWNE_NMF...');
% EM_GWNE_NMF  = WNMFLibrary.EM_GWNE_NMF( W , X , Ginit , Finit ,Iter_max_M_step,Iter_max_E_step,Nesterov_Max_Iterations,Nesterov_Min_Iterations  );

fprintf ('\n ### Matrix size is %d x %d ' ,m,n);
fprintf (' \n ### Nesterov_Max_Iterations= % d ,  Nesterov_Min_Iterations = %d, Iter_max_E_step = %d , Iter_max_M_step= %d ',Nesterov_Max_Iterations,Nesterov_Min_Iterations,Iter_max_E_step,Iter_max_M_step);
% fprintf (' \n ### EM_WMU_NMF, Time = %d , Initial value = %d, Objective value = %d  ',real(T_EM_WMU_NMF),real(I_EM_WMU_NMF),real(f_EM_WMU_NMF)); 
fprintf (' \n ### EM_WNE_NMF, Time = %d , Initial value = %d, Objective value = %d  ',real(T_EM_WNE_NMF),real(I_EM_WNE_NMF),real(f_EM_WNE_NMF)); 
fprintf (' \n ### EM_L1R_WNE_NMF, Time = %d , Initial value = %d, Objective value = %d  ',real(T_EM_L1R_WNE_NMF),real(I_EM_L1R_WNE_NMF),real(f_EM_L1R_WNE_NMF)); 
fprintf('\n');