%% Cleaning and loading everthing to the path
close all;
clear;
run_me_first;

%% Setting parameters for the expirement

% Parameters for KPWNMF
num_band=25; %    %number of bands
sz=[100,100];     %image size
smp_scenario=1; %sampling operator, 1=binary, 2=random projections, 3=filter based
num_obs_pxl=1; %  %number of exposures
WNMF_params.WNMF_Offset=sqrt(num_band); % Dimensions of the patch
WNMF_params.rank=1;% rank of the patch
WNMF_params.Iter_max_E_step=25; % Number of iteration for the E step in Expectation Maximization algorithm0
WNMF_params.Iter_max_M_step=1; % Number of iteration for the M step in Expectation Maximization algorithm
WNMF_params.Nesterov_Max_Iterations=1000;
WNMF_params.Nesterov_Min_Iterations=10;
WNMF_params.beta=0; %% Used to control graph regularization and L1,L2 norms with nmf
%  WNMF_params.Step_size=floor(WNMF_params.WNMF_Offset/2)-1; %%uncomment
% WNMF_params.Step_size=1; % Uncomment this line if you need overlaps
% between the patches
WNMF_params.Step_size=WNMF_params.WNMF_Offset;
WNMF_params.Scaling=false; % To apply sum to one constrain on the rows of G (Abandunces Matrix) in Nesterov
WNMF_params.I_WB_Initialization=true ; % To use the output of the WB algorithm as input
unmixing_rank=WNMF_params.rank; % Rank for the entire image
% unmixing=true; % put it to ture if you want to test unmixing after democaising
num_of_experiments=1;% The number of times to repeat the same experiment (The average of the results will be taken at the end)
WNMF_params.Kmeans=true;% True = Kmeans for clustering; false= VCA for selecting the endmemebers
WNMF_params.Kmeans_cut=4;   %null for no cut, 1 for cut on the mean, 2 for cut on the sqrt of the mean, 3 half of the mean, 4 median, 5 the first 100 element

unmixing=true; %Ture, if you want to test unmixing after democaising
WNMF_params.NesterovUnmixing=true; % True to use Nesterov as unmixing method for all the 2-stage approaces; False, nfindr with Least Squares will be used instead
WNMF_params.NesterovScaling=15; % To control sum-to-one constraint
With_Noise=false;
% Parameters for GRMR
GRMR_params.offset=5;
GRMR_params.maxIter=5; %was 20
GRMR_params.sgm2=1e1;
GRMR_params.gamma=0.2; % was 0.2 and modified by Kinan
GRMR_params.rank_sel=2;

% Run the expirement
if(With_Noise==true)
    [mean_PSNR,mean_SAM,exec_time,std_PSNSR,std_SAM,mean_MSE,mean_SIR,mean_MER,mean_U_PSNR,mean_SSIM,mean_RMSE,mean_MRSA]=Evaluate_on_simple_image_noisy(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params,WNMF_params,unmixing,num_of_experiments);
    
else
    [mean_PSNR,mean_SAM,exec_time,std_PSNSR,std_SAM,mean_MSE,mean_SIR,mean_MER,mean_U_PSNR,mean_SSIM,mean_RMSE,mean_MRSA]=Evaluate_on_single_method(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params,WNMF_params,unmixing,num_of_experiments);
end
% Showing the results
fprintf('\n');
fprintf('Method | U_PSNR | D_PSNR | SAM  | MSE       | RMSE | SSIM | SIR  | MER  | Time | MRSA \n');
fprintf('KPWNMF | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean_U_PSNR(1),mean_PSNR(1),mean_SAM(1),mean_MSE(1),mean_RMSE(1),mean_SSIM(1),mean_SIR(1),mean_MER(1),exec_time(6),mean_MRSA(1));
