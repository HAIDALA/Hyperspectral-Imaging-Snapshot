function [mean_PSNR]=evaluate_on_CAVE_Demosaicing(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params,WNMF_params,dataset_size)


% Test demosaicing performance on CAVE dataset
% Author: Kinan ABBAS
% Creation Date: Oct 11 2022


show_Figuers=true;
load('Cave_Dataset_Rank.mat');
fnames=dir('Data/complete_ms_data');

for zz=1:dataset_size
    %% Loading the cube
    zz=6;
    fname=fnames(zz+2).name;
    
    r=GetCaveImageRank(fname,Cave_names,Cave_Rank);
    fprintf('%s, %d out of %d\n',fname,zz,numel(fnames)-2);
%     fprintf('Rank is %d \n',r);
    %     Load hyperspectral cube
    I_HS=load_hypercube_CAVE(fname,sz);
    I_HS=I_HS(:,:,1:num_band);
    % Normalize the values to be between 0 and 255
    mx=max(max(max(I_HS)));
%     I_HS=I_HS./mx;
%     I_HS=I_HS*255;
%     I_HS=round(I_HS);
    %% Initializing some values
    [n1,n2,n3]=size(I_HS);
    m = n1*n2; n = num_band;
%     Ginit = rand(m,r)*10+0.1;
%     Ginit=ScaleRows(Ginit);
%     Finit = rand(r,n)+0.0001;
%     WNMF_params.Ginit=Ginit;
%     WNMF_params.Finit=Finit;    
    
    %% Applying the mosaic filter to aquire SSI image out of the full cube
    if num_band==25
        load('spectral_responses_5x5.mat');
    elseif num_band==16
        load('spectral_responses_4x4.mat');
        CentralWavelengths=CentralWavelength;
    else
        disp('Error');
    end
    
    temp2=sort( round(CentralWavelengths))-400;
    SpectralProfiles=SpectralProfiles(:,temp2);
    SpectralProfiles=rot90(SpectralProfiles);
    
    [n1,n2,n3]=size(I_HS);
    [SMP_seq,FilterPattern_lst]=make_sampling_operators2(n1,n2,n3,num_obs_pxl,num_band,smp_scenario,SpectralProfiles);
    [I_MOS_seq]=acquire_observations(I_HS,SMP_seq,num_obs_pxl);% The SSI Image
    SMP_SEQ=SMP_seq; % The Sampling Matrix
    
    
    %% Initialize some variables
    I_WB_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
    I_BTES_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
    I_ItSD_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
    I_PPID_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
    
    
    
    
    
    %% Weighted bilinear interpolation (WB)
    disp('Running WB');
    tic;
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        I_WB_tmp(:,:,:,pp)=run_WB(I_MOS,FilterPattern,num_band);
    end
    I_WB=mean(I_WB_tmp,4);
    WB_toc=toc;
    %% Graph and Rank Regularized low rank matrix approximation (GRMR)
    sgm2=GRMR_params.sgm2;
    offset=GRMR_params.offset;
    maxIter=GRMR_params.maxIter;
    sgm2=GRMR_params.sgm2;
    gamma=GRMR_params.gamma;
    rank_sel=GRMR_params.rank_sel;
    disp('Running GRMR');
    tic;
    I_GRMR_rec=run_GRMR_demosaick(I_MOS_seq,SMP_SEQ,num_band,offset,sgm2,maxIter,rank_sel,gamma,I_WB);
    GRMR_toc=toc;
    
    %%    Multispectral Demosaicing using Pseudo-panchromatic image (PPID)
    
    disp('Running PPID');
    tic
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        PPI=mean(squeeze(I_WB_tmp(:,:,:,pp)),3);
        I_PPID_tmp(:,:,:,pp)=run_PPID(I_MOS,FilterPattern,num_band,PPI);
    end
    I_PPID=mean(I_PPID_tmp,4);
    PPID_toc=toc;
    
    %%    KPWNMF
    FilterPattern_Final=zeros(n1,n2,num_obs_pxl);
    % convert FilterPattern_1st to a multidimension array
    for gg=1:num_obs_pxl
        cc=cell2mat(FilterPattern_lst(gg));
        FilterPattern_Final(:,:,gg)=cc;
    end
    disp('Running WNMF with K-means');
    tic;
    [I_WNMF_rec,~,~]=WNMF_Demosaicing(I_MOS_seq,SMP_seq,FilterPattern_Final,num_band,WNMF_params,num_obs_pxl,I_WB,smp_scenario);
    WNMF=I_WNMF_rec;
    WNMF_toc=toc;
    
    %%    VCA_PWNMF
%     disp('Running WNMF with VCA');
%     WNMF_params.Kmeans=false;
%     tic;
%     disp('Running Naive WNMF');
%     WNMF_params1=WNMF_params;
%     WNMF_params1.WNMF_Offset=sz(1); % 
%     WNMF_params1.rank=r;% rank of the patch
%     WNMF_params1.Iter_max_E_step=25; % Number of iteration for the E step in Expectation Maximization algorithm
%     WNMF_params1.Iter_max_M_step=1; % Number of iteration for the M step in Expectation Maximization algorithm
%     WNMF_params1.Nesterov_Max_Iterations=1000;
%     WNMF_params1.Nesterov_Min_Iterations=10;
%     WNMF_params1.global_rank=r;
%     WNMF_params1.Kmeans=false;
%     WNMF_params1.Scaling=false;
%     tic;
%     [I_WNMF_rec1,~,~]=Naive_WNMF(I_MOS_seq,SMP_seq,FilterPattern_Final,num_band,WNMF_params1,num_obs_pxl,I_WB,smp_scenario);  
%     WNMF1=I_WNMF_rec1;
%     WNMF_toc1=toc;
    
    %%  Binary Tree-Based Edge-Sensing (BTES)
    
    disp('Running BTES');
    tic;
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        I_BTES_tmp(:,:,:,pp)=run_BTES(I_MOS,FilterPattern,num_band,squeeze(I_WB_tmp(:,:,:,pp)));
    end
    I_BTES=mean(I_BTES_tmp,4);
    BTES_toc=WB_toc;
    
    %% Iterative Spectral Difference (ItSD)
    disp('Running ItSD');
    tic;
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        I_ItSD_tmp(:,:,:,pp)=run_ItSD(I_MOS,FilterPattern,num_band);
    end
    I_ItSD=mean(I_ItSD_tmp,4);

    ItSD_toc=toc;
    %
    
    
    %% PSNR for Demosaicing
    for band=1:num_band
        err_GRMR(zz,band)=psnr(squeeze(I_GRMR_rec(:,:,band)),squeeze(I_HS(:,:,band)));
        err_PPID(zz,band)=psnr(squeeze(I_PPID(:,:,band)),squeeze(I_HS(:,:,band)));
        err_WB(zz,band)=psnr(squeeze(I_WB(:,:,band)),squeeze(I_HS(:,:,band)));
        err_ItSD(zz,band)=psnr(squeeze(I_ItSD(:,:,band)),squeeze(I_HS(:,:,band)));
        err_BTES(zz,band)=psnr(squeeze(I_BTES(:,:,band)),squeeze(I_HS(:,:,band)));
        err_WNMF(zz,band)=psnr(squeeze(WNMF(:,:,band)),squeeze(I_HS(:,:,band)));
%         err_WNMF1(zz,band)=psnr(squeeze(WNMF1(:,:,band)),squeeze(I_HS(:,:,band)),256);
    end
       
    %% Final procesing
    % process to show the SSI image
    tmp=reshape(I_MOS_seq,[m,num_obs_pxl]);
    tmp2=reshape(SMP_seq,[m,n,num_obs_pxl]);
    FinalMatrix=applyBWR(tmp2,tmp);
    I_MOS_SSI=reshape(FinalMatrix,[n1,n2,num_band]);
    
    if(show_Figuers)
        
        h1=figure;
        h2=figure;
        h3=figure;
        h4=figure;
        h5=figure;
        h6=figure;
        h7=figure;
        h8=figure;
        
        
        for tt=1:1
            figure(h1);  imagesc(squeeze(I_HS(:,:,tt))); title('TRUE'); %colormap('gray');
            figure(h2);  imagesc(squeeze(I_MOS_SSI(:,:,tt))); title('Starting Image'); %colormap('gray');
            figure(h3);  imagesc(squeeze(I_GRMR_rec(:,:,tt))); title('GRMR'); %colormap('gray');
            figure(h4);  imagesc(squeeze(I_BTES(:,:,tt))); title('BTES'); %colormap('gray');
            figure(h5);  imagesc(squeeze(I_WB(:,:,tt))); title('WB'); %colormap('gray');
            figure(h6);  imagesc(squeeze(I_PPID(:,:,tt))); title('PPID'); %colormap('gray');
            figure(h7);  imagesc(squeeze(I_ItSD(:,:,tt))); title('ItSD'); %colormap('gray');
            figure(h8);  imagesc(squeeze(WNMF(:,:,tt))); title('WNMF'); %colormap('gray');
            pause(1);
            
        end
        
        
        
        
    end
%     exec_time=[GRMR_toc,BTES_toc,WB_toc,PPID_toc,ItSD_toc,WNMF_toc,WNMF_toc1];
    fprintf('Method  D_PSNR  \n');
    fprintf('GMRM   | %.1f   \n',mean(err_GRMR(zz,:)));
    fprintf('BTES   | %.1f   \n',mean(err_BTES(zz,:)));
    fprintf('WB     | %.1f   \n',mean(err_WB(zz,:)));
    fprintf('PPID   | %.1f   \n',mean(err_PPID(zz,:)));
    fprintf('ItSD   | %.1f   \n',mean(err_ItSD(zz,:)));
    fprintf('KPWNMF | %.1f   \n',mean(err_WNMF(zz,:)));
%     fprintf('Naive  | %.1f   \n',mean(err_WNMF1(zz,:)));
    
    %     mean_U_PSNR=[mean(mean(err_M_GRMR)),0,mean(mean(err_M_WB)),mean(mean(err_M_PPID)),mean(mean(err_M_ItSD)),mean(mean(err_M_WNMF)),mean(mean(err_M_WNMF1))];
end

% mean_PSNR=[mean(mean(err_GRMR)),mean(mean(err_BTES)),mean(mean(err_WB)),mean(mean(err_PPID)),mean(mean(err_ItSD)),mean(mean(err_WNMF)),mean(mean(err_WNMF1))];
mean_PSNR=[mean(mean(err_GRMR)),mean(mean(err_BTES)),mean(mean(err_WB)),mean(mean(err_PPID)),mean(mean(err_ItSD)),mean(mean(err_WNMF))];

% mean_SAM=[mean((err_SAM_GRMR)),err_SAM_GRMR,mean((err_SAM_WB)),mean((err_SAM_PPID)),mean((err_SAM_ItSD)),mean((err_SAM_WNMF)),mean((err_SAM_WNMF1))];
% std_PSNSR=[std(mean(err_GRMR)),0,std(mean(err_WB)),std(mean(err_PPID)),std(mean(err_ItSD)),std(mean(err_WNMF))];
% std_SAM=[std((err_SAM_GRMR)),0,std((err_SAM_WB)),std((err_SAM_PPID)),std((err_SAM_ItSD)),std((err_SAM_WNMF))];
% mean_MER=[mean(MER_GRMR),mean(MER_BTES),mean(MER_I_WB),mean(MER_PPID),mean(MER_ItSD),mean(MER_WNMF),mean(MER_WNMF1)];
% mean_SIR=[mean(SIR_GRMR),mean(SIR_BTES),mean(SIR_I_WB),mean(SIR_PPID),mean(SIR_ItSD),mean(SIR_WNMF),mean(SIR_WNMF1)];
% mean_U_PSNR=[mean(mean(err_M_GRMR)),mean(mean(err_M_BTES)),mean(mean(err_M_WB)),mean(mean(err_M_PPID)),mean(mean(err_M_ItSD)),mean(mean(err_M_WNMF)),mean(mean(err_M_WNMF1))];
% mean_SSIM =[mean((SSIM_GRMR)),mean((SSIM_BTES)),mean((SSIM_WB)),mean((SSIM_PPID)),mean((SSIM_ItSD)),mean((SSIM_WNMF)),mean((SSIM_WNMF1))];
% mean_RMSE =[mean((RMSE_GRMR)),mean((RMSE_BTES)),mean((RMSE_WB)),mean((RMSE_PPID)),mean((RMSE_ItSD)),mean((RMSE_WNMF)),mean((RMSE_WNMF1))];
% mean_MSE =[mean((MSE_GRMR)),mean(mean(MSE_BTES)),mean((MSE_WB)),mean((MSE_PPID)),mean((MSE_ItSD)),mean((MSE_WNMF)),mean((MSE_WNMF1))];
% mean_MRSA =[mean(MRSA_GRMR),mean(MRSA_BTES),mean(MRSA_WB),mean(MRSA_PPID),mean(MRSA_ItSD),mean(MRSA_WNMF),mean(MRSA_WNMF1)];