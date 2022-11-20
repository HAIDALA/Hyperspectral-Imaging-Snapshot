function [mean_PSNR,mean_SAM, exec_time,std_PSNSR,std_SAM,mean_MSE,mean_SIR,mean_MER,mean_U_PSNR,mean_SSIM,mean_RMSE,mean_MRSA]=evaluate_on_CAVE_Demosaicing_and_Unmixing(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params,WNMF_params,unmixing,dataset_size)


% Test demosaicing and unmixing on CAVE dataset
% Oct 5 2022 - Kinan ABBAS for ICASSP 2023 conference


load('Cave_Dataset_Rank.mat');
show_Figuers=true;

fnames=dir('Data/complete_ms_data');

for zz=1:dataset_size
    %% Loading the cube
    fname=fnames(zz+2).name;
    r=GetCaveImageRank(fname,Cave_names,Cave_Rank);
    fprintf('%s, %d out of %d\n',fname,zz,numel(fnames)-2);
    fprintf('Rank is %d \n',r);
    %     Load hyperspectral cube
    I_HS=load_hypercube_CAVE(fname,sz);
    
    % Normalize the values to be between 0 and 255
    mx=max(max(max(I_HS)));
    I_HS=I_HS./mx;
    I_HS=I_HS*255;
    I_HS=round(I_HS);
    %% Initializing some values
    [n1,n2,n3]=size(I_HS);
    m = n1*n2; n = num_band;
    Ginit = rand(m,r)*10+0.1;
    Ginit=ScaleRows(Ginit);
    Finit = rand(r,n)+0.0001;
    WNMF_params.Ginit=Ginit;
    WNMF_params.Finit=Finit;
    %% Unmix the full cube I_HS using Nfindr method, Hyperspectral imaging toolbox must be installled
    F_HS=ppi(I_HS,r)';
    %     abundanceMap = estimateAbundanceLS(I_HS,F_HS','Method','ncls');
    %     G_HS=reshape(abundanceMap,[n1*n2,r]);
    A=reshape(I_HS,[n1*n2,n3]);
    A1=[A,WNMF_params.NesterovScaling*ones(size(A,1),1)]; Final_F1=[F_HS,WNMF_params.NesterovScaling*ones(r,1)];
    [ ~ , G_HS , iter,elapse,HIS]=NeNMF_Fixed_W(A1',r,'MAX_ITER',1000,'MIN_ITER',10,'W_FFIXID',Final_F1','H_INIT',WNMF_params.Ginit','SCALLING',false);
    G_HS=G_HS';
    F_HS=F_HS(:,1:num_band);
    I_HS=I_HS(:,:,1:num_band);
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
    if(unmixing)
        if(WNMF_params.NesterovUnmixing==true)
            [G_WB,F_WB,iter_WB,elapse_WB,HIS_WB, ] = unmix(I_WB,r,Ginit,Finit,WNMF_params.NesterovScaling);
        else
            F_WB=ppi(I_WB,r)';
            abundanceMap = estimateAbundanceLS(I_WB,F_WB','Method','ncls');
            G_WB=reshape(abundanceMap,[n1*n2,r]);
        end
    end
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
    if(unmixing)
        if(WNMF_params.NesterovUnmixing==true)
            [G_GRMR,F_GRMR,iter_GRMR,elapse_GRMR,HIS_GRMR ] = unmix(I_GRMR_rec,r,Ginit,Finit,WNMF_params.NesterovScaling);
        else
            F_GRMR=ppi(I_GRMR_rec,r)';
            abundanceMap = estimateAbundanceLS(I_GRMR_rec,F_GRMR','Method','ncls');
            G_GRMR=reshape(abundanceMap,[n1*n2,r]);
            
        end
    end
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
    if(unmixing)
        if(WNMF_params.NesterovUnmixing==true)
            [G_PPID,F_PPID,iter_PPID,elapse_PPID,HIS_PPID ] = unmix(I_PPID,r,Ginit,Finit,WNMF_params.NesterovScaling);
        else
            F_PPID=ppi(I_PPID,r)';
            abundanceMap = estimateAbundanceLS(I_PPID,F_PPID','Method','ncls');
            G_PPID=reshape(abundanceMap,[n1*n2,r]);
        end
    end
    PPID_toc=toc;
    
    %%    KPWNMF
    FilterPattern_Final=zeros(n1,n2,num_obs_pxl);
    % convert FilterPattern_1st to a multidimension array
    WNMF_params.global_rank=r;
    for gg=1:num_obs_pxl
        cc=cell2mat(FilterPattern_lst(gg));
        FilterPattern_Final(:,:,gg)=cc;
    end
    disp('Running WNMF with K-means');
    tic;
    [I_WNMF_rec,F_pool,Final_norm]=WNMF_Demosaicing(I_MOS_seq,SMP_seq,FilterPattern_Final,num_band,WNMF_params,num_obs_pxl,I_WB,smp_scenario);
    [I_WNMF_Unmixing,G_WNMF,F_WNMF]=low_rank_completion_methods(I_MOS_seq,SMP_seq,num_band,WNMF_params,num_obs_pxl,I_WNMF_rec,F_pool,Final_norm);
    WNMF=I_WNMF_rec;
    WNMF_toc=toc;
    
    %%    VCA_PWNMF
    disp('Running WNMF with VCA');
    WNMF_params.Kmeans=false;
    tic;
    [I_WNMF_Unmixing1,G_WNMF1,F_WNMF1]=low_rank_completion_methods(I_MOS_seq,SMP_seq,num_band,WNMF_params,num_obs_pxl,I_WNMF_rec,F_pool,Final_norm);    
    WNMF1=I_WNMF_rec;
    WNMF_toc1=toc;
    
    %%  Binary Tree-Based Edge-Sensing (BTES)
    
    disp('Running BTES');
    tic;
    for pp=1:num_obs_pxl
        I_MOS=I_MOS_seq(:,:,pp);
        FilterPattern=cell2mat(FilterPattern_lst(pp));
        I_BTES_tmp(:,:,:,pp)=run_BTES(I_MOS,FilterPattern,num_band,squeeze(I_WB_tmp(:,:,:,pp)));
    end
    I_BTES=mean(I_BTES_tmp,4);
    %I_BTES=I_WB;
    if(unmixing)
        if(WNMF_params.NesterovUnmixing==true)
            [G_BTES,F_BTES,iter_BTES,elapse_BTES,HIS_BTES ] = unmix(I_BTES,r,Ginit,Finit,WNMF_params.NesterovScaling);
        else
            F_BTES=ppi(I_BTES,r)';
            abundanceMap = estimateAbundanceLS(I_BTES,F_BTES','Method','ncls');
            G_BTES=reshape(abundanceMap,[n1*n2,r]);
        end
    end
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
    if(unmixing)
        if(WNMF_params.NesterovUnmixing==true)
            [G_ItSD,F_ItSD,iter_ItSD,elapse_ItSD,HIS_ItSD ] = unmix(I_ItSD,r,Ginit,Finit,WNMF_params.NesterovScaling);
        else
            F_ItSD=ppi(I_ItSD,r)';
            abundanceMap = estimateAbundanceLS(I_ItSD,F_ItSD','Method','ncls');
            G_ItSD=reshape(abundanceMap,[n1*n2,r]);
        end
    end
    ItSD_toc=toc;
    %
    
    %% Calculate Metrics
    
    %% PSNR for Demosaicing
    for band=1:num_band
        err_GRMR(zz,band)=psnr(squeeze(I_GRMR_rec(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_PPID(zz,band)=psnr(squeeze(I_PPID(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_WB(zz,band)=psnr(squeeze(I_WB(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_ItSD(zz,band)=psnr(squeeze(I_ItSD(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_BTES(zz,band)=psnr(squeeze(I_BTES(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_WNMF(zz,band)=psnr(squeeze(WNMF(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_WNMF1(zz,band)=psnr(squeeze(WNMF1(:,:,band)),squeeze(I_HS(:,:,band)),256);
    end
    
    
    
    %% SAM for Demosaicing
    for xx=1:size(I_GRMR_rec,1)
        for yy=1:size(I_GRMR_rec,2)
            tmp0=squeeze(I_HS(xx,yy,:))+eps;
            
            tmp1=round(squeeze(I_GRMR_rec(xx,yy,:)))+eps;
            err_SAM_GRMR_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
            
            tmp1=round(squeeze(I_PPID(xx,yy,:)))+eps;
            err_SAM_PPID_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
            
            tmp1=round(squeeze(I_WB(xx,yy,:)))+eps;
            err_SAM_WB_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
            
            tmp1=round(squeeze(I_ItSD(xx,yy,:)))+eps;
            err_SAM_ItSD_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
            
            tmp1=round(squeeze(I_BTES(xx,yy,:)))+eps;
            err_SAM_BTES_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
            
            tmp1=round(squeeze(WNMF(xx,yy,:)))+eps;
            err_SAM_WNMF_temp(xx,yy)=real(hyperSam(tmp1,tmp0));
            
            tmp1=round(squeeze(WNMF1(xx,yy,:)))+eps;
            err_SAM_WNMF_temp1(xx,yy)=real(hyperSam(tmp1,tmp0));
            
        end
    end
    err_SAM_GRMR(zz)=mean(mean(err_SAM_GRMR_temp));
    err_SAM_PPID(zz)=mean(mean(err_SAM_PPID_temp));
    err_SAM_WB(zz)=mean(mean(err_SAM_WB_temp));
    err_SAM_ItSD(zz)=mean(mean(err_SAM_ItSD_temp));
    err_SAM_BTES(zz)=mean(mean(err_SAM_BTES_temp));
    err_SAM_WNMF(zz)=mean(mean(err_SAM_WNMF_temp));
    err_SAM_WNMF1(zz)=mean(mean(err_SAM_WNMF_temp1));
    
    
    %% Reproduce the images
    WNMF_Image=G_WNMF*F_WNMF;
    WNMF_Image=reshape(WNMF_Image,[sz(1),sz(2),num_band]);
    
    WNMF_Image1=G_WNMF1*F_WNMF1;
    WNMF_Image1=reshape(WNMF_Image1,[sz(1),sz(2),num_band]);
    
    PPID_Image=G_PPID*F_PPID;
    PPID_Image=reshape(PPID_Image,[sz(1),sz(2),num_band]);
    
    ItSD_Image=G_ItSD*F_ItSD;
    ItSD_Image=reshape(ItSD_Image,[sz(1),sz(2),num_band]);
    
    GRMR_Image=G_GRMR*F_GRMR;
    GRMR_Image=reshape(GRMR_Image,[sz(1),sz(2),num_band]);
    
    WB_Image=G_WB*F_WB;
    WB_Image=reshape(WB_Image,[sz(1),sz(2),num_band]);
    
    BTES_Image=G_BTES*F_BTES;
    BTES_Image=reshape(BTES_Image,[sz(1),sz(2),num_band]);
    %% PSNR after unmixing
    for band=1:num_band
        err_M_WNMF(zz,band)=psnr(squeeze(WNMF_Image(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_M_WNMF1(zz,band)=psnr(squeeze(WNMF_Image1(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_M_PPID(zz,band)=psnr(squeeze(PPID_Image(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_M_GRMR(zz,band)=psnr(squeeze(GRMR_Image(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_M_ItSD(zz,band)=psnr(squeeze(ItSD_Image(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_M_WB(zz,band)=psnr(squeeze(WB_Image(:,:,band)),squeeze(I_HS(:,:,band)),256);
        err_M_BTES(zz,band)=psnr(squeeze(BTES_Image(:,:,band)),squeeze(I_HS(:,:,band)),256);
    end
    %
    
    %% MSE
    
    for band=1:num_band
        MSE_WNMF(zz,band)=MSE(squeeze(WNMF_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        MSE_WNMF1(zz,band)=MSE(squeeze(WNMF_Image1(:,:,band)),squeeze(I_HS(:,:,band)));
        MSE_PPID(zz,band)=MSE(squeeze(PPID_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        MSE_GRMR(zz,band)=MSE(squeeze(GRMR_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        MSE_ItSD(zz,band)=MSE(squeeze(ItSD_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        MSE_WB(zz,band)=MSE(squeeze(WB_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        MSE_BTES(zz,band)=MSE(squeeze(BTES_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        
    end
    
    %% RMSE
    for band=1:num_band
        RMSE_WNMF(zz,band)=RMSE(squeeze(WNMF_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        RMSE_WNMF1(zz,band)=RMSE(squeeze(WNMF_Image1(:,:,band)),squeeze(I_HS(:,:,band)));
        RMSE_PPID(zz,band)=RMSE(squeeze(PPID_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        RMSE_GRMR(zz,band)=RMSE(squeeze(GRMR_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        RMSE_ItSD(zz,band)=RMSE(squeeze(ItSD_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        RMSE_WB(zz,band)=RMSE(squeeze(WB_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        RMSE_BTES(zz,band)=RMSE(squeeze(BTES_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        
    end
    
    %% SSIM
    for band=1:num_band
        SSIM_WNMF(zz,band)=ssim(WNMF_Image(:,:,band),I_HS(:,:,band));
        SSIM_WNMF1(zz,band)=ssim(WNMF_Image1(:,:,band),I_HS(:,:,band));
        SSIM_PPID(zz,band)=ssim(PPID_Image(:,:,band),I_HS(:,:,band));
        SSIM_GRMR(zz,band)=ssim(GRMR_Image(:,:,band),I_HS(:,:,band));
        SSIM_ItSD(zz,band)=ssim(ItSD_Image(:,:,band),I_HS(:,:,band));
        SSIM_WB(zz,band)=ssim(WB_Image(:,:,band),I_HS(:,:,band));
        SSIM_BTES(zz,band)=ssim(BTES_Image(:,:,band),I_HS(:,:,band));
    end
    
    
    %% MRSA (mean -removed spectral angle)
    [Vout, MRSA_WNMF1_temp] = performances(F_WNMF1, F_HS);
    [Vout, MRSA_WNMF_temp] = performances(F_WNMF, F_HS);
    [Vout, MRSA_GRMR_temp] = performances(F_GRMR, F_HS);
    [Vout, MRSA_WB_temp] = performances(F_WB, F_HS);
    [Vout, MRSA_ItSD_temp] = performances(F_ItSD, F_HS);
    [Vout, MRSA_PPID_temp] = performances(F_PPID, F_HS);
    [Vout, MRSA_BTES_temp] = performances(F_BTES, F_HS);
    
    MRSA_WNMF1(zz)=mean(MRSA_WNMF1_temp);
    MRSA_WNMF(zz)=mean(MRSA_WNMF_temp);
    MRSA_GRMR(zz)=mean(MRSA_GRMR_temp);
    MRSA_WB(zz)=mean(MRSA_WB_temp);
    MRSA_ItSD(zz)=mean(MRSA_ItSD_temp);
    MRSA_PPID(zz)=mean(MRSA_PPID_temp);
    MRSA_BTES(zz)=mean(MRSA_BTES_temp);
    
    
    %% SIR
    disp('Calculating SIR..');
    %         [G_HS,F_HS,iter_HS,elapse_HS,HIS_HS ] = unmix(I_HS,r,Ginit,Finit);
    [SDR_WNMF,SIR_WNMF_temp,SAR_WNMF,perm_WNMF] = bss_eval_sources(F_WNMF,F_HS);
    [SDR_WNMF1,SIR_WNMF1_temp,SAR_WNMF1,perm_WNMF1] = bss_eval_sources(F_WNMF1,F_HS);
    [SDR_I_WB,SIR_I_WB_temp,SAR_I_WB,perm_I_WB] = bss_eval_sources(F_WB,F_HS);
    [SDR_GRMR,SIR_GRMR_temp,SAR_GRMR,perm_GRMR] = bss_eval_sources(F_GRMR,F_HS);
    [SDR_PPID,SIR_PPID_temp,SAR_PPID,perm_PPID] = bss_eval_sources(F_PPID,F_HS);
    [SDR_BTES,SIR_BTES_temp,SAR_BTES,perm_BTES] = bss_eval_sources(F_BTES,F_HS);
    [SDR_ItSD,SIR_ItSD_temp,SAR_ItSD,perm_ItSD] = bss_eval_sources(F_ItSD,F_HS);
    SIR_WNMF(zz)=mean(SIR_WNMF_temp);
    SIR_WNMF1(zz)=mean(SIR_WNMF1_temp);
    SIR_I_WB(zz)=mean(SIR_I_WB_temp);
    SIR_PPID(zz)=mean(SIR_PPID_temp);
    SIR_GRMR(zz)=mean(SIR_GRMR_temp);
    SIR_BTES(zz)=mean(SIR_BTES_temp);
    SIR_ItSD(zz)=mean(SIR_ItSD_temp);
    
    
    %% MER
    disp('Calculating MER..');
    
    [MER_WNMF(zz)] = mean(bss_eval_mix(G_WNMF,G_HS));
    [MER_WNMF1(zz)] = mean(bss_eval_mix(G_WNMF1,G_HS));
    [MER_I_WB(zz)] = mean(bss_eval_mix(G_WB,G_HS));
    [MER_GRMR(zz)] = mean(bss_eval_mix(G_GRMR,G_HS));
    [MER_PPID(zz)] = mean(bss_eval_mix(G_PPID,G_HS));
    [MER_ItSD(zz)] = mean(bss_eval_mix(G_ItSD,G_HS));
    [MER_BTES(zz)] = mean(bss_eval_mix(G_BTES,G_HS));
    
    
    
    
    
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
    exec_time=[GRMR_toc,BTES_toc,WB_toc,PPID_toc,ItSD_toc,WNMF_toc,WNMF_toc1];
    fprintf('Method | U_PSNR | D_PSNR | SAM  | MSE       | RMSE | SSIM | SIR  | MER  | Time | MRSA \n');
    fprintf('GMRM   | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean(err_M_GRMR(zz,:)),mean(err_GRMR(zz,:)),err_SAM_GRMR(zz),mean(MSE_GRMR(zz,:)),mean(RMSE_GRMR(zz,:)),mean(SSIM_GRMR(zz,:)),SIR_GRMR(zz),MER_GRMR(zz),exec_time(1),MRSA_GRMR(zz));
    fprintf('BTES   | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean(err_M_BTES(zz,:)),mean(err_BTES(zz,:)),err_SAM_BTES(zz),mean(MSE_BTES(zz,:)),mean(RMSE_BTES(zz,:)),mean(SSIM_BTES(zz,:)),SIR_BTES(zz),MER_BTES(zz),exec_time(2),MRSA_BTES(zz));
    fprintf('WB     | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean(err_M_WB(zz,:)), mean(err_WB(zz,:)), err_SAM_WB(zz),mean(MSE_WB(zz,:)),mean(RMSE_WB(zz,:)),        mean(SSIM_WB(zz,:)),SIR_I_WB(zz),MER_I_WB(zz),exec_time(3),MRSA_WB(zz));
    fprintf('PPID   | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean(err_M_PPID(zz,:)),mean(err_PPID(zz,:)),err_SAM_PPID(zz),mean(MSE_PPID(zz,:)),mean(RMSE_PPID(zz,:)),mean(SSIM_PPID(zz,:)),SIR_PPID(zz),MER_PPID(zz),exec_time(4),MRSA_PPID(zz));
    fprintf('ItSD   | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean(err_M_ItSD(zz,:)),mean(err_ItSD(zz,:)),err_SAM_ItSD(zz),mean(MSE_ItSD(zz,:)),mean(RMSE_ItSD(zz,:)),mean(SSIM_ItSD(zz,:)),SIR_ItSD(zz),MER_ItSD(zz),exec_time(5),MRSA_ItSD(zz));
    fprintf('KPWNMF | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean(err_M_WNMF(zz,:)),mean(err_WNMF(zz,:)),err_SAM_WNMF(zz),mean(MSE_WNMF(zz,:)),mean(RMSE_WNMF(zz,:)),mean(SSIM_WNMF(zz,:)),SIR_WNMF(zz),MER_WNMF(zz),exec_time(6),MRSA_WNMF(zz));
    fprintf('PWNMF  | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean(err_M_WNMF1(zz,:)),mean(err_WNMF1(zz,:)),err_SAM_WNMF1(zz),mean(MSE_WNMF1(zz,:)),mean(RMSE_WNMF1(zz,:)),mean(SSIM_WNMF1(zz,:)),SIR_WNMF1(zz),MER_WNMF1(zz),exec_time(7),MRSA_WNMF1(zz));
    
    %     mean_U_PSNR=[mean(mean(err_M_GRMR)),0,mean(mean(err_M_WB)),mean(mean(err_M_PPID)),mean(mean(err_M_ItSD)),mean(mean(err_M_WNMF)),mean(mean(err_M_WNMF1))];
end

exec_time=[GRMR_toc,BTES_toc,WB_toc,PPID_toc,ItSD_toc,WNMF_toc,WNMF_toc1];
mean_PSNR=[mean(mean(err_GRMR)),mean(mean(err_BTES)),mean(mean(err_WB)),mean(mean(err_PPID)),mean(mean(err_ItSD)),mean(mean(err_WNMF)),mean(mean(err_WNMF1))];
mean_SAM=[mean((err_SAM_GRMR)),err_SAM_GRMR,mean((err_SAM_WB)),mean((err_SAM_PPID)),mean((err_SAM_ItSD)),mean((err_SAM_WNMF)),mean((err_SAM_WNMF1))];
std_PSNSR=[std(mean(err_GRMR)),0,std(mean(err_WB)),std(mean(err_PPID)),std(mean(err_ItSD)),std(mean(err_WNMF))];
std_SAM=[std((err_SAM_GRMR)),0,std((err_SAM_WB)),std((err_SAM_PPID)),std((err_SAM_ItSD)),std((err_SAM_WNMF))];
mean_MER=[mean(MER_GRMR),mean(MER_BTES),mean(MER_I_WB),mean(MER_PPID),mean(MER_ItSD),mean(MER_WNMF),mean(MER_WNMF1)];
mean_SIR=[mean(SIR_GRMR),mean(SIR_BTES),mean(SIR_I_WB),mean(SIR_PPID),mean(SIR_ItSD),mean(SIR_WNMF),mean(SIR_WNMF1)];
mean_U_PSNR=[mean(mean(err_M_GRMR)),mean(mean(err_M_BTES)),mean(mean(err_M_WB)),mean(mean(err_M_PPID)),mean(mean(err_M_ItSD)),mean(mean(err_M_WNMF)),mean(mean(err_M_WNMF1))];
mean_SSIM =[mean((SSIM_GRMR)),mean((SSIM_BTES)),mean((SSIM_WB)),mean((SSIM_PPID)),mean((SSIM_ItSD)),mean((SSIM_WNMF)),mean((SSIM_WNMF1))];
mean_RMSE =[mean((RMSE_GRMR)),mean((RMSE_BTES)),mean((RMSE_WB)),mean((RMSE_PPID)),mean((RMSE_ItSD)),mean((RMSE_WNMF)),mean((RMSE_WNMF1))];
mean_MSE =[mean((MSE_GRMR)),mean(mean(MSE_BTES)),mean((MSE_WB)),mean((MSE_PPID)),mean((MSE_ItSD)),mean((MSE_WNMF)),mean((MSE_WNMF1))];
mean_MRSA =[mean(MRSA_GRMR),mean(MRSA_BTES),mean(MRSA_WB),mean(MRSA_PPID),mean(MRSA_ItSD),mean(MRSA_WNMF),mean(MRSA_WNMF1)];