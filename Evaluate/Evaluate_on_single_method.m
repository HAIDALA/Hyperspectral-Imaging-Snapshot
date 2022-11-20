function [mean_PSNR,mean_SAM, exec_time,std_PSNSR,std_SAM,mean_MSE,mean_SIR,mean_MER,mean_U_PSNR,mean_SSIM,mean_RMSE,mean_MRSA]=Evaluate_on_single_method(num_band,sz,smp_scenario,num_obs_pxl,GRMR_params,WNMF_params,unmixing,num_of_experiments)


% Test demosaicing and unmixing on complex Image
% Author: Kinan ABBAS
% Creation Date: Oct 14 2022
% Last modification date: Oct 14 2022


show_Figuers=true;
for zz=1:num_of_experiments
    %% Loading the cube
    load('Water_Metal_Concret.mat')
%     load('abundance_simple_image.mat');
    load('abundance_complex_image.mat');
    r=3% The rank of the image is fixed 
    Selected_Wavelengths=[12,50,93,147,180,201,263,304,317,329,395,464,569,600,632,650,662,665,700,710,718,783,799,840,845];
    if(num_band==25)
        F_HS=F_HS(1:r,Selected_Wavelengths);

    else
        F_HS=F_HS(1:r,Selected_Wavelengths(1:16));
    end
    I_HS_Temp=G_HS*F_HS;
    I_HS=reshape(I_HS_Temp,[sz(1),sz(2),num_band]);
%     I_HS=I_HS*255;
%     I_HS=round(I_HS);
    %% Initializing some values
    [n1,n2,n3]=size(I_HS);
    m = n1*n2; n = num_band;
    Ginit = rand(m,r)*10+0.1;
    Ginit=ScaleRows(Ginit);
    Finit = rand(r,n)+0.0001;
    WNMF_params.Ginit=Ginit;
    WNMF_params.Finit=Finit;
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
    [I_WNMF_rec,F_pool,Final_norm,I_WNMF_err]=WNMF_Demosaicing_With_Analysis(I_MOS_seq,SMP_seq,FilterPattern_Final,num_band,WNMF_params,num_obs_pxl,I_WB,smp_scenario);
    [I_WNMF_Unmixing,G_WNMF,F_WNMF]=low_rank_completion_methods(I_MOS_seq,SMP_seq,num_band,WNMF_params,num_obs_pxl,I_WNMF_rec,F_pool,Final_norm);
    WNMF=I_WNMF_rec;
    WNMF_toc=toc;

    


    
    %% Calculate Metrics
    
    %% PSNR for Demosaicing
    for band=1:num_band

        err_WNMF(zz,band)=psnr(squeeze(WNMF(:,:,band)),squeeze(I_HS(:,:,band)),256);
    end
    
    
    
    %% SAM for the restored spectra
   
    for i=1:r
        for j=1:r
            err_SAM_WNMF_temp(i,j)=sam(F_HS(i,:),F_WNMF(j,:));

        end

        err_SAM_WNMF_t(i)=min(err_SAM_WNMF_temp(i,:));
    end
    
    err_SAM_WNMF(zz)=mean(err_SAM_WNMF_t);
    
    %% Reproduce the images
%     WNMF_Image=G_WNMF*F_WNMF;
%     WNMF_Image=reshape(WNMF_Image,[sz(1),sz(2),num_band]);
    WNMF_Image=I_WNMF_Unmixing;
    WNMF_Image=reshape(WNMF_Image,[sz(1),sz(2),num_band]);
    
    %% PSNR after unmixing
    for band=1:num_band
        err_M_WNMF(zz,band)=psnr(squeeze(WNMF_Image(:,:,band)),squeeze(I_HS(:,:,band)),256);
    end
    %
    
    %% MSE
    
    for band=1:num_band
        MSE_WNMF(zz,band)=MSE(squeeze(WNMF_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        
    end
    
    %% RMSE
    for band=1:num_band
        RMSE_WNMF(zz,band)=RMSE(squeeze(WNMF_Image(:,:,band)),squeeze(I_HS(:,:,band)));
        
    end
    
    %% SSIM
    for band=1:num_band
        SSIM_WNMF(zz,band)=ssim(WNMF_Image(:,:,band),I_HS(:,:,band));
    end
    
    
    %% MRSA (mean -removed spectral angle)
    [Vout, MRSA_WNMF_temp] = performances(F_WNMF, F_HS);   
    MRSA_WNMF(zz)=mean(MRSA_WNMF_temp);
    
    
    %% SIR
    disp('Calculating SIR..');
    %         [G_HS,F_HS,iter_HS,elapse_HS,HIS_HS ] = unmix(I_HS,r,Ginit,Finit);
    [SDR_WNMF,SIR_WNMF_temp,SAR_WNMF,perm_WNMF] = bss_eval_sources(F_WNMF,F_HS);
    SIR_WNMF(zz)=mean(SIR_WNMF_temp);

    
    
    %% MER
    disp('Calculating MER..');
    
    [MER_WNMF(zz)] = mean(bss_eval_mix(G_WNMF,G_HS)); 
    
    %% Final procesing
    % process to show the SSI image
    tmp=reshape(I_MOS_seq,[m,num_obs_pxl]);
    tmp2=reshape(SMP_seq,[m,n,num_obs_pxl]);
    FinalMatrix=applyBWR(tmp2,tmp);
    I_MOS_SSI=reshape(FinalMatrix,[n1,n2,num_band]);
    if(num_band==25)
        load('Data\pure_patches_25_complex_image.mat');
        for i=1:5:sz(1)-4;
            for j=1:5:sz(2)-4;
                tmp=I_MOS_seq(i:1:i+4,j:1:j+4);
                if(isequal(tmp,pure_patch_1)|isequal(tmp,pure_patch_2)|isequal(tmp,pure_patch_3))
                    I_MOS_seq(i:1:i+4,j:1:j+4)=200;
                end
            end
        end
    end
    
    
    if(show_Figuers)
        
        h1=figure;
        h2=figure;
        
        
        
        for tt=1:1

            figure(h1);  imagesc(I_MOS_seq,[0,255]); %colormap('gray');
            figure(h2);  imagesc(I_WNMF_err(:,:,1),[0,255])
            
        end
        
        
        
    end
    exec_time=[0,0,0,0,0,WNMF_toc];
    fprintf('Naive  | %.1f   | %.1f   | %.2f | %.2f      | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f \n',mean(err_M_WNMF(zz,:)),mean(err_WNMF(zz,:)),err_SAM_WNMF(zz),mean(MSE_WNMF(zz,:)),mean(RMSE_WNMF(zz,:)),mean(SSIM_WNMF(zz,:)),SIR_WNMF(zz),MER_WNMF(zz),exec_time(6),MRSA_WNMF(zz));
    %     mean_U_PSNR=[mean(mean(err_M_GRMR)),0,mean(mean(err_M_WB)),mean(mean(err_M_PPID)),mean(mean(err_M_ItSD)),mean(mean(err_M_WNMF)),mean(mean(err_M_WNMF1))];
end

exec_time=[0,0,0,0,0,WNMF_toc];
mean_PSNR=[mean(mean(err_WNMF)),0];
mean_SAM=[mean((err_SAM_WNMF)),0];
std_PSNSR=[std(mean(err_WNMF)),0];
std_SAM=[std((err_SAM_WNMF)),0];
mean_MER=[mean(MER_WNMF),0];
mean_SIR=[mean(SIR_WNMF),0];
mean_U_PSNR=[mean(mean(err_M_WNMF)),0];
mean_SSIM =[mean((SSIM_WNMF)),0];
mean_RMSE =[mean((RMSE_WNMF)),0];
mean_MSE =[mean((MSE_WNMF)),0];
mean_MRSA =[mean(MRSA_WNMF),0];