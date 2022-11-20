function [Unmixing_Image, Final_G , Final_F]=low_rank_completion_methods(I_MOS_seq,SMP_seq,num_band,WNMF_params,num_obs_pxl,I_WNMF_rec,F_pool,Final_norm)

%Input:
    % I_MOS_seq:  The SSI image
    % SMP_seq: The sampling matrix
    % num_band: The number of wavelenghts
    % WNMF_params: paramters required to run the algorithm:
        % WNMF_params.Kmeans: true to apply Kmeans, fales to apply VCA
        % WNMF_params.Kmeans_cut: to determine the threshold for kmeans/VCA,
        % possible values  %null for no cut, 1 for cut on the mean, 2 for cut on the sqrt of the mean, 3 half of the mean, 4 median, 5 the first 100 element
        % WNMF_params.global_rank: the rank of the full image for clustering.
    % num_obs_pxl: Number of scans of the scene.
    % I_WNMF_rec: The completed images using WNMF method (3d datacube)
    % F_pool: The pool of spectra collected from the rank-1 patches
    % Final_norm: The approximation error for every spectrum in the F_pool


%Output:
    % Final_G: The abundaces matrix
    % Final_F: The endmembers matrix
    % Unmixing_Image: The estimated image after multiplying Final_F*Final_G

% Author: Kinan ABBAS
% Creation Date: OCT 5 2022
% Last update date: OCT 10 2022
    
[n1,n2,n3]=size(I_MOS_seq);
%% Filtering the pool of the spectra
F_pool_filtered=[];
if(isempty(WNMF_params.Kmeans_cut)||WNMF_params.Kmeans_cut==0)
    mean_F_pool=500000000000000000000;
elseif WNMF_params.Kmeans_cut==1
    mean_F_pool=mean(Final_norm);
elseif WNMF_params.Kmeans_cut==2
    mean_F_pool=mean(Final_norm);
    if(mean_F_pool<1)
        mean_F_pool=mean_F_pool*mean_F_pool;
    else
        mean_F_pool=sqrt(mean_F_pool);
    end
elseif WNMF_params.Kmeans_cut==3
    mean_F_pool=mean(Final_norm);
    if(mean_F_pool<1)
        mean_F_pool=mean_F_pool*2;
    else
        mean_F_pool=mean_F_pool/2;
    end
elseif WNMF_params.Kmeans_cut==4
    mean_F_pool=median(Final_norm);
    
elseif WNMF_params.Kmeans_cut==5
    Final_norm_temp=sort(Final_norm);
    mean_F_pool=Final_norm_temp(400);
    
end
j=1;
for i=1:size(Final_norm,1)
    if(Final_norm(i)<(mean_F_pool))
        F_pool_filtered(j,:)=F_pool(i,:);
        j=j+1;
    end
end
%% Run the clustring stage to estimate the final endmembers Final_F variable
disp('Clustering');
% KPWNMF
if(WNMF_params.Kmeans==true||isempty(WNMF_params.Kmeans)) % For running K-means
    disp('KPWNMF');
    [~,C]=kmeans(F_pool_filtered,WNMF_params.global_rank);
    Final_F=C;
% VCA_PWNMF    
else
    disp('VCA_PWNMF');
    [~, K] = VCA(F_pool_filtered','Endmembers',WNMF_params.global_rank,'verbose','off');
    Final_F=F_pool_filtered([K],:);
    
end

%% Estimating the abundances
disp('Estimating G');
% Unfolding the 3d datacube
A=reshape(I_WNMF_rec,[n1*n2,num_band]);

% Unfolding the sampling matrix
B=reshape(SMP_seq,[n1*n2,num_band,size(SMP_seq,4)]);

% Process the multiple scans case for the sampling matrix
smp_temp=zeros(n1*n2,n3,num_obs_pxl);
for ii=1:num_obs_pxl
    smp_temp=smp_temp(:,:,1)+B(:,:,ii);
end
B=squeeze(smp_temp(:,:,1));

% Estimating G using NMF
%Sum to one constraint
% A1=[A,WNMF_params.NesterovScaling*ones(size(A,1),1)]; Final_F1=[Final_F,WNMF_params.NesterovScaling*ones(WNMF_params.global_rank,1)];
% [ ~ , Final_G , iter,elapse,HIS]=NeNMF_Fixed_W(A1',WNMF_params.global_rank,'MAX_ITER',1000,'MIN_ITER',10,'W_FFIXID',Final_F1','H_INIT',WNMF_params.Ginit','SCALLING',false);
% Final_G=Final_G' ;

% Estimating G using matlab function, Hyperspectral Imaging Toolbox must be installed first 
abundanceMap = estimateAbundanceLS(I_WNMF_rec,Final_F','Method','ncls');
Final_G=reshape(abundanceMap,[n1*n2,WNMF_params.global_rank]);

oness=ones(size(A,1),size(A,2));
Unmixing_Image=B.*A + (oness-B).*(Final_G*Final_F);

end




