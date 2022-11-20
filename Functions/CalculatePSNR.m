function [error] = CalculatePSNR(I1,I2,num_band)

%Author Kinan  ABBAS
% Creation Date: 19 Oct 2022
% Purpose: calculate PSNR in case of Parfor loop
% error=zeros(num_band);
for band=1:num_band
    error(band)=psnr(squeeze(I1(:,:,band)),squeeze(I2(:,:,band)),256);
end

end

