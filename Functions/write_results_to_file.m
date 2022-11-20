function write_results_to_file(metric,num_bands,tables)
% Author: Kinan ABBAS
% Creation date: Oct 11 2022
tables.SIR.Properties.VariableNames=["SNR","GRMR","BTES","WB","PPID","ItSD","KPWNMF","VCA_PWNMF","Naive"];
tables.MER.Properties.VariableNames=["SNR","GRMR","BTES","WB","PPID","ItSD","KPWNMF","VCA_PWNMF","Naive"];
tables.PSNR.Properties.VariableNames=["SNR","GRMR","BTES","WB","PPID","ItSD","KPWNMF","VCA_PWNMF","Naive"];
tables.SAM.Properties.VariableNames=["SNR","GRMR","BTES","WB","PPID","ItSD","KPWNMF","VCA_PWNMF","Naive"];
timestamp=datetime('now','TimeZone','local','Format','d_MMM_y_HH_MM_SS');
filename='Results/'+string(metric)+'_'+string(num_bands)+'_bands_'+string(timestamp)+'.xls';
writetable(tables.SIR,filename,'Sheet','SIR');
writetable(tables.MER,filename,'Sheet','MER');
writetable(tables.SAM,filename,'Sheet','SAM');
writetable(tables.PSNR,filename,'Sheet','PSNR');
end

