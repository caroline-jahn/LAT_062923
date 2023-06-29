
task = 'Learning_Attentional_Templates';
fsroot = fullfile(filesep,'Volumes','buschman'); %fullfile(filesep,'Volumes','buschman'), 'Z:' fullfile(filesep,'jukebox','buschman')

SD_list=[1.8 2 2.2 2.5 2.8 3 3.2 3.5 3.8 4 4.5 4.8 5];
arrayID_list=[4,9,13,15,19,21,23,27,29,31,33,35,37,39,41,43,45,47]; %19,,27,47,49 BUGS to fix

for i_sd=1:length(SD_list)
    
    SD = SD_list(i_sd)

anova_sig=0.05;

subtask='delsacc';

    nSD   = SD; %how many sd to use as cutoff?
	sdstr = sprintf('%.3dsd', round(10*nSD));


deststem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,sdstr);
destpath = fullfile(fsroot,deststem);

load(fullfile(destpath,'Stats_restricted_FEF'),'Stats');

Total_LIP(i_sd)=Stats.LIP.delay_on;
Total_FEF(i_sd)=Stats.FEF.delay_on;
Total_PFC(i_sd)=Stats.PFC.delay_on;

end

