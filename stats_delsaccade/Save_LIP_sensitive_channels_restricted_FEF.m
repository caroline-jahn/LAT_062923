%% Save LIP sensitive channels
% choose nSD and report LIP_both sensitive channels in ns6directory

fsroot='/Volumes/buschman';
task='Learning_Attentional_Templates';

%load ns6directory
destpath='/Volumes/buschman/Projects/Learning_Attentional_Templates/Data/General';
load(fullfile(destpath,'NS6Directory_AC_restricted_FEF.mat'),'ns6directory');

%load summary stats for nSD chosen
subtask='delsacc';

nSD   = 2.5; %how many sd to use as cutoff? => choose the one with most delay on response?
sdstr = sprintf('%.3dsd', round(10*nSD));


dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,sdstr);
dirpath = fullfile(fsroot,dirstem);

load(fullfile(dirpath,'StimTypeFR_AllChannels_restricted_FEF'),'StimTypeFR_AllChannels');

%save in ns6directory which LIP channels were delay on sensitive
for i=1:length(StimTypeFR_AllChannels)
    
    ns6directory(StimTypeFR_AllChannels(i).ArrayID).LIP_sensitive=StimTypeFR_AllChannels(i).LIP.Delay_on_005_both; %delayed saccade
    ns6directory(StimTypeFR_AllChannels(i).ArrayID-1).LIP_sensitive=StimTypeFR_AllChannels(i).LIP.Delay_on_005_both; % explore exploit => arrayID -1
    if ~isempty(StimTypeFR_AllChannels(i).FEF.ActiveChannels)
        ns6directory(StimTypeFR_AllChannels(i).ArrayID).FEF_sensitive=StimTypeFR_AllChannels(i).FEF.Delay_on_005_both; %delayed saccade
        ns6directory(StimTypeFR_AllChannels(i).ArrayID-1).FEF_sensitive=StimTypeFR_AllChannels(i).FEF.Delay_on_005_both; % explore exploit => arrayID -1
    else
        ns6directory(StimTypeFR_AllChannels(i).ArrayID).FEF_sensitive=[]; %delayed saccade
        ns6directory(StimTypeFR_AllChannels(i).ArrayID-1).FEF_sensitive=[]; % explore exploit => arrayID -1
    end
    ns6directory(StimTypeFR_AllChannels(i).ArrayID).PFC_sensitive=StimTypeFR_AllChannels(i).PFC.Delay_on_005_both; %delayed saccade
    ns6directory(StimTypeFR_AllChannels(i).ArrayID-1).PFC_sensitive=StimTypeFR_AllChannels(i).PFC.Delay_on_005_both; % explore exploit => arrayID -1
    
end

save(fullfile(destpath,'NS6Directory_SC_restricted_FEF.mat'),'ns6directory')