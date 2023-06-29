%function [StimTypeFR_AllChannels, Stats] = compute_stats(fsroot,task,arrayID_list,SD)

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
    mkdir(destpath);
    
    %load ns6directory which contains the list of active channels for each ROI
    dirstem = fullfile('Projects',task,'Data','General');
    dirpath = fullfile(fsroot,dirstem);
    
    load(fullfile(dirpath,'NS6Directory_AC_restricted_FEF.mat'),'ns6directory');
    
    %combine all StimTypeFR in one structure
    StimTypeFR_AllChannels=tempCombineSegments_restricted_FEF(fsroot,task,subtask,arrayID_list,SD);
    
    %compute the anova
    %for each session
    for i=1:length(arrayID_list)
        StimTypeFR_AllChannels(i).FEF.ActiveChannels=ns6directory(arrayID_list(i)).FEF;
        StimTypeFR_AllChannels(i).PFC.ActiveChannels=ns6directory(arrayID_list(i)).PFC;
        StimTypeFR_AllChannels(i).LIP.ActiveChannels=ns6directory(arrayID_list(i)).LIP;
        
        %for each channel
        for chi=1:length(StimTypeFR_AllChannels(i).FEF.ActiveChannels)
            channel_id=0;
            for j=1:length(StimTypeFR_AllChannels(i).StimTypeFR)
                if StimTypeFR_AllChannels(i).StimTypeFR(j).ID==StimTypeFR_AllChannels(i).FEF.ActiveChannels(chi)
                    channel_id=j;
                end
            end
            if channel_id>0
                StimTypeFR_AllChannels(i).FEF.Sample_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,2),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
                StimTypeFR_AllChannels(i).FEF.Delay_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,3),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
                StimTypeFR_AllChannels(i).FEF.Target_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,4),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
                StimTypeFR_AllChannels(i).FEF.Response_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,5),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
            else
                ssprintf('missing %d in %s for %s ',chi, arrayID_list(i), SD)
                StimTypeFR_AllChannels(i).FEF.Sample_on(chi)=NaN;
                StimTypeFR_AllChannels(i).FEF.Delay_on(chi)=NaN;
                StimTypeFR_AllChannels(i).FEF.Target_on(chi)=NaN;
                StimTypeFR_AllChannels(i).FEF.Response_on(chi)=NaN;
            end
            
            % channels sig
            StimTypeFR_AllChannels(i).FEF.Sample_on_005=StimTypeFR_AllChannels(i).FEF.ActiveChannels(find(StimTypeFR_AllChannels(i).FEF.Sample_on<anova_sig));
            StimTypeFR_AllChannels(i).FEF.Delay_on_005=StimTypeFR_AllChannels(i).FEF.ActiveChannels(find(StimTypeFR_AllChannels(i).FEF.Delay_on<anova_sig));
            StimTypeFR_AllChannels(i).FEF.Target_on_005=StimTypeFR_AllChannels(i).FEF.ActiveChannels(find(StimTypeFR_AllChannels(i).FEF.Target_on<anova_sig));
            StimTypeFR_AllChannels(i).FEF.Response_on_005=StimTypeFR_AllChannels(i).FEF.ActiveChannels(find(StimTypeFR_AllChannels(i).FEF.Response_on<anova_sig));
            % channels sig + coupled channels
            StimTypeFR_AllChannels(i).FEF.Sample_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).FEF.ActiveChannels, StimTypeFR_AllChannels(i).FEF.Sample_on_005);
            StimTypeFR_AllChannels(i).FEF.Delay_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).FEF.ActiveChannels, StimTypeFR_AllChannels(i).FEF.Delay_on_005);
            StimTypeFR_AllChannels(i).FEF.Target_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).FEF.ActiveChannels, StimTypeFR_AllChannels(i).FEF.Target_on_005);
            StimTypeFR_AllChannels(i).FEF.Response_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).FEF.ActiveChannels, StimTypeFR_AllChannels(i).FEF.Response_on_005);
        end
        
        for chi=1:length(StimTypeFR_AllChannels(i).PFC.ActiveChannels)
            channel_id=0;
            for j=1:length(StimTypeFR_AllChannels(i).StimTypeFR)
                if StimTypeFR_AllChannels(i).StimTypeFR(j).ID==StimTypeFR_AllChannels(i).PFC.ActiveChannels(chi)
                    channel_id=j;
                end
            end
            if channel_id>0
                StimTypeFR_AllChannels(i).PFC.Sample_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,2),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
                StimTypeFR_AllChannels(i).PFC.Delay_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,3),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
                StimTypeFR_AllChannels(i).PFC.Target_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,4),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
                StimTypeFR_AllChannels(i).PFC.Response_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,5),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
            else
                sprintf('missing %d in %s for %s ',chi, arrayID_list(i), SD)
                StimTypeFR_AllChannels(i).PFC.Sample_on(chi)=NaN;
                StimTypeFR_AllChannels(i).PFC.Delay_on(chi)=NaN;
                StimTypeFR_AllChannels(i).PFC.Target_on(chi)=NaN;
                StimTypeFR_AllChannels(i).PFC.Response_on(chi)=NaN;
            end
            % channels sig
            StimTypeFR_AllChannels(i).PFC.Sample_on_005=StimTypeFR_AllChannels(i).PFC.ActiveChannels(find(StimTypeFR_AllChannels(i).PFC.Sample_on<anova_sig));
            StimTypeFR_AllChannels(i).PFC.Delay_on_005=StimTypeFR_AllChannels(i).PFC.ActiveChannels(find(StimTypeFR_AllChannels(i).PFC.Delay_on<anova_sig));
            StimTypeFR_AllChannels(i).PFC.Target_on_005=StimTypeFR_AllChannels(i).PFC.ActiveChannels(find(StimTypeFR_AllChannels(i).PFC.Target_on<anova_sig));
            StimTypeFR_AllChannels(i).PFC.Response_on_005=StimTypeFR_AllChannels(i).PFC.ActiveChannels(find(StimTypeFR_AllChannels(i).PFC.Response_on<anova_sig));
            % channels sig + coupled channels
            StimTypeFR_AllChannels(i).PFC.Sample_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).PFC.ActiveChannels, StimTypeFR_AllChannels(i).PFC.Sample_on_005);
            StimTypeFR_AllChannels(i).PFC.Delay_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).PFC.ActiveChannels, StimTypeFR_AllChannels(i).PFC.Delay_on_005);
            StimTypeFR_AllChannels(i).PFC.Target_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).PFC.ActiveChannels, StimTypeFR_AllChannels(i).PFC.Target_on_005);
            StimTypeFR_AllChannels(i).PFC.Response_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).PFC.ActiveChannels, StimTypeFR_AllChannels(i).PFC.Response_on_005);
            
            
        end
        
        for chi=1:length(StimTypeFR_AllChannels(i).LIP.ActiveChannels)
            channel_id=0;
            for j=1:length(StimTypeFR_AllChannels(i).StimTypeFR)
                if StimTypeFR_AllChannels(i).StimTypeFR(j).ID==StimTypeFR_AllChannels(i).LIP.ActiveChannels(chi)
                    channel_id=j;
                end
            end
            if channel_id>0
            StimTypeFR_AllChannels(i).LIP.Sample_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,2),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
            StimTypeFR_AllChannels(i).LIP.Delay_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,3),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
            StimTypeFR_AllChannels(i).LIP.Target_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,4),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
            StimTypeFR_AllChannels(i).LIP.Response_on(chi)=anovan(StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,5),StimTypeFR_AllChannels(i).StimTypeFR(channel_id).Spike_count(:,1),'display','off');
            else
                sprintf('missing %d in %s for %s ',chi, arrayID_list(i), SD)
                StimTypeFR_AllChannels(i).LIP.Sample_on(chi)=NaN;
                StimTypeFR_AllChannels(i).LIP.Delay_on(chi)=NaN;
                StimTypeFR_AllChannels(i).LIP.Target_on(chi)=NaN;
                StimTypeFR_AllChannels(i).LIP.Response_on(chi)=NaN;
            end
            % channels sig
            StimTypeFR_AllChannels(i).LIP.Sample_on_005=StimTypeFR_AllChannels(i).LIP.ActiveChannels(find(StimTypeFR_AllChannels(i).LIP.Sample_on<anova_sig));
            StimTypeFR_AllChannels(i).LIP.Delay_on_005=StimTypeFR_AllChannels(i).LIP.ActiveChannels(find(StimTypeFR_AllChannels(i).LIP.Delay_on<anova_sig));
            StimTypeFR_AllChannels(i).LIP.Target_on_005=StimTypeFR_AllChannels(i).LIP.ActiveChannels(find(StimTypeFR_AllChannels(i).LIP.Target_on<anova_sig));
            StimTypeFR_AllChannels(i).LIP.Response_on_005=StimTypeFR_AllChannels(i).LIP.ActiveChannels(find(StimTypeFR_AllChannels(i).LIP.Response_on<anova_sig));
            % channels sig + coupled channels
            StimTypeFR_AllChannels(i).LIP.Sample_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).LIP.ActiveChannels, StimTypeFR_AllChannels(i).LIP.Sample_on_005);
            StimTypeFR_AllChannels(i).LIP.Delay_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).LIP.ActiveChannels, StimTypeFR_AllChannels(i).LIP.Delay_on_005);
            StimTypeFR_AllChannels(i).LIP.Target_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).LIP.ActiveChannels, StimTypeFR_AllChannels(i).LIP.Target_on_005);
            StimTypeFR_AllChannels(i).LIP.Response_on_005_both=add_coupled_channels(StimTypeFR_AllChannels(i).LIP.ActiveChannels, StimTypeFR_AllChannels(i).LIP.Response_on_005);
            
            
        end
        
    end
    
    save(fullfile(destpath,'StimTypeFR_AllChannels_restricted_FEF'),'StimTypeFR_AllChannels');
    
    
    %%Compute for all sessions
    
    
    Stats.LIP.total=0;
    Stats.LIP.delay_on=0;
    Stats.LIP.response_on=0;
    Stats.LIP.delay_on_both=0;
    Stats.LIP.response_on_both=0;
    
    for i=1:length(arrayID_list)
        Stats.LIP.total=Stats.LIP.total+length(StimTypeFR_AllChannels(i).LIP.ActiveChannels);
        Stats.LIP.delay_on=Stats.LIP.delay_on+length(StimTypeFR_AllChannels(i).LIP.Delay_on_005);
        Stats.LIP.response_on=Stats.LIP.response_on+length(StimTypeFR_AllChannels(i).LIP.Response_on_005);
        Stats.LIP.delay_on_both=Stats.LIP.delay_on_both+length(~isnan(StimTypeFR_AllChannels(i).LIP.Delay_on_005_both));
        Stats.LIP.response_on_both=Stats.LIP.response_on_both+length(~isnan(StimTypeFR_AllChannels(i).LIP.Response_on_005_both));
        
    end
    
    Stats.FEF.total=0;
    Stats.FEF.delay_on=0;
    Stats.FEF.response_on=0;
    Stats.FEF.delay_on_both=0;
    Stats.FEF.response_on_both=0;
    
    
    for i=1:length(arrayID_list)
        if ~isempty(StimTypeFR_AllChannels(i).FEF.ActiveChannels)
            Stats.FEF.total=Stats.FEF.total+length(StimTypeFR_AllChannels(i).FEF.ActiveChannels);
            Stats.FEF.delay_on=Stats.FEF.delay_on+length(StimTypeFR_AllChannels(i).FEF.Delay_on_005);
            Stats.FEF.response_on=Stats.FEF.response_on+length(StimTypeFR_AllChannels(i).FEF.Response_on_005);
            Stats.FEF.delay_on_both=Stats.FEF.delay_on_both+length(~isnan(StimTypeFR_AllChannels(i).FEF.Delay_on_005_both));
            Stats.FEF.response_on_both=Stats.FEF.response_on_both+length(~isnan(StimTypeFR_AllChannels(i).FEF.Response_on_005_both));
        end
    end
    
    Stats.PFC.total=0;
    Stats.PFC.delay_on=0;
    Stats.PFC.response_on=0;
    Stats.PFC.delay_on_both=0;
    Stats.PFC.response_on_both=0;
    
    
    for i=1:length(arrayID_list)
        Stats.PFC.total=Stats.PFC.total+length(StimTypeFR_AllChannels(i).PFC.ActiveChannels);
        Stats.PFC.delay_on=Stats.PFC.delay_on+length(StimTypeFR_AllChannels(i).PFC.Delay_on_005);
        Stats.PFC.response_on=Stats.PFC.response_on+length(StimTypeFR_AllChannels(i).PFC.Response_on_005);
        Stats.PFC.delay_on_both=Stats.PFC.delay_on_both+length(~isnan(StimTypeFR_AllChannels(i).PFC.Delay_on_005_both));
        Stats.PFC.response_both=Stats.PFC.response_on_both+length(~isnan(StimTypeFR_AllChannels(i).PFC.Response_on_005_both));
        
    end
    
    
    save(fullfile(destpath,'Stats_restricted_FEF'),'Stats');
    
end












