function extract_data_for_belief_analysis_Reset_RW_restricted_FEF(fsroot,event,window_start,window_duration)

arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44];

load(fullfile(fsroot,'/Projects/Learning_Attentional_Templates/Data/General/NS6Directory_SC_restricted_FEF.mat'))

switch event
    case 'target' 
    incode=2;
    case 'response'    
        incode=3;
    case 'fixation'
        incode=1;
    case 'reward'
        incode=4;
end

    window_end=window_start+window_duration;

window_step=1;
nb_windows=1;

switch_type='Reset';
N_channels=6;
model_type='RW';

save_path=fullfile('/Volumes/buschman/Projects/Learning_Attentional_Templates/Analysis/Electrophy_analysis/exploreexploit/Restricted_FEF/Reset_RW_model',sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_duration,window_start,window_end));
mkdir(save_path);

%% raster script:

for arrayID_ind=1:length(arrayID_list)
    
    arrayID=arrayID_list(arrayID_ind)
    monkey=ns6directory(arrayID).Subject;
    date=ns6directory(arrayID).FolderStem(end-5:end);
        
    model_path=fullfile(fsroot,'/Projects/Learning_Attentional_Templates/Analysed_data');
    model_name=sprintf('%s_%s_%d_channels_VBMC',switch_type,model_type,N_channels);
    
    this_session_path=fullfile(model_path,monkey,date,model_name);
    load(this_session_path);
    
    if ~isnan(ns6directory(arrayID).LIP_sensitive)
        channel_id_list=ns6directory(arrayID).LIP_sensitive;
    else
        channel_id_list=ns6directory(arrayID).PFC;
    end
    
    channel_id=channel_id_list(1);
    
    [trials, cutTrials, Spikes_time_serie, ~, ~]=init_multi_raster_restricted_FEF(arrayID);
    list_conditions=[1 -1];
    
    for index_in_sts=1:length(Spikes_time_serie)
        if channel_id==Spikes_time_serie(index_in_sts).ID
            channel_nb=index_in_sts;
            continue
        end
    end
    
    %% Add behavioral info
    Behavior=NaN(length(trials),5);
    for ti=1:length(trials)
        %if correct trial
        if  ismember(trials(ti).Behavior.StopCondition,list_conditions) && ~isnan(cutTrials(ti).Timing.NS4.TimeStampSec(1)) %only correct trials => 5 stages and response is in direction .  %ismember(cache.lists.StopCondition(ti),opts.Behavior.StopConditionCorrect)
            %assign stim direction
            Behavior(ti,1)=trials(ti).Behavior.ChosenColor; %response = location chosen (1 to 4)
            Behavior(ti,2)=trials(ti).Behavior.ResponseError;
            Behavior(ti,3)=trials(ti).Conditions.BlockNo; %first put block number and find switch point after
            Behavior(ti,4)=trials(ti).Conditions.BlockNo; %first put block number and find switch point after
            Behavior(ti,5)=trials(ti).Conditions.BestColor; %best Nb_block (in rad)
            Behavior(ti,6)=trials(ti).Conditions.TargetColorToLocIndex(trials(ti).Behavior.Response); %response = location chosen (1 to 4)
            for j_s=1:4
                if sum(trials(ti).Conditions.TargetColorToLocIndex(:))==5+j_s
                    Behavior(ti,7)=j_s;
                end
            end
            Behavior(ti,8:11)=NaN;
            for j_s=1:3
                Behavior(ti,7+trials(ti).Conditions.TargetColorToLocIndex(j_s))=trials(ti).Conditions.TargetColorAngles(j_s);
            end
            Behavior(ti,12)=trials(ti).Reward.NumRewards;
        end
    end
    
    %find block switches
    correct_trials=NaN(length(trials),1);
    count_trial=1;
    for ti=1:length(trials)
        if ~isnan(Behavior(ti,3)) && Behavior(ti,4)<max(Behavior(:,4)) 
            correct_trials(ti)= count_trial;
            count_trial=count_trial+1;
        end
    end
    
    block_sw=zeros(max(Behavior(:,3)),1);
    block_sw(1)=1;
    for i=2:length(block_sw)
        block_sw(i)=find(Behavior(:,3)==i,1,'first');
        block_sw(i)=correct_trials(block_sw(i));
    end
    
    for ti=1:length(trials)
        %if correct trial
        if ~isnan(Behavior(ti,3)) 
            Behavior(ti,3)= correct_trials(ti) - block_sw(Behavior(ti,3)) + 1;
        end
    end
    for ti=1:length(trials)    
            if Behavior(ti,4)==1 || isnan(Behavior(ti,4)) %first block
                Behavior(ti,13)=NaN;
            else
                Behavior(ti,13)=Behavior(find(Behavior(:,4)==Behavior(ti,4)-1,1,'first'),5);
            end
    end

    
        %% Create the conditions here

        angle_binned=0:2*pi/Model_predictions.model_specificities.N_bins:2*pi-2*pi/Model_predictions.model_specificities.N_bins;
        
        for ti=1:length(trials)
            if ~isnan(Behavior(ti,3)) && Behavior(ti,4)<max(Behavior(:,4))
                %use weight of the previous trial for hcoice in current
                %trial
                Conditions.weight(:,correct_trials(ti))=Model_predictions.model_outputs.Weight(:,correct_trials(ti));
                Conditions.belief(:,correct_trials(ti))=Model_predictions.model_outputs.Value_for_choice(:,correct_trials(ti));
                Conditions.down_sampled_belief(:,correct_trials(ti))=Model_predictions.model_outputs.down_sampled_Value_for_choice(:,correct_trials(ti));
                Conditions.switch(correct_trials(ti))=Model_predictions.model_outputs.Switch(:,correct_trials(ti));
                 Conditions.accuracy(correct_trials(ti))=Model_predictions.behavior.Best_chosen(correct_trials(ti));
                if Model_predictions.model_outputs.Precision(:,correct_trials(ti))>0
                    Conditions.mean_belief(correct_trials(ti))=angle(sum(Model_predictions.model_outputs.Value_for_choice(:,correct_trials(ti)).*exp(1i.*angle_binned')));
                else
                    Conditions.mean_belief(correct_trials(ti))=NaN;
                end
                Conditions.belief_precision(correct_trials(ti))=Model_predictions.model_outputs.Precision(:,correct_trials(ti));
                Conditions.trial_in_block(correct_trials(ti))=Behavior(ti,3);
                Conditions.block_nb(correct_trials(ti))=Behavior(ti,4);
                Conditions.rule(correct_trials(ti))=Behavior(ti,5);
                Conditions.old_rule(correct_trials(ti))=Behavior(ti,13);
                Conditions.chosen_color(correct_trials(ti))=Behavior(ti,1);
                Conditions.chosen_color_error(correct_trials(ti))=Behavior(ti,2);
                Conditions.choice(correct_trials(ti))=Behavior(ti,6);
                Conditions.choice_context(correct_trials(ti))=Behavior(ti,7);
                Conditions.stim(:,correct_trials(ti))=Behavior(ti,8:11);
                Conditions.RPE(correct_trials(ti))=Model_predictions.model_outputs.RPE(correct_trials(ti));
                Conditions.Reward(correct_trials(ti))=Behavior(ti,12)/max(Behavior(:,12));
                Conditions.RT(correct_trials(ti))=Model_predictions.model_input.in(correct_trials(ti),18);
                % reatribute value and P_choice to location
                for k=1:4 %each location
                    if Model_predictions.model_input.in(correct_trials(ti),13)==k %best at this location
                        Conditions.P_choice(k,correct_trials(ti))=Model_predictions.model_outputs.Choice_prediction(1,correct_trials(ti));                    
                    end
                    if Model_predictions.model_input.in(correct_trials(ti),14)==k %best at this location
                        Conditions.P_choice(k,correct_trials(ti))=Model_predictions.model_outputs.Choice_prediction(2,correct_trials(ti));        
                    end
                    if Model_predictions.model_input.in(correct_trials(ti),15)==k %best at this location
                        Conditions.P_choice(k,correct_trials(ti))=Model_predictions.model_outputs.Choice_prediction(3,correct_trials(ti));        
                    end                    
                end
            end
        end
        
    %LIP
    
    channel_id_list=ns6directory(arrayID).LIP_sensitive;
    
    if ~isnan(ns6directory(arrayID).LIP_sensitive)
          count=1;          
        for i=1:length(channel_id_list)
            
            channel_id=channel_id_list(i);
            
            for cellnb=1:9
                   
                flag=1;            
                if ~isempty(find([Spikes_time_serie.ID]==channel_id,1))
                    channel_nb=find([Spikes_time_serie.ID]==channel_id);
                else
                   flag=0; 
                end
                
                if ~isempty(Spikes_time_serie(channel_nb).Ts{cellnb,1}) && flag==1
                                        
                    for n_w=1:nb_windows
                        window=[window_start+(n_w-1)*window_step, window_start+window_duration+(n_w-1)*window_step];
                        [spkcount_tmp] = get_spike_count_ind(Spikes_time_serie(channel_nb).Ts, cutTrials, cellnb, incode, correct_trials, window);
                        LIP(count,:,n_w)=spkcount_tmp(:);
                        LIP_wv(count,:)=mean(Spikes_time_serie(channel_nb).Wv{cellnb,1},1);                        
                    end
                        count=count+1;                        
                end
            end
        end
        
        
        clear spkcount_tmp channel_id_list
        
    else
        
        LIP=NaN;
    end
    
        %FEF
        
        channel_id_list=ns6directory(arrayID).FEF;
        
    if ~isnan(ns6directory(arrayID).FEF)
          count=1;          
        for i=1:length(channel_id_list)
            
            channel_id=channel_id_list(i);
            
            for cellnb=1:9
                
                 flag=1;            
                if ~isempty(find([Spikes_time_serie.ID]==channel_id,1))
                    channel_nb=find([Spikes_time_serie.ID]==channel_id);
                else
                   flag=0; 
                end
                
                if ~isempty(Spikes_time_serie(channel_nb).Ts{cellnb,1}) && flag==1
                                       
                   for n_w=1:nb_windows
                        window=[window_start+(n_w-1)*window_step, window_start+window_duration+(n_w-1)*window_step];
                        [spkcount_tmp] = get_spike_count_ind(Spikes_time_serie(channel_nb).Ts, cutTrials, cellnb, incode, correct_trials, window);
                        FEF(count,:,n_w)=spkcount_tmp(:);
                        FEF_wv(count,:)=mean(Spikes_time_serie(channel_nb).Wv{cellnb,1},1);
                    end
                       count=count+1;                        
                end
            end
        end
        
        
        clear spkcount_tmp channel_id_list
        
    else
        
        FEF=NaN;
        
    end

        %PFC
        
        channel_id_list=ns6directory(arrayID).PFC;
        
    if ~isnan(ns6directory(arrayID).PFC)
          count=1;          
        for i=1:length(channel_id_list)
            
            channel_id=channel_id_list(i);
            
            for cellnb=1:9
                
                 flag=1;            
                if ~isempty(find([Spikes_time_serie.ID]==channel_id,1))
                    channel_nb=find([Spikes_time_serie.ID]==channel_id);
                else
                   flag=0; 
                end
                
                if ~isempty(Spikes_time_serie(channel_nb).Ts{cellnb,1}) && flag==1
                                       
                    for n_w=1:nb_windows
                        window=[window_start+(n_w-1)*window_step, window_start+window_duration+(n_w-1)*window_step];
                        [spkcount_tmp] = get_spike_count_ind(Spikes_time_serie(channel_nb).Ts, cutTrials, cellnb, incode, correct_trials, window);
                        PFC(count,:,n_w)=spkcount_tmp(:);
                        PFC_wv(count,:)=mean(Spikes_time_serie(channel_nb).Wv{cellnb,1},1);    
                    end
                        count=count+1;                        
                end
            end
        end
        
        clear spkcount_tmp channel_id_list
        
    else
        
        PFC=NaN;
        
    end

    save_name=sprintf('ID_%d.mat',arrayID);
       
    save(fullfile(save_path,save_name),'Conditions','LIP','FEF','PFC','LIP_wv','FEF_wv','PFC_wv')
    
    clear Conditions LIP FEF PFC Behavior ind* correct_trials
    
end
