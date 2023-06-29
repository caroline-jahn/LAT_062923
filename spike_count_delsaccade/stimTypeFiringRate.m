function stimTypeFiringRate(fsroot, task, arrayID, SD)

%% Create a strcuture that contains the FR for each channels and each stage

%StimTypeFR.Spike_count for correct trials only :   col 1 : target location (1 to 8)
%                                                   col 2 : sample on (stage 2)
%                                                   col 3 : delay on (from stage 3 to stage 4)
%                                                   col 4 : target on (stage 4)
%                                                   col 5 : response (stage 4 in NEV + RT)
% ALL FR ARE NORMALISED BY WINDOW DURATION !!!!!

%StimTypeFR.Window_duration for correct trials only : col 1 : sample on
%                                                     col 2 : delay on (from stage 3 to stage 4)
%                                                     col 3 : target on 
%                                                     col 4 : response 

%% Param and path

window_sample_on=[0 100]; %in ms
window_target_on=[0 200];
window_reponse_on=[-50 150];

window_sample_on=window_sample_on/1000; %in sec
window_target_on=window_target_on/1000;
window_reponse_on=window_reponse_on/1000;


addpath(genpath(fullfile(filesep,fsroot,'Users','Caroline',task,'Electrophy_analysis','spike_count_delsaccade','NPMK')));

dirstem = fullfile('Projects',task,'Data','General');
dirpath = fullfile(fsroot,dirstem);
load(fullfile(dirpath,'NS6Directory.mat'));

entry = ns6directory(arrayID);

subtask='delsacc'; %% Make sure that there is only ONE relevant behavior file in bhvg folder

    nSD   = SD; %how many sd to use as cutoff?
	sdstr = sprintf('%.3dsd', round(10*nSD));



%% data path
datpath=fullfile(fullfile(fsroot,entry.FolderStem));

destpath = fullfile(datpath,'MUA','stimTypeFR');
mkdir(destpath);
matds = dir(fullfile(datpath, 'bhv' ,sprintf('%s_%s_*bhv.mat', entry.Subject, subtask)));  %%% IS SUBTASK DEPENDENT here Delayed saccade task !!!!!
cutds = dir(fullfile(datpath,'CutTrials',sprintf('CutTrials_%s_%s*.mat',  entry.Subject, subtask)));
muads = dir(fullfile(datpath, 'MUA', sdstr, sprintf('*%s*chan*.mat', subtask)));
rawds = dir(fullfile(datpath, entry.FileName));
if length(matds) ~= 1 || length(cutds) ~= 1 || isempty(muads) || length(rawds) ~= 1
    fprintf('Not appropriate number of bhv, cut, NS6, or MUA files or directories\n');
    return;
end


%% Read files [trials, cutTrials, Ts, header]
load(fullfile(datpath,'bhv',matds.name),'trials','cache','opts');
load(fullfile(datpath,'CutTrials',cutds.name));
Ts = cell(1,length(muads));
for chi=1:length(muads)
    tempTs = load(fullfile(datpath,'MUA',sdstr,muads(chi).name),'Ts');
    Ts{chi} = tempTs.Ts{1};
end
head = openNSx(fullfile(datpath,rawds.name),'noread'); % correct Ts in sec usinf sampling freq
for chi=1:length(muads)
    Ts_sec{chi}=Ts{chi}/head.MetaTags.SamplingFreq; % Because of sampling frequencey, could use Pavlos' method and find spike(ind) is the timestamp of the segment in which that spike has occured + the offset of that spike factoring out previous segments
end
%struct for storing times and FRs for each stim type for each stim for
%each channel
stimTypeFR = repmat(struct('ID', NaN, 'Spike_count', NaN(length(trials),5),'Window_duration',ones(length(trials),4)),1,length(Ts)); %check dimension
for chi=1:length(Ts)
    stimTypeFR(chi).ID = str2num(muads(chi).name(end-12:end-10));  %ID =channel number
end

%% For each channel
for chi=1:length(Ts)
    %for each trial
    for ti=1:length(trials)
        %if correct trial
        if ismember(cache.lists.StopCondition(ti),opts.Behavior.StopConditionCorrect) && ~isnan(cutTrials(ti).Timing.NS4.TimeStampSec(1)) %only correct trials => 5 stages and response is in direction
            %assign stim direction
            stimTypeFR(chi).Spike_count(ti,1)=trials(ti).Conditions.TargetLocation;
            
            %computes spike count (ignoring the sampling problem) and
            %windows
            %Sample_on = stage 2
            stimTypeFR(chi).Spike_count(ti,2)=size(find(Ts_sec{chi}>cutTrials(ti).Timing.NS4.TimeStampSec(2)+window_sample_on(1) & Ts_sec{chi}<cutTrials(ti).Timing.NS4.TimeStampSec(2)+window_sample_on(2)),1);
            stimTypeFR(chi).Window_duration(ti,1)=window_sample_on(2)-window_sample_on(1);
            stimTypeFR(chi).Spike_count(ti,2)=stimTypeFR(chi).Spike_count(ti,2)/stimTypeFR(chi).Window_duration(ti,1);
            %Delay_on = stage 3 until stage 4
            stimTypeFR(chi).Spike_count(ti,3)=size(find(Ts_sec{chi}>cutTrials(ti).Timing.NS4.TimeStampSec(3) & Ts_sec{chi}<cutTrials(ti).Timing.NS4.TimeStampSec(4)),1);
            stimTypeFR(chi).Window_duration(ti,2)=cutTrials(ti).Timing.NS4.TimeStampSec(4)-cutTrials(ti).Timing.NS4.TimeStampSec(3);
            stimTypeFR(chi).Spike_count(ti,3)=stimTypeFR(chi).Spike_count(ti,3)/stimTypeFR(chi).Window_duration(ti,2);
            %Target on = stage 4
            stimTypeFR(chi).Spike_count(ti,4)=size(find(Ts_sec{chi}>cutTrials(ti).Timing.NS4.TimeStampSec(4)+window_target_on(1) & Ts_sec{chi}<cutTrials(ti).Timing.NS4.TimeStampSec(4)+window_target_on(2)),1);
            stimTypeFR(chi).Window_duration(ti,3)=window_target_on(2)-window_target_on(1);
            stimTypeFR(chi).Spike_count(ti,4)=stimTypeFR(chi).Spike_count(ti,4)/stimTypeFR(chi).Window_duration(ti,3);
            %Response => time recorded in NEV
            stimTypeFR(chi).Spike_count(ti,5)=size(find(Ts_sec{chi}>cutTrials(ti).Timing.NEV.RESPONSE_ON_Sec+window_reponse_on(1) & Ts_sec{chi}<cutTrials(ti).Timing.NEV.RESPONSE_ON_Sec+window_reponse_on(2)),1);
            stimTypeFR(chi).Window_duration(ti,4)=window_reponse_on(2)-window_reponse_on(1);
            stimTypeFR(chi).Spike_count(ti,5)=stimTypeFR(chi).Spike_count(ti,5)/stimTypeFR(chi).Window_duration(ti,4);
        end
    end
end



    save(fullfile(destpath,sprintf('%s_%s_segmented.mat', rawds.name(1:end-4),sdstr)),'stimTypeFR');

end