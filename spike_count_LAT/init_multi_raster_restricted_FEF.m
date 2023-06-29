function [trials, cutTrials, Spikes_time_serie, opts, specs]=init_multi_raster_restricted_FEF(arrayID)



task = 'Learning_Attentional_Templates';
fsroot = fullfile(filesep,'Volumes','buschman'); %fullfile(filesep,'Volumes','buschman'), 'Z:' fullfile(filesep,'jukebox','buschman')
subtask='exploreexploit'; %% Make sure that there is only ONE relevant behavior file in bhv folder


%% Path

dirstem = fullfile('Projects',task,'Data','General');
dirpath = fullfile(fsroot,dirstem);
load(fullfile(dirpath,'NS6Directory_SC_restricted_FEF.mat')); %carries info about LIP sensitive channels

entry = ns6directory(arrayID);



%% data path
datpath=fullfile(fullfile(fsroot,entry.FolderStem));

matds = dir(fullfile(datpath, 'bhv' ,sprintf('%s_18*bhv.mat', entry.Subject)));  %%% IS SUBTASK DEPENDENT here NOT Delayed saccade task !!!!!
cutds = dir(fullfile(datpath,'CutTrials',sprintf('CutTrials_%s_18*.mat',  entry.Subject)));
spikeds = dir(fullfile(datpath, 'Spikes','*chan*-srt.mat'));
if length(matds) ~= 1 || length(cutds) ~= 1 || isempty(spikeds)
    fprintf('Not appropriate number of bhv, cut, spikes files or directories\n');
    return;
end


%% Read files [trials, cutTrials, Ts] and organise per unit
load(fullfile(datpath,'bhv',matds.name),'trials', 'opts', 'specs');
load(fullfile(datpath,'CutTrials',cutds.name));
for chi=1:length(spikeds)
    Spikes_time_serie(chi).Ts = cell(9,1);
    Spikes_time_serie(chi).Wv = cell(9,1);
end %max 9 units per channel, completely arbitrary

Nb_units=zeros(1,length(spikeds));
for chi=1:length(spikeds)
    TempTs = load(fullfile(datpath,'Spikes',spikeds(chi).name),'id','ts','wv');
    Nb_units(chi)=max(TempTs.id);
    for n_unit=1:max(TempTs.id)
        Spikes_time_serie(chi).Ts{n_unit} = TempTs.ts(TempTs.id==n_unit);
        Spikes_time_serie(chi).Wv{n_unit} = TempTs.wv(TempTs.id==n_unit,:);
    end
    Spikes_time_serie(chi).ID = str2num(spikeds(chi).name(end-10:end-8));
end

