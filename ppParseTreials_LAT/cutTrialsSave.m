function cutTrialsSave(fsroot, task, arrayID)

%     fsroot = 'Z:'; %'Z:';
%                    %fullfile(filesep,'Volumes','buschman');
%                    %fullfile(filesep,'jukebox','buschman');
%     task = 'GatingInWorkingMemory';

    addpath(genpath(fullfile(filesep,fsroot,'Users','Caroline',task,'Electrophy_analysis','ppParseTrials_exploreeploit','NPMK')));

    chanListCompile(fsroot, task);

    dirstem = fullfile('Projects',task,'Data','General');
    dirpath = fullfile(fsroot,dirstem);
    load(fullfile(dirpath,'NS6Directory_sem.mat'));
    
    entry = ns6directory(arrayID);
    
    subtask='exploreexploit'; %% Make sure that there is only ONE relevant behavior file in bhvg folder

    
    %% data path
    datpath=fullfile(fullfile(fsroot,entry.FolderStem));
    
    %find paths and load NEV, NS4, and mat files
    destpath = fullfile(datpath,'CutTrials');
    mkdir(destpath);
    nevName=[entry.FileName(1:end-3) 'nev'];
    nevds = dir(fullfile(datpath,nevName));
    ns4Name=[entry.FileName(1:end-3) 'ns4'];
    ns4ds = dir(fullfile(datpath,ns4Name));
    matds = dir(fullfile(datpath, 'bhv' ,sprintf('%s_18*bhv.mat', entry.Subject)));  %%% IS SUBTASK DEPENDENT here Explore Exploit task !!!!! 
    %exception catch
    if length(nevds) ~= 1 || length(matds) ~= 1 || length(ns4ds) ~= 1
        fprintf('%s Session %s has exceptional number of NEV, NS4, or mat files\n',entry.FileName);
        return;
    end
    NEV = openNEV(fullfile(datpath,nevds.name));
    NS4 = openNSx(fullfile(datpath,ns4ds.name));
    load(fullfile(datpath, 'bhv', matds.name),'trials','specs');

    %cut trials
    [cutTrials, cutSpecs] = cutNEVTrial(NEV,trials,specs,entry.Semaphore); %NEV = Blackrocktime
    if arrayID>1
    [cutTrials, cutSpecs] = cutNS4Trial(NS4,cutTrials,cutSpecs,arrayID); %NS4 = Actual time of dislay (based on photosensor) Need arrayID to debug 46
    end

%     %add other trial info
%     [cutTrials] = add_info(cutTrials,trials);

    %save cut
    save(fullfile(destpath,sprintf('CutTrials_%s.mat',matds.name(1:end-4))),'cutTrials','cutSpecs'); 
    
end