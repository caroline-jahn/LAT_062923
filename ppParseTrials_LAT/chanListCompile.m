function chanListCompile(fsroot, task)

subjects = {'Beaker' 'Scooter'};
% task = 'Learning_Attentional_Templates';
% fsroot = fullfile(filesep,'Volumes','buschman'); %fullfile(filesep,'Volumes','buschman'), 'Z:' fullfile(filesep,'jukebox','buschman')

deststem = fullfile('Projects',task,'Data','General');
destpath = fullfile(fsroot,deststem);
mkdir(destpath);

entry = struct('FileName', NaN, ...
               'FolderStem', NaN, ...
               'Subject', NaN, ...
               'ChannelCount', NaN, ...
               'ChannelID', NaN, ...
               'Semaphore', NaN);
ns6directory = repmat(entry,1,0);

for si=1:length(subjects),
    
    subject = subjects{si};
    srcstem = fullfile('Projects',task,'Data',subject);
    srcpath = fullfile(fsroot,srcstem);
    ds = dir(srcpath);
    
    for di=1:length(ds),
        
        if ~isempty(str2num(ds(di).name)) && length(ds(di).name) == 6,
            srcdatpath = fullfile(srcpath,ds(di).name);
            srcsortpath = fullfile(srcdatpath,'CellSorting');
            mkdir(srcsortpath);
            
            
            ns6ds = dir(fullfile(srcdatpath,'*.ns6'));
            for nsi=1:length(ns6ds),
                dat = openNSx(fullfile(srcdatpath,ns6ds(nsi).name),'noread');
                tmp = entry;
                tmp.FileName = ns6ds(nsi).name;
                tmp.FolderStem = fullfile(srcstem,ds(di).name);
                tmp.Subject = subject;
                tmp.ChannelCount = dat.MetaTags.ChannelCount;
                tmp.ChannelID = dat.MetaTags.ChannelID;
                ns6directory(end+1) = tmp;
            end
        end
    end
end

    ns6directory = makesemaphore(ns6directory); 


save(fullfile(destpath, 'NS6Directory_sem'), 'ns6directory');

end