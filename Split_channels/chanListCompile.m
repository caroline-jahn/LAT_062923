subjects = {'Beaker' 'Scooter'};
task = 'Learning_Attentional_Templates';
fsroot = fullfile(filesep,'Volumes','buschman'); %fullfile(filesep,'Volumes','buschman'), 'Z:' fullfile(filesep,'jukebox','buschman')

deststem = fullfile('Projects',task,'Data','General');
destpath = fullfile(fsroot,deststem);
mkdir(destpath);

entry = struct('FileName', NaN, ...
               'FolderStem', NaN, ...
               'Subject', NaN, ...
               'ChannelCount', NaN, ...
               'ChannelID', NaN);
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
            
            %fix all channel names
%             subd = [dir(fullfile(srcsortpath,'*ns6')); ...
%                     dir(fullfile(srcsortpath,'*plx')); ...
%                     dir(fullfile(srcsortpath,'*mat')); ...
%                     dir(fullfile(srcsortpath,'*png'))];
%             if ~isempty(subd),
%                 for sdi=1:length(subd),
%                     %fix all ns6, plx, png, and mat
%                     chanstr = subd(sdi).name(1:(strfind(subd(sdi).name,'chan')+3));
%                     idstr = sprintf('%.3d',str2num(subd(sdi).name((strfind(subd(sdi).name,'chan')+4):(strfind(subd(sdi).name,'.')-1))));
%                     extstr = subd(sdi).name(strfind(subd(sdi).name,'.'):end);
%                     newname = strcat(chanstr,idstr,extstr);
%                     oldfile = fullfile(srcsortpath,subd(sdi).name);
%                     newfile = fullfile(srcsortpath,newname);
%                     if ~strcmp(oldfile,newfile),
%                         movefile(fullfile(srcsortpath,subd(sdi).name),fullfile(srcsortpath,newname));
%                     end
%                 end
%             end
            
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

save(fullfile(destpath, 'NS6Directory'), 'ns6directory');

%13:38 saveChNSx -> :27
%13:50 openNSx -> :27