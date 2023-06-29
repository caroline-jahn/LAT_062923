function chanSplit(fsroot, task, arrayID)

%     fsroot = 'Z:'; %'Z:';
%                    %fullfile(filesep,'Volumes','buschman');
%                    %fullfile(filesep,'jukebox','buschman');
%     task = 'GatingInWorkingMemory';

    dirstem = fullfile('Projects',task,'Data','General');
    dirpath = fullfile(fsroot,dirstem);
    load(fullfile(dirpath,'NS6Directory.mat'));

    entry = ns6directory(arrayID);

    getme = fullfile(fullfile(fsroot,entry.FolderStem,entry.FileName));
    meta = openNSx('noread',getme);

    for mx=1:meta.MetaTags.ChannelCount,
%         try,
            tmp = openNSx ('read',getme,sprintf('%s%i','c:',mx));
            fprintf('Attempting Channel %d/%d (%d)\n', mx, meta.MetaTags.ChannelCount, tmp.MetaTags.ChannelID);
            chanNSxSave(tmp, fullfile(fsroot,entry.FolderStem,'CellSorting', ...
                        sprintf('%s_chan%.3d.ns6',tmp.MetaTags.Filename,tmp.MetaTags.ChannelID)));
            clear tmp;
%          catch,
%              fprintf('Failed\n');
%          end
    end
    clear meta;

end