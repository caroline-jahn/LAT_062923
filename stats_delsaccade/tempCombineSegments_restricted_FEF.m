function StimTypeFR_AllChannels=tempCombineSegments_restricted_FEF(fsroot,task,subtask,arrayID_list, SD)

deststem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask);
destpath = fullfile(fsroot,deststem);
mkdir(destpath);

    nSD   = SD; %how many sd to use as cutoff?
	sdstr = sprintf('%.3dsd', round(10*nSD));


load(fullfile(fsroot,'Projects',task,'Data','General','NS6Directory_AC_restricted_FEF.mat')); 

for index_arrayID=1:length(arrayID_list)
    
    entry=ns6directory(arrayID_list(index_arrayID));
    
    load(fullfile(fsroot,entry.FolderStem,'MUA','stimTypeFR', sprintf('%s_%s_segmented.mat',entry.FileName(1:end-4),sdstr)),'stimTypeFR');

     StimTypeFR_AllChannels(index_arrayID).ArrayID=arrayID_list(index_arrayID);
     StimTypeFR_AllChannels(index_arrayID).FileName=entry.FileName(1:end-4);
     StimTypeFR_AllChannels(index_arrayID).StimTypeFR=stimTypeFR;
     
     clear stimTypeFR

end

end