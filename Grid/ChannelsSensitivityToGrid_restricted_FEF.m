%% ChannelsSensitivityToGrid
%From selective channels on each days give you the correponding grid location
%Warning: grids' orientations have to be the same over the days of
%recordings OR need to had a rotation (not coded here)

addpath('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Electrophy_Analysis/Grid');
addpath('/Volumes/buschman/Projects/Learning_Attentional_Templates/Data/');

load(fullfile(destpath,'NS6Directory_AC_restricted_FEF.mat'));
%choose the nSD you want to analyse
load('/Volumes/buschman/Projects/Learning_Attentional_Templates/Analysis/Electrophy_analysis/delsacc/022sd/StimTypeFR_AllChannels.mat')

%% OPTIONS
PLOT=1;

%% Grids
%Check with Matt's scipt to localise

%Parital grid
LIP_b=[3 6; 4 5; 4 6; 4 7; 5 4; 5 5; 5 6; 3 5; 4 4];
LIP_s=[1 6; 2 4; 2 5; 2 6; 3 2; 3 3; 3 4; 3 5; 4 1; 4 2; 4 3];



%% Sessions for Beaker

arrayID_list_beaker=[4 9 13 15 19 21 23 27];
dates_list_beaker={'101918','102518','102918','103018','110618','110918','111018','111518'}; %'102618', %,'111418','111518' %,'111418' '101718',

for nday=1:length(arrayID_list_beaker)
    for i=1:length(StimTypeFR_AllChannels)
        if StimTypeFR_AllChannels(i).ArrayID==arrayID_list_beaker(nday)
            index_nday=i;
            continue
        end
    end
    active_channels_parietal_beaker{1,nday}=StimTypeFR_AllChannels(index_nday).LIP.ActiveChannels;
    sensitive_channels_parietal_beaker{1,nday}=StimTypeFR_AllChannels(index_nday).LIP.Delay_on_005_both;
end


%% Sessions for Scooter

arrayID_list_scooter=[29 31 33 35 37 39 41 43 45 47];
dates_list_scooter={'110218','110318','110718','110818','111218','111318','111618','111718','111918','112618'}; %%,'113018'

for nday=1:length(arrayID_list_scooter)
    for i=1:length(StimTypeFR_AllChannels)
        if StimTypeFR_AllChannels(i).ArrayID==arrayID_list_scooter(nday)
            index_nday=i;
            continue
        end
    end
    active_channels_parietal_scooter{1,nday}=StimTypeFR_AllChannels(index_nday).LIP.ActiveChannels;
    sensitive_channels_parietal_scooter{1,nday}=StimTypeFR_AllChannels(index_nday).LIP.Delay_on_005_both;
end



%% Inititalise variables

%Beaker
units_P1_b=NaN(6,8,size(dates_list_beaker,2)); %all recorded channels
units_P2_b=NaN(6,8,size(dates_list_beaker,2));

active_units_P1_b=NaN(6,8,size(dates_list_beaker,2)); % active channels
active_units_P2_b=NaN(6,8,size(dates_list_beaker,2));

sensitive_units_P1_b=NaN(6,8,size(dates_list_beaker,2)); % selective during delay on
sensitive_units_P2_b=NaN(6,8,size(dates_list_beaker,2));

depth_P_b=zeros(6,8,size(dates_list_beaker,2));

%Scooter
units_P1_s=NaN(6,8,size(dates_list_scooter,2)); %channel matrix (electrod number, 3rd dim = 2 electrodes/position)
units_P2_s=NaN(6,8,size(dates_list_scooter,2));

active_units_P1_s=NaN(6,8,size(dates_list_scooter,2));
active_units_P2_s=NaN(6,8,size(dates_list_scooter,2));

sensitive_units_P1_s=NaN(6,8,size(dates_list_scooter,2));
sensitive_units_P2_s=NaN(6,8,size(dates_list_scooter,2));

depth_P_s=zeros(6,8,size(dates_list_scooter,2));

%Color map of brain areas on grid
Parietal_area_b=NaN(6,8);
Parietal_area_s=NaN(6,8);

%Total units
total_active_units_P_b=zeros(6,8);
total_active_units_P_s=zeros(6,8);

total_sensitive_units_P_b=zeros(6,8);
total_sensitive_units_P_s=zeros(6,8);



%% Beaker

cd('/Volumes/buschman/Projects/Learning_Attentional_Templates/Data/ElectrodeTracking');

for nday=1:size(dates_list_beaker,2)
    
    load(['LearnAttTemplate_Beaker_' dates_list_beaker{1,nday} '_ElectrodeTracks_Lowered.mat'],'elec_track')
    
    [   units_P1_b, units_P2_b,...
        active_units_P1_b, active_units_P2_b,...
        sensitive_units_P1_b,sensitive_units_P2_b] = make_sensitive_parietal_grid( elec_track, nday, ...
        active_channels_parietal_beaker{1,nday}, sensitive_channels_parietal_beaker{1,nday}, ...
        units_P1_b, units_P2_b,...
        active_units_P1_b, active_units_P2_b,...
        sensitive_units_P1_b, sensitive_units_P2_b);
    
    [ depth_P_b] = compute_depth_parietal( elec_track, nday, depth_P_b);
    
    clear elec_track
    
end

%Compute total active channels for each position on grid
for i=1:6
    for j=1:8
        total_active_units_P_b(i,j)=nansum(active_units_P1_b(i,j,:),3)+nansum(active_units_P2_b(i,j,:),3);
        total_sensitive_units_P_b(i,j)=nansum(sensitive_units_P1_b(i,j,:),3)+nansum(sensitive_units_P2_b(i,j,:),3);
    end
end

active_depth_P_b=depth_P_b;
sensitive_depth_P_b=depth_P_b;


%Remove depths for non recorded or no active or non sensitive days and positions
for i=1:size(depth_P_b,1)
    for j=1:size(depth_P_b,2)
        for k=1:size(depth_P_b,3)
            if depth_P_b(i,j,k)==0 || depth_P_b(i,j,k)==-1 || depth_P_b(i,j,k)==1
                depth_P_b(i,j,k)=NaN;
                active_depth_P_b(i,j,k)=NaN;
                sensitive_depth_P_b(i,j,k)=NaN;
            end
            if isnan(active_units_P1_b(i,j,k)) && isnan(active_units_P2_b(i,j,k))
                active_depth_P_b(i,j,k)=NaN;
                sensitive_depth_P_b(i,j,k)=NaN;
            end
            if isnan(sensitive_units_P1_b(i,j,k)) && isnan(sensitive_units_P2_b(i,j,k))
                sensitive_depth_P_b(i,j,k)=NaN;
            end
            
        end
    end
end

%Compute mean and std
mean_depth_P_b=nanmean(depth_P_b,3);
std_depth_P_b=nanstd(depth_P_b,0,3);
mean_active_depth_P_b=nanmean(active_depth_P_b,3);
std_active_depth_P_b=nanstd(active_depth_P_b,0,3);
mean_sensitive_depth_P_b=nanmean(sensitive_depth_P_b,3);
std_sensitive_depth_P_b=nanstd(sensitive_depth_P_b,0,3);

cd ('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Electrophy_Analysis/Grid')

%% Scooter

cd('/Volumes/buschman/Projects/Learning_Attentional_Templates/Data/ElectrodeTrackingScooter');

for nday=1:size(dates_list_scooter,2)
    
    load(['LearnAttTemplate_scooter_' dates_list_scooter{1,nday} '_ElectrodeTracks_Lowered.mat'],'elec_track')
    
    [   units_P1_s, units_P2_s,...
        active_units_P1_s, active_units_P2_s,...
        sensitive_units_P1_s,sensitive_units_P2_s] = make_sensitive_parietal_grid( elec_track, nday, ...
        active_channels_parietal_scooter{1,nday}, sensitive_channels_parietal_scooter{1,nday}, ...
        units_P1_s, units_P2_s,...
        active_units_P1_s, active_units_P2_s,...
        sensitive_units_P1_s, sensitive_units_P2_s);
    
    [ depth_P_s] = compute_depth_parietal( elec_track, nday, depth_P_s);
    
    clear elec_track
    
end

%Compute total active channels for each position on grid
for i=1:6
    for j=1:8
        total_active_units_P_s(i,j)=nansum(active_units_P1_s(i,j,:),3)+nansum(active_units_P2_s(i,j,:),3);
        total_sensitive_units_P_s(i,j)=nansum(sensitive_units_P1_s(i,j,:),3)+nansum(sensitive_units_P2_s(i,j,:),3);
    end
end

active_depth_P_s=depth_P_s;
sensitive_depth_P_s=depth_P_s;


%Remove depths for non recorded or no active or non sensitive days and positions
for i=1:size(depth_P_s,1)
    for j=1:size(depth_P_s,2)
        for k=1:size(depth_P_s,3)
            if depth_P_s(i,j,k)==0  || depth_P_s(i,j,k)==-1 || depth_P_s(i,j,k)==1
                depth_P_s(i,j,k)=NaN;
                active_depth_P_s(i,j,k)=NaN;
                sensitive_depth_P_s(i,j,k)=NaN;
            end
            if isnan(active_units_P1_s(i,j,k)) && isnan(active_units_P2_s(i,j,k))
                active_depth_P_s(i,j,k)=NaN;
                sensitive_depth_P_s(i,j,k)=NaN;
            end
            if isnan(sensitive_units_P1_s(i,j,k)) && isnan(sensitive_units_P2_s(i,j,k))
                sensitive_depth_P_s(i,j,k)=NaN;
            end
            
        end
    end
end

%Compute mean and std
mean_depth_P_s=nanmean(depth_P_s,3);
std_depth_P_s=nanstd(depth_P_s,0,3);
mean_active_depth_P_s=nanmean(active_depth_P_s,3);
std_active_depth_P_s=nanstd(active_depth_P_s,0,3);
mean_sensitive_depth_P_s=nanmean(sensitive_depth_P_s,3);
std_sensitive_depth_P_s=nanstd(sensitive_depth_P_s,0,3);

cd ('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Electrophy_Analysis/Grid')


%% Brain region attribution


%Create map with colours correponding to regions

for i=1:size(LIP_b,1)
    Parietal_area_b(LIP_b(i,1),LIP_b(i,2))=2;
end

for i=1:size(LIP_s,1)
    Parietal_area_s(LIP_s(i,1),LIP_s(i,2))=2;
end

%% Summary plot

if PLOT==1
    
    figure;
    subplot(4,2,1);
    imagesc(total_active_units_P_s);
    title('Scooter parietal grid activity');
    colorbar;
    caxis([0 20]);
    subplot(4,2,2);
    imagesc(total_active_units_P_b);
    title('Beaker parietal grid activity');
    colorbar;
    caxis([0 20]);
    
    subplot(4,2,3);
    imagesc(total_sensitive_units_P_s);
    title('Scooter parietal grid selectivity');
    colorbar;
    caxis([0 20]);
    subplot(4,2,4);
    imagesc(total_sensitive_units_P_b);
    title('Beaker parietal grid selectivity');
    colorbar;
    caxis([0 20]);
    
    subplot(4,2,5);
    imagesc(-mean_active_depth_P_s);
    title('Scooter parietal active grid mean depth');
    colorbar;
    caxis([0 200]);
    subplot(4,2,6);
    imagesc(-mean_active_depth_P_b);
    title('Beaker parietal active grid mean depth');
    colorbar;
    caxis([0 200]);
    
    subplot(4,2,7);
    imagesc(-mean_sensitive_depth_P_s);
    title('Scooter parietal selective grid mean depth');
    colorbar;
    caxis([0 200]);
    subplot(4,2,8);
    imagesc(-mean_sensitive_depth_P_b);
    title('Beaker parietal selective grid mean depth');
    colorbar;
    caxis([0 200]);
    
end













