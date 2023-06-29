%% ChannelsActivityToGrid
%From active channels on each days give you the correpondind grid location
%Warning: grids' orientations have to be the same over the deys of
%recordings OR need to had a rotation (not coded here)

addpath('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Electrophy_Analysis/Grid');
addpath('/Volumes/buschman/Projects/Learning_Attentional_Templates/Data/');

%% OPTIONS
PLOT=1;

%% Sessions for Beaker

dates_list_beaker={'101718','101918','102518','102918','103018','110618','110918','111018','111518'}; %'102618', %,'111418','111518' %,'111418'

active_channels_list_frontal_beaker={{[1 2 4 5 7 10 12 13 16 18 19 22 23 24 26 28 31]},...
    {[1 4 6 7 10 12 13 16 18 19 22 23 24 25 28 29 32]}, ...
    {[2 3 8 10 11 13 14 16 17 18 19 20 23 24 25 26 28 30 32]},... %{[2 4 5 8 9 10 11 12 13 14 16 17 18 19 20 22 23 24 25 26 29 31 32]},...
    {[2 3 5 6 7 8 9 10 12 13 15 16 18 19 21 22 23 26 27 29 30 32]},...
    {[1 3 4 5 7 8 10 12 15 16 18 20 21 22 23 24 32]},...
    {[2 4 6 7 9 10 13 14 16 18 20 21 22 23 24 26 28 30 31 34 37 40 41 49]},...
    {[1 2 4 5 6 10 11 12 13 14 15 16 17 18 21 23 24 25 26 27 28 30 32 33 34 35 37 38 39 40 41 42]},...
    {[1 2 4 5 6 10 13 15 16 17 18 22 23 24 25 26 27 30 32 33 34 35 37 38 39 40 41 42]},... %{[1 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 20 21 23 26 27 28 30 31 33 34 36 37 39 40]},...
    {[3 4 5 6 7 12 13 14 16 17 18 20 22 23 24 26 27 29 30 33 34 36 39 40]}};

active_channels_list_parietal_beaker={{[34 35 37 38 43 44 47 48]},...
    {[34 35 37 38 41 42 43]}, ...
    {[33 37 38 41 42 43 44 45 46 47 48]},... %{[33 34 37 39 42 43 44 46 47 48]},...
    {[33 34 39 40 41 42 43 44 45 46 47 48]},...
    {[34 38 40 41 42 43]},...
    {[51 52 54 55 56 58 61 62 64 67 68]},...
    {[65 66 67 70 71 74 75 76 80 81 84]},...
    {[65 66 68 69 70 71 74 75 76 80 81]},... %{[65 69 71 72 73 74 75 76 77 78 79 80 82 83]},...
    {[65 69 71 72 73 74 79 80 81 83]}};


%% Sessions for Scooter

dates_list_scooter={'110218','110318','110718','110818','111218','111318','111618','111718','111918','112618','113018'}; %%

active_channels_list_frontal_scooter={{[1 2 4 5 6 7 8 11 12 13 14 15 28 49]},...
    {[1 2 4 6 7 8 9 10 11 12 13 18 19 20 21 22 23 25 26 28 30 31 32 33 34 35 36 37]},...
    {[1 4 5 6 8 9 10 11 12 13 15 16 17 20 21 22 24 25 26 31 33 34 35 36]},...
    {[2 4 5 6 7 9 10 11 12 13 15 16 17 20 21 22 24 25 26 31 34 37 38 40]},...
    {[3 5 6 7 8 9 10 11 12 13 14 16 17 19 21 24 25 27 28 29 32 35 36]},...
    {[1 3 5 6 7 8 9 10 11 12 13 14 16 17 19 21 22 23 24 25 26 27 28 29 32 35 36],}...
    {[3 5 6 7 8 10 11 12 14 15 16 17 19 20 21 22 24 25 27 31 33 35 36 38],}...
    {[1 3 5 6 7 8 9 10 11 12 14 15 16 19 20 21 22 24 25 26 27 31 33 34 38 39]},...
    {[1 3 5 6 7 8 9 10 11 12 13 14 15 19 20 21 22 23 26 27 31 35 38 39]},...
    {[2 3 5 7 8 9 11 12 15 18 28 29 31]}, ...
    {[1 2 7 8 9 11 12 13 15 16 22 28 31 36 40]}};

active_channels_list_parietal_scooter={{[54 55 66 67 68 69 70 71 72 73 74 75 76 79 80]},...
    {[50 53 65 66 67 68 69 70 71 74 75 76 78 79 80]}, ...
    {[65 67 68 69 70 81 82 85 87 89 90 92 95]},...
    {[70 85 86 87 88 89 91]},...
    {[65 66 68 72 81 82 85 86 88 89 90 91 92 94 96]},...
    {[66 67 68 71 72 82 85 86 87 88 89 90 91 92 94 96]},...
    {[65 69 70 72 75 81 82 84 85 86 87 89 90 92 93 94],}...
    {[65 69 72 75 76 81 82 83 84 85 86 88 90 91 94]},...
    {[66 69 70 76 82 83 84 86 87 89 91 95 96]},...
    {[100 102 104 105 106 112 124]}, ...
    {[97 99 100 105 111 117 123 124]}};


%% Grids
%Check with Matt's scipt to localise 

%Parital grid
LIP_b=[3 6; 4 5; 4 6; 4 7; 5 4; 5 5; 5 6; 3 5; 4 4]; %last two might be MIP but check with reponse during delay
LIP_s=[1 6; 2 4; 2 5; 2 6; 3 2; 3 3; 3 4; 3 5; 4 1; 4 2; 4 3];

%OtherP_b=[3 5; 4 4];

%Frontal grid
FEF_b=[2 5; 3 5; 3 6; 4 6];
FEF_s=[4 3; 4 4; 4 5; 3 5; 2 5];

PFC_b=[3 2; 4 2; 2 3; 3 3; 2 4; 4 3; 3 4; 4 4;  4 5];
PFC_s=[1 3; 2 3; 3 3; 3 2; 4 2; 3 4; 3 1; 4 1; 2 4];

OtherF_b=[5 3; 5 4; 5 5; 3 7; 4 7; 2 6; 5 6];
OtherF_s=[5 4; 5 5; 5 3; 3 6; 4 6];


%% Inititalise variables

%Beaker
units_P1_b=NaN(6,8,size(dates_list_beaker,2));
units_P2_b=NaN(6,8,size(dates_list_beaker,2));

units_F1_b=NaN(6,8,size(dates_list_beaker,2));
units_F2_b=NaN(6,8,size(dates_list_beaker,2));

active_units_P1_b=NaN(6,8,size(dates_list_beaker,2));
active_units_P2_b=NaN(6,8,size(dates_list_beaker,2));

active_units_F1_b=NaN(6,8,size(dates_list_beaker,2));
active_units_F2_b=NaN(6,8,size(dates_list_beaker,2));

total_active_units_F_b=NaN(6,8);
total_active_units_P_b=NaN(6,8);

depth_F_b=zeros(6,8,size(dates_list_beaker,2));
depth_P_b=zeros(6,8,size(dates_list_beaker,2));

%Scooter
units_P1_s=NaN(6,8,size(dates_list_scooter,2)); %channel matrix (electrod number, 3rd dim = 2 electrodes/position)
units_P2_s=NaN(6,8,size(dates_list_scooter,2));

units_F1_s=NaN(6,8,size(dates_list_scooter,2));
units_F2_s=NaN(6,8,size(dates_list_scooter,2));

active_units_P1_s=NaN(6,8,size(dates_list_scooter,2));  %activity channel (1 in active cell corresponding to active channel)
active_units_P2_s=NaN(6,8,size(dates_list_scooter,2));

active_units_F1_s=NaN(6,8,size(dates_list_scooter,2));
active_units_F2_s=NaN(6,8,size(dates_list_scooter,2));

total_active_units_F_s=NaN(6,8); %sum over the 2 positions and over days of recording
total_active_units_P_s=NaN(6,8);

depth_F_s=zeros(6,8,size(dates_list_scooter,2));
depth_P_s=zeros(6,8,size(dates_list_scooter,2));

%Active channels count
LIP=zeros(1,2);
% OtherP=zeros(1,2);

FEF=zeros(1,2);
PFC=zeros(1,2);
OtherF=zeros(1,2);

%Color map of brain areas on grid
Parietal_area_b=NaN(6,8);
Parietal_area_s=NaN(6,8);

Frontal_area_b=NaN(6,8);
Frontal_area_s=NaN(6,8);


%% Beaker

cd('/Volumes/buschman/Projects/Learning_Attentional_Templates/Data/ElectrodeTracking');

for nday=1:size(dates_list_beaker,2)
    
    load(['LearnAttTemplate_Beaker_' dates_list_beaker{1,nday} '_ElectrodeTracks_Lowered.mat'],'elec_track')
    
    [ units_F1_b,units_F2_b, units_P1_b, units_P2_b,...
        active_units_F1_b, active_units_F2_b,...
        active_units_P1_b,active_units_P2_b,...
        FEF_active_channels_beaker{1,nday}, PFC_active_channels_beaker{1,nday}, LIP_active_channels_beaker{1,nday}] = make_active_grid( elec_track, nday, ...
        active_channels_list_frontal_beaker, active_channels_list_parietal_beaker, ...
        units_F1_b, units_F2_b, units_P1_b, units_P2_b,...
        active_units_F1_b, active_units_F2_b,...
        active_units_P1_b, active_units_P2_b,...
        FEF_b, PFC_b, LIP_b);
    
    [ depth_P_b, depth_F_b ] = compute_depth( elec_track, nday, depth_P_b, depth_F_b);
    
    clear elec_track
    
end

%Compute total active channels for each position on grid
for i=1:6
    for j=1:8
        total_active_units_F_b(i,j)=nansum(active_units_F1_b(i,j,:),3)+nansum(active_units_F2_b(i,j,:),3);
        total_active_units_P_b(i,j)=nansum(active_units_P1_b(i,j,:),3)+nansum(active_units_P2_b(i,j,:),3);
    end
end

%Remove depths for non recorded days and positions
for i=1:size(depth_F_b,1)
    for j=1:size(depth_F_b,2)
        for k=1:size(depth_F_b,3)
            if depth_F_b(i,j,k)==0
                depth_F_b(i,j,k)=NaN;
            end
            if depth_P_b(i,j,k)==0
                depth_P_b(i,j,k)=NaN;
            end
        end
    end
end

%Compute mean and std
mean_depth_F_b=nanmean(depth_F_b,3);
std_depth_F_b=nanstd(depth_F_b,0,3);
mean_depth_P_b=nanmean(depth_P_b,3);
std_depth_P_b=nanstd(depth_P_b,0,3);

%Make list of active units for each brain region



cd ('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Electrophy_Analysis/Grid')

%% Scooter

cd('/Volumes/buschman/Projects/Learning_Attentional_Templates/Data/ElectrodeTrackingScooter');

for nday=1:size(dates_list_scooter,2)
    
    load(['LearnAttTemplate_Scooter_' dates_list_scooter{1,nday} '_ElectrodeTracks_Lowered.mat'],'elec_track')
    
    [units_F1_s,units_F2_s, units_P1_s, units_P2_s,...
        active_units_F1_s, active_units_F2_s,...
        active_units_P1_s,active_units_P2_s,...
        FEF_active_channels_scooter{1,nday}, PFC_active_channels_scooter{1,nday}, LIP_active_channels_scooter{1,nday}] = make_active_grid( elec_track, nday, ...
        active_channels_list_frontal_scooter, active_channels_list_parietal_scooter, ...
        units_F1_s, units_F2_s, units_P1_s, units_P2_s,...
        active_units_F1_s, active_units_F2_s,...
        active_units_P1_s, active_units_P2_s,...
        FEF_s, PFC_s, LIP_s);
    
    [ depth_P_s, depth_F_s ] = compute_depth( elec_track, nday, depth_P_s, depth_F_s);
    
    clear elec_track
    
end

for i=1:6
    for j=1:8
        total_active_units_F_s(i,j)=nansum(active_units_F1_s(i,j,:),3)+nansum(active_units_F2_s(i,j,:),3);
        total_active_units_P_s(i,j)=nansum(active_units_P1_s(i,j,:),3)+nansum(active_units_P2_s(i,j,:),3);
    end
end

cd ('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Electrophy_Analysis')

%Remove depths for non recorded days and positions
for i=1:size(depth_F_s,1)
    for j=1:size(depth_F_s,2)
        for k=1:size(depth_F_s,3)
            if depth_F_s(i,j,k)==0
                depth_F_s(i,j,k)=NaN;
            end
            if depth_P_s(i,j,k)==0
                depth_P_s(i,j,k)=NaN;
            end
        end
    end
end

%Compute mean and std
mean_depth_F_s=nanmean(depth_F_s,3);
std_depth_F_s=nanstd(depth_F_s,0,3);
mean_depth_P_s=nanmean(depth_P_s,3);
std_depth_P_s=nanstd(depth_P_s,0,3);

%% Brain region attribution

%WARNING: mirror transformation compared to png images

%Parietal
for i=1:size(LIP_b,1)
    LIP(1,1)=LIP(1,1)+total_active_units_P_b(LIP_b(i,1),LIP_b(i,2));
end
for i=1:size(LIP_s,1)
    LIP(1,2)=LIP(1,2)+total_active_units_P_s(LIP_s(i,1),LIP_s(i,2));
end

% for i=1:size(OtherP_b,1)
%     OtherP(1,1)=OtherP(1,1)+total_active_units_P_b(OtherP_b(i,1),OtherP_b(i,2));
% end

%Frontal
for i=1:size(FEF_b,1)
    FEF(1,1)=FEF(1,1)+total_active_units_F_b(FEF_b(i,1),FEF_b(i,2));
end
for i=1:size(FEF_s,1)
    FEF(1,2)=FEF(1,2)+total_active_units_F_s(FEF_s(i,1),FEF_s(i,2));
end

for i=1:size(PFC_b,1)
    PFC(1,1)=PFC(1,1)+total_active_units_F_b(PFC_b(i,1),PFC_b(i,2));
end
for i=1:size(PFC_s,1)
    PFC(1,2)=PFC(1,2)+total_active_units_F_s(PFC_s(i,1),PFC_s(i,2));
end

for i=1:size(OtherF_b,1)
    OtherF(1,1)=OtherF(1,1)+total_active_units_F_b(OtherF_b(i,1),OtherF_b(i,2));
end
for i=1:size(OtherF_s,1)
    OtherF(1,2)=OtherF(1,1)+total_active_units_F_s(OtherF_s(i,1),OtherF_s(i,2));
end

%Create map with colours correponding to regions

for i=1:size(LIP_b,1)
    Parietal_area_b(LIP_b(i,1),LIP_b(i,2))=2;
end
% for i=1:size(OtherP_b,1)
%     Parietal_area_b(OtherP_b(i,1),OtherP_b(i,2))=1;
% end
for i=1:size(LIP_s,1)
    Parietal_area_s(LIP_s(i,1),LIP_s(i,2))=2;
end

for i=1:size(FEF_b,1)
    Frontal_area_b(FEF_b(i,1),FEF_b(i,2))=3;
end
for i=1:size(FEF_s,1)
    Frontal_area_s(FEF_s(i,1),FEF_s(i,2))=3;
end
for i=1:size(PFC_b,1)
    Frontal_area_b(PFC_b(i,1),PFC_b(i,2))=2;
end
for i=1:size(PFC_s,1)
    Frontal_area_s(PFC_s(i,1),PFC_s(i,2))=2;
end
for i=1:size(OtherF_b,1)
    Frontal_area_b(OtherF_b(i,1),OtherF_b(i,2))=1;
end
for i=1:size(OtherF_s,1)
    Frontal_area_s(OtherF_s(i,1),OtherF_s(i,2))=1;
end

%% Summary plot

if PLOT==1
    
    figure;
    subplot(4,2,1);
    imagesc(Frontal_area_b);
    title('Beaker frontal grid areas');
    colorbar;
    caxis([0 3]);
    subplot(4,2,2);
    imagesc(Parietal_area_b);
    title('Beaker parietal grid areas');
    colorbar;
    caxis([0 2]);
    subplot(4,2,3);
    imagesc(total_active_units_F_b);
    title('Beaker frontal grid activity');
    colorbar;
    caxis([0 20]);
    subplot(4,2,4);
    imagesc(total_active_units_P_b);
    title('Beaker parietal grid activity');
    colorbar;
    caxis([0 20]);
    subplot(4,2,5);
    imagesc(-mean_depth_F_b);
    title('Beaker frontal grid mean depth');
    colorbar;
    caxis([0 200]);
    subplot(4,2,6);
    imagesc(-mean_depth_P_b);
    title('Beaker parietal grid mean depth');
    colorbar;
    caxis([0 200]);
    subplot(4,2,7);
    imagesc(std_depth_F_b);
    title('Beaker frontal grid std depth');
    colorbar;
    caxis([0 5]);
    subplot(4,2,8);
    imagesc(std_depth_P_b);
    title('Beaker parietal grid std depth');
    colorbar;
    caxis([0 5]);
    
    figure;
    subplot(4,2,1);
    imagesc(Frontal_area_s);
    title('Scooter frontal grid areas');
    colorbar;
    caxis([0 3]);
    subplot(4,2,2);
    imagesc(Parietal_area_s);
    title('Scooter parietal grid areas');
    colorbar;
    caxis([0 2]);
    subplot(4,2,3);
    imagesc(total_active_units_F_s);
    title('Scooter frontal grid activity');
    colorbar;
    caxis([0 20]);
    subplot(4,2,4);
    imagesc(total_active_units_P_s);
    title('Scooter parietal grid activity');
    colorbar;
    caxis([0 20]);
    subplot(4,2,5);
    imagesc(-mean_depth_F_s);
    title('Scooter frontal grid mean depth');
    colorbar;
    caxis([0 200]);
    subplot(4,2,6);
    imagesc(-mean_depth_P_s);
    title('Scooter parietal grid mean depth');
    colorbar;
    caxis([0 200]);
    subplot(4,2,7);
    imagesc(std_depth_F_s);
    title('Scooter frontal grid std depth');
    colorbar;
    caxis([0 5]);
    subplot(4,2,8);
    imagesc(std_depth_P_s);
    title('Scooter parietal grid std depth');
    colorbar;
    caxis([0 5]);
    
end

%% Save active channels for ROI is ns6directory

dates_list_beaker_for_ns6={'181017','181019','181025','181029','181030','181106','181109','181110','181115'}; %'102618', %,'111418','111518' %,'111418'
dates_list_scooter_for_ns6={'181102','181103','181107','181108','181112','181113','181116','181117','181119','181126','181130'}; %%



    destpath='/Volumes/buschman/Projects/Learning_Attentional_Templates/Data/General';
    load(fullfile(destpath,'NS6Directory_sem.mat'));
    
    
    %Give active channel list for each regions to each date
    %Beaker
    for nday=1:size(dates_list_beaker_for_ns6,2)
        for i=1:length(ns6directory)
            if strcmp(dates_list_beaker_for_ns6(nday),ns6directory(i).FolderStem(end-5:end))==1 && strcmp(ns6directory(i).Subject,'Beaker')==1 %same date and same subject => should work for both delsaccade and ecxploreexploit
                ns6directory(i).FEF=FEF_active_channels_beaker{1,nday};
                ns6directory(i).PFC=PFC_active_channels_beaker{1,nday};
                ns6directory(i).LIP=LIP_active_channels_beaker{1,nday};
            end
        end
    end
                
    %Scooter
    for nday=1:size(dates_list_scooter_for_ns6,2)
        for i=1:length(ns6directory)
            if strcmp(dates_list_scooter_for_ns6(nday),ns6directory(i).FolderStem(end-5:end))==1 && strcmp(ns6directory(i).Subject,'Scooter')==1 %same date and same subject => should work for both delsaccade and ecxploreexploit
                ns6directory(i).FEF=FEF_active_channels_scooter{1,nday};
                ns6directory(i).PFC=PFC_active_channels_scooter{1,nday};
                ns6directory(i).LIP=LIP_active_channels_scooter{1,nday};
            end
        end
    end

save(fullfile(destpath, 'NS6Directory_AC_restricted_FEF'), 'ns6directory');




















