function [ units_frontal1,units_frontal2, units_parietal1, units_parietal2,...
    active_units_frontal1,active_units_frontal2,...
    active_units_parietal1, active_units_parietal2,...
    FEF_active_channels, PFC_active_channels, LIP_active_channels] = make_active_grid( elec_track, n_day, ...
    active_channels_list_frontal, active_channels_list_parietal, ...
    units_frontal1, units_frontal2, units_parietal1,units_parietal2,...
    active_units_frontal1, active_units_frontal2, active_units_parietal1, active_units_parietal2,...
    FEF, PFC, LIP)

%Initialise lists
FEF_active_channels=[];
PFC_active_channels=[];
LIP_active_channels=[];


%Create the grid for each session

    %First fill GridElectrodChannels in the channel matrices on the top
    %left corner
    for i=1:size(elec_track.GridElectrodeChannels{1,1},1)
        for j=1:size(elec_track.GridElectrodeChannels{1,1},2)
            units_frontal1(i,j,n_day)=elec_track.GridElectrodeChannels{1,1}(i,j,1);
            units_frontal2(i,j,n_day)=elec_track.GridElectrodeChannels{1,1}(i,j,2);
        end
    end
    for i=1:size(elec_track.GridElectrodeChannels{1,2},1)
        for j=1:size(elec_track.GridElectrodeChannels{1,2},2)
            units_parietal1(i,j,n_day)=elec_track.GridElectrodeChannels{1,2}(i,j,1);
            units_parietal2(i,j,n_day)=elec_track.GridElectrodeChannels{1,2}(i,j,2);
        end
    end
    
    %Second check that the GridShape matrices and the template
    %have the same dimension, if not add rows/cols of zeros where
    %needed
    %Frontal grid
[units_frontal1] = reshape_unit_grid( elec_track, 1, units_frontal1, n_day);
[units_frontal2] = reshape_unit_grid( elec_track, 1, units_frontal2, n_day);
    %Parietal grid
[units_parietal1] = reshape_unit_grid( elec_track, 2, units_parietal1, n_day);
[units_parietal2] = reshape_unit_grid( elec_track, 2, units_parietal2, n_day);
    
    %Third remove cells that are not on the Shape matrix (from history if
    %not starting from a parameter file??)
    units_frontal1(:,:,n_day)=units_frontal1(:,:,n_day).*elec_track.GridShape{1,1}(:,:);
    units_frontal2(:,:,n_day)=units_frontal2(:,:,n_day).*elec_track.GridShape{1,1}(:,:);
    units_parietal1(:,:,n_day)=units_parietal1(:,:,n_day).*elec_track.GridShape{2,1}(:,:);
    units_parietal2(:,:,n_day)=units_parietal2(:,:,n_day).*elec_track.GridShape{2,1}(:,:);
    
    %Fourth find the active channels and put 1 in corresponding cell
    
    for i=1:size(units_frontal1,1)
        for j=1:size(units_frontal1,2)
            
            for k=1:size(active_channels_list_frontal{1,n_day}{1,1},2)
                if units_frontal1(i,j,n_day)==active_channels_list_frontal{1,n_day}{1,1}(1,k)
                    active_units_frontal1(i,j,n_day)=1;
                elseif units_frontal2(i,j,n_day)==active_channels_list_frontal{1,n_day}{1,1}(1,k)
                    active_units_frontal2(i,j,n_day)=1;
                end
            end
            
            for k=1:size(active_channels_list_parietal{1,n_day}{1,1},2)
                if units_parietal1(i,j,n_day)==active_channels_list_parietal{1,n_day}{1,1}(1,k)
                    active_units_parietal1(i,j,n_day)=1;
                elseif units_parietal2(i,j,n_day)==active_channels_list_parietal{1,n_day}{1,1}(1,k)
                    active_units_parietal2(i,j,n_day)=1;
                end
            end
            
        end
    end
    
    % Fifth, create list of active channel for each region
    for i=1:size(FEF,1)
        if active_units_frontal1(FEF(i,1),FEF(i,2),n_day)==1
            FEF_active_channels=[FEF_active_channels, units_frontal1(FEF(i,1),FEF(i,2),n_day)];
        end
        if active_units_frontal2(FEF(i,1),FEF(i,2),n_day)==1
            FEF_active_channels=[FEF_active_channels, units_frontal2(FEF(i,1),FEF(i,2),n_day)];
        end
    end
     for i=1:size(PFC,1)
        if active_units_frontal1(PFC(i,1),PFC(i,2),n_day)==1
            PFC_active_channels=[PFC_active_channels, units_frontal1(PFC(i,1),PFC(i,2),n_day)];
        end
        if active_units_frontal2(PFC(i,1),PFC(i,2),n_day)==1
            PFC_active_channels=[PFC_active_channels, units_frontal2(PFC(i,1),PFC(i,2),n_day)];
        end
     end

     for i=1:size(LIP,1)
        if active_units_parietal1(LIP(i,1),LIP(i,2),n_day)==1
            LIP_active_channels=[LIP_active_channels, units_parietal1(LIP(i,1),LIP(i,2),n_day)];
        end
        if active_units_parietal2(LIP(i,1),LIP(i,2),n_day)==1
            LIP_active_channels=[LIP_active_channels, units_parietal2(LIP(i,1),LIP(i,2),n_day)];
        end
    end

    
end


