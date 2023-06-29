function [units_parietal1, units_parietal2,...
    active_units_parietal1, active_units_parietal2,...
    sensitive_units_parietal1, sensitive_units_parietal2] = make_sensitive_parietal_grid( elec_track, n_day, ...
    active_channels_list_parietal, sensitive_channels_list_parietal, ...
    units_parietal1, units_parietal2,...
    active_units_parietal1, active_units_parietal2,...
    sensitive_units_parietal1, sensitive_units_parietal2)

%Create the grid for each session

%First fill GridElectrodChannels in the channel matrices on the top
%left corner
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
%Parietal grid
[units_parietal1] = reshape_unit_grid( elec_track, 2, units_parietal1, n_day);
[units_parietal2] = reshape_unit_grid( elec_track, 2, units_parietal2, n_day);

%Third remove cells that are not on the Shape matrix (from history if
%not starting from a parameter file??)
units_parietal1(:,:,n_day)=units_parietal1(:,:,n_day).*elec_track.GridShape{2,1}(:,:);
units_parietal2(:,:,n_day)=units_parietal2(:,:,n_day).*elec_track.GridShape{2,1}(:,:);

%Fourth find the active channels and put 1 in corresponding cell

for i=1:size(units_parietal1,1)
    for j=1:size(units_parietal1,2)
        %active channels
        for k=1:size(active_channels_list_parietal,2)
            if units_parietal1(i,j,n_day)==active_channels_list_parietal(1,k)
                active_units_parietal1(i,j,n_day)=1;
            elseif units_parietal2(i,j,n_day)==active_channels_list_parietal(1,k)
                active_units_parietal2(i,j,n_day)=1;
            end
        end
        
        for k=1:size(sensitive_channels_list_parietal,2)
            if units_parietal1(i,j,n_day)==sensitive_channels_list_parietal(1,k)
                sensitive_units_parietal1(i,j,n_day)=1;
            elseif units_parietal2(i,j,n_day)==sensitive_channels_list_parietal(1,k)
                sensitive_units_parietal2(i,j,n_day)=1;
            end
        end
        
        
    end
end

end


