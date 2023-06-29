function [ depth_parietal, depth_frontal ] = compute_depth( elec_track, n_day, depth_parietal, depth_frontal )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

        %Fifth and unrelated, compute the depth of the recording for each channel and each
    %day byt adding the values in Movements
   for i=1:size(elec_track.GridElectrodeChannels{1,1},1)
        for j=1:size(elec_track.GridElectrodeChannels{1,1},2)
            for k=1:size(elec_track.Movements,1)
                if elec_track.GridShape{1,1}(i,j)==1 && elec_track.Movements(k,1)==1 && elec_track.Movements(k,2)==i && elec_track.Movements(k,3)==j
                    depth_frontal(i,j,n_day)=depth_frontal(i,j)+elec_track.Movements(k,4);
                end
            end
        end
   end
    
      for i=1:size(elec_track.GridElectrodeChannels{1,2},1)
        for j=1:size(elec_track.GridElectrodeChannels{1,2},2)
            for k=1:size(elec_track.Movements,1)
                if elec_track.GridShape{2,1}(i,j)==1 && elec_track.Movements(k,1)==2 && elec_track.Movements(k,2)==i && elec_track.Movements(k,3)==j
                    depth_parietal(i,j,n_day)=depth_parietal(i,j)+elec_track.Movements(k,4);
                end
            end
        end
      end
    

end

