function [units_frontal1_s, buffer_frontal1_s ] = reshape_unit_grid( elec_track, grid, units_frontal1_s, n_day)
%check that the GridShape matrices and the template
    %have the same dimension, if not add rows/cols of zeros where
    %needed
        
    %Row first
    i=1;
    if sum(elec_track.GridShape{grid,1}(1,:))==0 && sum(units_frontal1_s(i,:,n_day))>0 %check if you need to move something row-wise
        while sum(elec_track.GridShape{1,1}(i,:))==0 %row first
            buffer_frontal1_s=units_frontal1_s(:,:,n_day);
            units_frontal1_s(i,:,n_day)=buffer_frontal1_s(end,:);
            for j=i:size(units_frontal1_s,1)-1
                units_frontal1_s(j,:,n_day)=buffer_frontal1_s(j+1,:,n_day);
            end
            i=i+1;
        end
    end
    
    %Col second
    i=1;
    if sum(elec_track.GridShape{1,1}(:,1))==0 && sum(units_frontal1_s(:,i,n_day))>0  %check if you need to move something col-wise
        while sum(elec_track.GridShape{1,1}(:,i))==0 %col second
            buffer_frontal1_s=units_frontal1_s(:,:,n_day);
            units_frontal1_s(i,:,n_day)=buffer_frontal1_s(:,end);
            for j=i:size(units_frontal1_s,2)-1
                units_frontal1_s(:,j,n_day)=buffer_frontal1_s(:,j+1,n_day);
            end
            i=i+1;
        end
    end
    
end

