function [tempTrials] = add_info(tempTrials, trials)

%%
% Add necessary time info to trials
% Inp[ut is aster NS4 correction but uses NEV timing 
%For instance add response info by combining NEV and RT info
% Add comments %%% Do sanity checks %%%

%%

%%
for ti=1:length(trials)
    
        if length(tempTrials(ti).Timing.NEV.TimeStampSec)>=4       
            tempTrials(ti).Timing.NEV.RESPONSE_ON_Sec = trials(ti).Behavior.ReactionTime + tempTrials(ti).Timing.NEV.TimeStampSec(4); 
        end
    
end

