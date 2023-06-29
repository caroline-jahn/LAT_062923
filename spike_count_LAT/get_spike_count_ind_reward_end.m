function [spkcount] = get_spike_count_ind_reward_end(Ts, cutTrials, cellnb, incode,  ind, win)
%only difference is that we have a different window for each trial

%% get the time serie
spktimes=Ts{cellnb};

StartTime=spktimes(1);
EndTime=spktimes(end);



%% load global variables and corresponding data from file


%% Prepares intermediate variables for plot

trialcnt=length(cutTrials);

evt=NaN(trialcnt,1);
win=win./1000;

for ti=1:trialcnt
    %if correct trial
    if ind(ti)>0 && ~isnan(cutTrials(ti).Timing.NS4.TimeStampSec(1))
        if cutTrials(ti).Timing.NS4.TimeStampSec(1)>StartTime && cutTrials(ti).Timing.NS4.TimeStampSec(4)<EndTime %only for correct trials and cell recorded during the whole trial
        %incode=event time
        if incode==1
            evt(ti)=cutTrials(ti).Timing.NS4.TimeStampSec(1);
        elseif  incode==2
            evt(ti)=cutTrials(ti).Timing.NS4.TimeStampSec(2);
        elseif incode==3
            evt(ti)=cutTrials(ti).Timing.NEV.RESPONSE_ON_Sec;
        elseif incode==4
            evt(ti)=cutTrials(ti).Timing.NS4.TimeStampSec(3);
        end
        
        end
    end
end


spkcount=NaN(sum(~isnan(ind)),1);


for i=1:trialcnt
    if ~isnan(ind(i))
       spkcount(ind(i),1)=length(find(spktimes>win(i,1)+evt(i) & spktimes<win(i,2)+evt(i)))/(win(i,2)-win(i,1));
    end
end

for i=1:size(spkcount,1)-500
    if sum(spkcount(i:i+499))==0
        spkcount(i:end)=NaN;
        break
    end
end

% y=spkcount;
% x=1:length(ind);
% x=zscore(x);
% beta(:)=glmfit(x,y);
% y(:)=y(:)-beta(2).*x(:);
% y(:)=(y-nanmean(y))./nanstd(y);
% 
% spkcount=y;
% 






end

