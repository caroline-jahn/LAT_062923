function [tempTrialsNew, tempSpecs] = cutNS4Trial(NS4, tempTrialsOld, tempSpecs, arrayID)

if iscell(NS4.Data)
    bitcode = [];
    for i=1:length(NS4.Data)
        bitcode = [bitcode NS4.Data{i}(4,:)];
    end
    timeStampSec = zeros(1,length(bitcode));
    for i=1:length(NS4.Data)
        inds = sum(NS4.MetaTags.DataPoints(1:i-1))+(1:NS4.MetaTags.DataPoints(i));
        timeStampSec(inds) = NS4.MetaTags.Timestamp(i)/NS4.MetaTags.TimeRes+((1:NS4.MetaTags.DataPoints(i))/NS4.MetaTags.SamplingFreq);
    end
else
    bitcode = NS4.Data(4,:);
    timeStampSec = (1:NS4.MetaTags.DataPoints)/NS4.MetaTags.SamplingFreq;
end
timeStamp = timeStampSec*NS4.MetaTags.SamplingFreq;

%%
[counts, bins] = hist(double(bitcode),9);
[~, inds] = sort(counts,'descend');
binsoi = sort(bins(inds(1:2)));
binsdiff = abs(diff(binsoi));
downThr = round((binsoi(1)+binsdiff/4)/100)*100;
upThr = round((binsoi(2)-binsdiff/4)/100)*100;
% Arbitrary thresholds for bitcode change detection
% upThr = 8000;
% downThr = 6000;

pos = 0;
state = 0;
xInd = []; %inds of state transition
xStamp = []; %timeStamps of state transition
xStampSec = []; %time of state transition
ySignal = []; %bitcode value at state transition
stateList = []; %transition state [1 means transition up, 0 means transition down]
for i=1:length(bitcode)
    if (~state & bitcode(i) > upThr) | (state & bitcode(i) < downThr),
        pos = pos+1;
        state = ~state;
        xInd(end+1) = i;
        xStamp(end+1) = timeStamp(i);
        xStampSec(end+1) = timeStampSec(i);
        ySignal(end+1) = bitcode(i);
        stateList(end+1) = state;
    end
end

%%
onsetList = zeros(1,length(tempTrialsOld));
timeDiffRange = [-.02 .1]; %consider -20 +100ms window of tolerance for Encode to bitcode change
tempTrialsNew = tempTrialsOld;
for ti=1:length(tempTrialsOld),

fprintf('%d\n',ti);
    tempTrialsNew(ti).Timing = tempTrialsOld(ti).Timing;
    tempTrialsNew(ti).Timing.NS4 = struct('TimeStamp', NaN, 'TimeStampSec', NaN);
    fixOnset = tempTrialsOld(ti).Timing.NEV.TimeStampSec(1); %find FIX_ON time from encodes
    trialEnd = tempTrialsOld(ti).Timing.NEV.TimeStampSec(end); %find trial end time from encodes
    
    if ~isnan(fixOnset) && ~isnan(trialEnd),
        fixInds = find(timeStampSec >= fixOnset+timeDiffRange(1) & timeStampSec <= fixOnset+timeDiffRange(2)); %range of timestamps to look for earliest upstate
        endInds = find(timeStampSec >= trialEnd+timeDiffRange(1) & timeStampSec <= trialEnd+timeDiffRange(2)); %range of timestamps to look for earliest upstate
        
        
        %if whole trial is included in recording
        if ~isempty(fixInds) && ~isempty(endInds),

            tempTrialsNew(ti).Timing.NS4.TimeStamp = zeros(1,length(tempTrialsOld(ti).Timing.NEV.StageSequence));
            tempTrialsNew(ti).Timing.NS4.TimeStampSec = zeros(1,length(tempTrialsOld(ti).Timing.NEV.StageSequence));
     
    
            %%% PK
            
            if arrayID==47 %% Debug
                xInd(7)=fixInds(7);
            end
            

            onsetInd = find(ismember(xInd,fixInds) & stateList,1,'first');
            onsetList(ti) = xInd(onsetInd);
            %assume that all stage changes of trial will follow the
            %fixonset
            range = onsetInd:onsetInd+length(tempTrialsOld(ti).Timing.NEV.StageSequence)-1;
            tempoffset = repmat([-.004 -.0065],1,round(length(range)/2));
            tempoffset = tempoffset(1:length(range));
            tempTrialsNew(ti).Timing.NS4.TimeStampSec = xStampSec(range)-tempoffset;
            tempTrialsNew(ti).Timing.NS4.TimeStamp = xStamp(range)+(tempoffset*NS4.MetaTags.SamplingFreq);
            %%% pk
            %OSCAR
            %Hardcode a -4ms, -40 stamps diff for up moves
            %and -4.5ms, -45 stamps diff for down moves
            %GONZO
            %Hardcode a -4ms, -40 stamps diff for up moves
            %and -6.5ms, -65 stamps diff for down moves
            
        end
    end
end

tempSpecs.Recording.NS4SampleRes = NS4.MetaTags.SamplingFreq;
%%

visualize = 0;
if visualize == 1,
    figure;
    h = axes;
    segsize = 5*10^6;
    c = [.9 .9 .1; .6 .1 .7; .2 .4 .8];
    colormap(c);

    for si=1:ceil(length(bitcode)/segsize),
        di = (si-1)*segsize+[1:segsize];
        di = di(ismember(di,1:length(bitcode)));
        plot(h,di,bitcode(di));
        hold on;
        xi = xInd(ismember(xInd,di));
        yi = ySignal(ismember(xInd,di));
        ci = stateList(ismember(xInd,di))+1;
        ci(ismember(xi,onsetList)) = 3;
        scatter(h,xi,yi,100,ci,'filled');
        title('Press Enter');
        hold off;
        i = input('');
    end
elseif visualize == 2,
    figure;
    h = axes;
    ci = [1 0 0];
    for ti=1:length(tempTrialsNew),
        timestamprange = [tempTrialsNew(ti).Timing.NS4.TimeStamp(1)-.7*tempSpecs.Recording.NS4SampleRes ...
                          tempTrialsNew(ti).Timing.NS4.TimeStamp(end)+.7*tempSpecs.Recording.NS4SampleRes];
        tinds = ismember(timeStamp,timestamprange(1):timestamprange(2));

        tchinds = zeros(1,length(tempTrialsNew(ti).Timing.NS4.TimeStamp));
        for chi=1:length(tchinds),
            [~, pos] = min(abs(timeStamp-tempTrialsNew(ti).Timing.NS4.TimeStamp(chi)));
            tchinds(chi) = pos;
        end

        plot(h,find(tinds),bitcode(tinds));
        hold on;
        xi = timeStamp(tchinds);
        yi = bitcode(tchinds);
        scatter(h,xi,yi,100,ci,'filled');
        title('Press Enter');
        hold off;
        i = input('');
    end
end