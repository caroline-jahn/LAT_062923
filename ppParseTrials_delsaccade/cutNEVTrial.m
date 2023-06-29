function [tempTrials, tempSpecs] = cutNEVTrial(NEV, trials, specs, Semaphore)

%%
% VERY Buggy encodes that occur on next timestamp !!!
% Problematic encodes sequence
% Parses trial in FIX_ON [first fixation] CUE_ON DELAY_ON SAMPLE_ON format
% Add comments %%% Do sanity checks !!!

%%
nevF = NEV.MetaTags.SampleRes; %sampling frequency of nev signal
encode = NEV.Data.SerialDigitalIO.UnparsedData'; %encode values of events
timeStamp = NEV.Data.SerialDigitalIO.TimeStamp; %timestamp serial numbers of events
timeStampSec = NEV.Data.SerialDigitalIO.TimeStampSec; %time in secs of events
%timeStamp = timeStampSec*nevF

scEncodesList = 40:49; %[specs.Encodes.<TRIAL_STOP_CONDITION_ENCODE>]
%startEncodesList = specs.Encodes.COND_SEMAPHORE;
startEncodesList = Semaphore;
trialNumOffset = specs.Encodes.TRIAL_NUM_OFFSET;
numsteps = 3;

%%
parse = zeros(length(trials),numsteps); %points of trials (beginning, trialno, sc)
prevTrialEndIndex = 1;
%go through entire encode array and assign encodes to trials in parse array
while prevTrialEndIndex < length(encode) - length(startEncodesList)
    
    %Look for the exact semaphore
    Test_semaphore=zeros(length(startEncodesList),1);
    for i=1:length(startEncodesList)
         if encode(prevTrialEndIndex+i-1)==startEncodesList(i) 
             Test_semaphore(i)=1;
         end
    end
      
    if sum(Test_semaphore)==length(startEncodesList)
        curTrialEndIndex = prevTrialEndIndex;
        
        while curTrialEndIndex < length(encode) && ~ismember(encode(curTrialEndIndex),scEncodesList)
            curTrialEndIndex = curTrialEndIndex+1;
        end
        
        if curTrialEndIndex <= length(encode)
            indsRange = prevTrialEndIndex:curTrialEndIndex;
        else
            break
        end
        
        encodeSegment = encode(indsRange);
        timeStampSegment = timeStamp(indsRange);
        startEncodesInds = ismember(encodeSegment,startEncodesList);
        trialNumsInds = encodeSegment > trialNumOffset;
        trialNums = encodeSegment(trialNumsInds)-trialNumOffset;
        if ismember(encode(curTrialEndIndex),scEncodesList) && ...                          %if segment ends with stop condition
                ~isempty(startEncodesInds) && ...                                                %if there are start semaphores
                numel(unique(trialNums)) == 1 && ...                                             %if single trial num exists  all(arrayfun(@(x) all(x > find(startEncodesInds)),find(trialNumsInds))) && ...   %if all trial nums follow all start semaphores
                ~any(diff(timeStampSegment) == 1)                                                %if there are no bug encodes (next timestamp)
            trialNum = unique(trialNums);
            parse(trialNum,:) = [indsRange(find(startEncodesInds,1,'first')) ...
                indsRange(find(trialNumsInds,1,'first')) ...
                curTrialEndIndex];
        end
        
        prevTrialEndIndex = curTrialEndIndex;
        
    elseif prevTrialEndIndex > length(encode)
        break
    else
        prevTrialEndIndex = prevTrialEndIndex + 1;
    end
    
    
    
end

%%
for ti=1:length(trials)
    
    tempTrials(ti).Timing.NEV = struct('StageSequence', NaN, 'TimeStamp', NaN, 'TimeStampSec', NaN);
    
    if all(parse(ti,:) > 0)
        timeStampInds = parse(ti,2):parse(ti,3);
        trialSegment.encode = encode(timeStampInds);
        trialSegment.timeStamp = timeStamp(timeStampInds);
        trialSegment.timeStampSec = timeStampSec(timeStampInds);
        
        %construct stage encode sequence for trial . %% THIS IS TASK
        %SPECIFIC see specs.Task.StageNames
        trialStageSegment.encode = zeros(1,trials(ti).Behavior.StageNo);
         for si=1:trials(ti).Behavior.StageNo
            if si == 1, trialStageSegment.encode(si) = specs.Encodes.FIX_ON;
            elseif si == 2, trialStageSegment.encode(si) = specs.Encodes.SAMPLE_ON;
            elseif si == 3, trialStageSegment.encode(si) = specs.Encodes.DELAY_ON;
            elseif si == 4, trialStageSegment.encode(si) = specs.Encodes.TARGET_ON;
            elseif si == 5, trialStageSegment.encode(si) = specs.Encodes.FEEDBACK_ON;    
            end
        end
        
        %init stage times
        trialStageSegment.timeStamp = zeros(1,trials(ti).Behavior.StageNo);
        trialStageSegment.timeStampSec = zeros(1,trials(ti).Behavior.StageNo);
        
        %find last fix stage for assigning fix timestamp
        si=1;
        lastFix = find(trialSegment.encode == specs.Encodes.FIX_ON,1,'first');
        trialStageSegment.timeStamp(si) = trialSegment.timeStamp(lastFix);
        trialStageSegment.timeStampSec(si) = trialSegment.timeStampSec(lastFix);
        
        si = 2;
        segmentIndex = lastFix;
        %while not finished with all stages
        while si <= trials(ti).Behavior.StageNo
            segmentIndex = segmentIndex+1;
            %if found next stage within segment
            if segmentIndex>length(trialSegment.encode) || si>length(trialStageSegment.encode)
                %                 si
                %                 segmentIndex
                si=si+1;
            else
                
                if trialSegment.encode(segmentIndex) == trialStageSegment.encode(si),
                    %assign and go to next stage
                    trialStageSegment.timeStamp(si) = trialSegment.timeStamp(segmentIndex);
                    trialStageSegment.timeStampSec(si) = trialSegment.timeStampSec(segmentIndex);
                    si = si+1;
                end
            end
        end
        
        tempTrials(ti).Timing.NEV.StageSequence = trialStageSegment.encode;
        tempTrials(ti).Timing.NEV.TimeStamp = trialStageSegment.timeStamp;
        tempTrials(ti).Timing.NEV.TimeStampSec = trialStageSegment.timeStampSec;
    end
    
end

tempSpecs.Recording.NEVSampleRes = nevF;