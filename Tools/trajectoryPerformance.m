function [dialTime, transTime, totalTime, pathEff] = ...
    trajectoryPerformance(cursorPos, targetPos, acquisitionRad, trialEpochs, loopTime)
    
    %Computes four summary statistics for each trial to summarize cursor
    %control performance (dial-in time, translation time, total trial time,
    %and path efficiency). 
    
    nTrials = size(trialEpochs,1);
    dialTime = zeros(nTrials,1);
    transTime = zeros(nTrials,1);
    pathEff = zeros(nTrials,1);
    
    targetDist = sqrt(sum((targetPos-cursorPos).^2,2));
    inTarget = targetDist <= acquisitionRad;
    totalTime = (trialEpochs(:,2)-trialEpochs(:,1)+1)*loopTime;
    
    for t=1:nTrials
        loopIdx = trialEpochs(t,1):trialEpochs(t,2);
        touchIdx = find(inTarget(loopIdx),1,'first');
        if isempty(touchIdx)
            dialTime(t) = NaN;
            transTime(t) = totalTime(t);
        else
            dialTime(t) = (length(loopIdx)-touchIdx+1)*loopTime;
            transTime(t) = (touchIdx-1)*loopTime;
        end
        pathEff(t) = targetDist(loopIdx(1)) / sum(sqrt(sum(diff(cursorPos(loopIdx,:)).^2,2)));
    end
end