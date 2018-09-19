function [ internalStates ] = getInternalModelState( effectorStates, feedbackSteps, alpha, beta, timeStep, controlVectors, offsetConvention )
    %This function computes the user's internal model estimate of cursor
    %position and velocity. At each time step, it starts from a delayed
    %cursor position and velocity ("knownState") and then steps forward using a forward
    %model (implemented in cursorForwardFcn) to compute the "internalState". 
    
    %The offset convention variable is used to describe whether the decoded
    %velocity is applied to update the position at the same time step on which it is decoded (offsetConvention=0), or
    %whether it applies to the next time step (offsetConvention=1). 
    
    internalStates = zeros(size(effectorStates));
    for s = 1:size(internalStates,1)
        knownStateIdx = s - feedbackSteps - offsetConvention;
        if knownStateIdx < 1
            knownStateIdx = 1;
        end
        
        knownState = effectorStates(knownStateIdx,:);
        if feedbackSteps > 0 && knownStateIdx > 1
            internalStates(s,:) = cursorForwardFcn(knownState, controlVectors((knownStateIdx+1):(s-1),:), alpha, beta, timeStep);
        else
            internalStates(s,:) = knownState;
        end
    end
end

function [ cursorState ] = cursorForwardFcn( startState, controlVectors, alpha, beta, timeStep )
    %The first half of the cursorState variables  are position, the latter half are
    %velocities.
    cursorState = startState;
    for c=1:size(controlVectors,1)
        cursorState(((end/2)+1):end) = cursorState(((end/2)+1):end)*alpha + (1-alpha)*beta*controlVectors(c,:);
        cursorState(1:(end/2)) = cursorState(1:(end/2)) + cursorState(((end/2)+1):end)*timeStep;
    end
end

