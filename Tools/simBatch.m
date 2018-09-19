function [ out ] = simBatch( opts, targPos, startPos )
    %out = simBatch( opts, targPos, startPos ) simulates a batch of movements
    %and returns the trajectories in a struct.
    %
    %opts is a struct of simulation parameters that can be created with
    %makeBciSimOptions(). See makeBciSimOptions for explanations of each
    %parameters and default values.
    %
    %targPos is an N x D matrix of target positions, where N = number of
    %trials to simulated and D = number of dimensions (i.e. 2 for a 2D
    %cursor control task).
    %
    %startPos is either an N x D matrix of starting positions or a 1 x D
    %row vector. If startPos is a single row, then each movement will begin
    %at the end point of the previous one. If startPos has as many rows as targPos, the
    %cursor will be reset to startPos for every trial.
    
    nTrials = size(targPos,1);
    resetCursor = size(startPos,1)==size(targPos,1);
    
    %information stored reach-wise
    out.movTime = zeros(nTrials,1);
    out.reachEpochs = zeros(nTrials, 2);
    
    %information stored loop-wise
    maxLoopsPerTrial = 1+ceil(opts.trial.maxTrialTime/opts.loopTime);
    maxLoops = maxLoopsPerTrial * nTrials;
    globalLoopIdx = 1;
    out.pos = zeros(maxLoops, opts.plant.nDim);
    out.vel = zeros(maxLoops, opts.plant.nDim);
    out.posHat = zeros(maxLoops, opts.plant.nDim);
    out.velHat = zeros(maxLoops, opts.plant.nDim);
    out.targPos = zeros(maxLoops, opts.plant.nDim);
    out.controlVec = zeros(maxLoops, opts.plant.nDim);
    out.decVec = zeros(maxLoops, opts.plant.nDim);
    
    %initialize variables that won't change from movement to movement
    simBci(opts, 'init');
    
    %prepare runOpts struct that will change from movement to movement
    runOpts.noiseMatrix = opts.noiseMatrix';
    runOpts.noiseIdx = 1;
    runOpts.initC = zeros(opts.plant.nDim, opts.forwardModel.delaySteps + 1);
    runOpts.initX = zeros(2*opts.plant.nDim, opts.forwardModel.delaySteps + 1);
    runOpts.initX(1:opts.plant.nDim,end) = startPos(1,:);
    runOpts.targRad = opts.trial.targRad;
    
    %simulate one reach at a time, and record movement data and basic
    %performance metrics
    for r=1:nTrials
        runOpts.targetPos = targPos(r,:);
        if resetCursor
            %initialize forward model history to zero if the cursor gets
            %reset
            runOpts.initC = zeros(opts.plant.nDim, opts.forwardModel.delaySteps + 1);
            runOpts.initX = zeros(2*opts.plant.nDim, opts.forwardModel.delaySteps + 1);
            
            %reset cursor position
            runOpts.initX(1:opts.plant.nDim,end) = startPos(r,:);
        end
        
        %simulate the movement
        [xMatrix, xHatMatrix, uMatrix, cMatrix, simLoopIdx] = simBci(runOpts, 'run');
        xMatrix = xMatrix';
        xHatMatrix = xHatMatrix';
        uMatrix = uMatrix';
        cMatrix = cMatrix';
        
        %number of loops the movement lasted
        nLoops = simLoopIdx - opts.forwardModel.delaySteps;
        
        %advance the noise index forward so the next movement has different
        %noise
        runOpts.noiseIdx = runOpts.noiseIdx + nLoops;
        if runOpts.noiseIdx > size(runOpts.noiseMatrix,2)
            runOpts.noiseIdx = 1;
        end

        %store information from this movement
        out.movTime(r) = nLoops * opts.loopTime;

        if ~resetCursor && r>1
            loopIdxToUse = (opts.forwardModel.delaySteps+2):simLoopIdx;
        else
            loopIdxToUse = (opts.forwardModel.delaySteps+1):simLoopIdx;
        end
        nLoops = length(loopIdxToUse);
        entryIdx = globalLoopIdx:(globalLoopIdx+nLoops-1);
        
        pos = xMatrix(loopIdxToUse,1:opts.plant.nDim);
        out.pos(entryIdx,:) = pos;
        out.vel(entryIdx,:) = xMatrix(loopIdxToUse,(opts.plant.nDim+1):(2*opts.plant.nDim));
        out.posHat(entryIdx,:) = xHatMatrix(loopIdxToUse,1:opts.plant.nDim);
        out.velHat(entryIdx,:) = xHatMatrix(loopIdxToUse,(opts.plant.nDim+1):(2*opts.plant.nDim));
        out.targPos(entryIdx,:) = repmat(runOpts.targetPos,length(loopIdxToUse),1);
        out.controlVec(entryIdx,:) = cMatrix(loopIdxToUse,:);
        out.decVec(entryIdx,:) = uMatrix(loopIdxToUse,:);
        out.reachEpochs(r,:) = [globalLoopIdx, globalLoopIdx + nLoops - 1];
        
        globalLoopIdx = globalLoopIdx + nLoops;
 
        %prepare to start the next movement with the final state of this movement
        if ~resetCursor
            lastIdx = globalLoopIdx - 1;
            tmpIdx = (lastIdx-opts.forwardModel.delaySteps):lastIdx;
            negativeIdx = tmpIdx<1;
            if sum(negativeIdx)>0
                runOpts.initC = [zeros(size(cMatrix,2),sum(negativeIdx)), out.controlVec(tmpIdx(~negativeIdx),:)'];
                runOpts.initX = [zeros(size(xMatrix,2),sum(negativeIdx)), [out.pos(tmpIdx(~negativeIdx),:), out.vel(tmpIdx(~negativeIdx),:)]'];               
            else
                runOpts.initC = out.controlVec(tmpIdx,:)';
                runOpts.initX = [out.pos(tmpIdx,:), out.vel(tmpIdx,:)]';
            end
        end
    end
    
    keepIdx = 1:(globalLoopIdx-1);
    out.pos = out.pos(keepIdx,:);
    out.vel = out.vel(keepIdx,:);
    out.posHat = out.posHat(keepIdx,:);
    out.velHat = out.velHat(keepIdx,:);
    out.targPos = out.targPos(keepIdx,:);
    out.controlVec = out.controlVec(keepIdx,:);
    out.decVec = out.decVec(keepIdx,:);
end

