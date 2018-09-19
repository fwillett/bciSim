function [ modelOut ] = fitPiecewiseModel( opts )
    
    %Get a vector of sample indices on which we actually fit the control
    %policy model. Typically these are indices during which the trials
    %occurred, excluding the user's reaction time.
    rIdx = expandEpochIdx(opts.reachEpochsToFit);
    
    %The offset convention variable is used to describe whether the decoded
    %velocity is applied to update the position at the same time step on which it is decoded (offsetConvention=0), or
    %whether it applies to the next time step (offsetConvention=1). 
    opts.offsetConvention = 1;
    
    %Concatenate velocity and position together as a state vector.
    posIdx = 1:size(opts.pos,2);
    velIdx = (size(opts.pos,2)+1):(2*size(opts.pos,2));
    effectorStates = [opts.pos, opts.vel];

    %Start by using purely delayed states as an initial guess for the
    %internal model estimates.
    delayedStates = [zeros(opts.feedbackDelaySteps+opts.offsetConvention,size(effectorStates,2)); effectorStates(1:(end-opts.feedbackDelaySteps),:)];
    internalStates = delayedStates;
    
    %Here, we iteratively fit a control policy model and update the internal model
    %estimates of cursor position and velocity.
    nIterations = 4;

    for n=1:nIterations
        %Fit the piecewise control policy model.
        [model, predVals] = fitPW( opts.decoded_u(rIdx,:), internalStates(rIdx,posIdx), internalStates(rIdx,velIdx), opts.targPos(rIdx,:), opts.modelOpts);
           
        %Update our best guess of the user's internal model estimates.
        cVecPredicted = zeros(size(opts.decoded_u));
        cVecPredicted(rIdx,:) = predVals;
        internalStates = getInternalModelState(effectorStates, opts.feedbackDelaySteps, opts.filtAlpha, opts.filtBeta, opts.timeStep, ...
            cVecPredicted, opts.offsetConvention );
    end
    
    %Fit an autoregressive noise model if directed.
    maxLags = ceil(0.400 / opts.timeStep);
    if opts.fitNoiseModel
        arModel = fitARNoiseModel( opts.decoded_u - cVecPredicted, opts.reachEpochsToFit, maxLags );
        modelOut.noiseModel = arModel;
    end
        
    %Fit a signal-dependent noise model if directed.
    if opts.fitSDN && opts.fitNoiseModel
        cVecMag = sqrt(sum(cVecPredicted.^2,2));
        maxMag = max(cVecMag);
        
        binEdges = linspace(0,maxMag,21);
        binStd = zeros(length(binEdges)-1,1);
        binCenters = binEdges(1:(end-1)) + (binEdges(2:end)-binEdges(1:(end-1)))/2;
        noiseTimeSeries = opts.decoded_u - cVecPredicted;
        
        for b=1:(length(binEdges)-1)
            binIdx = intersect(rIdx, find(cVecMag>binEdges(b) & cVecMag<=binEdges(b+1)));
            binStd(b) = std(noiseTimeSeries(binIdx));
        end
        binStd(isnan(binStd)) = nanmean(binStd);
        stdSmoothed = filtfilt(ones(5,1)/5, 1, binStd);

        ts = generateNoiseFromModel( 10000, modelOut.noiseModel );
        modelOut.sdnX = binCenters;
        modelOut.sdnY = stdSmoothed' / mean(std(ts));
    else
        modelOut.sdnX = 1;
        modelOut.sdnY = 1;
    end
    
    %Group output variables into modelOut struct.
    modelOut.controlModel = model;
    modelOut.internalModelEstimates = internalStates;
    modelOut.modeledControlVector = cVecPredicted;
    
    %Make a struct for simulating movements with this model.
    simOpts = makeBciSimOptions( );
    simOpts.loopTime = opts.timeStep;
    simOpts.plant.alpha = opts.filtAlpha;
    simOpts.plant.beta = opts.filtBeta;
    simOpts.forwardModel.delaySteps = opts.feedbackDelaySteps;
    simOpts.forwardModel.forwardSteps = opts.feedbackDelaySteps;
    simOpts.control.fTargX = modelOut.controlModel.fTargX;
    simOpts.control.fTargY = modelOut.controlModel.fTargY;
    simOpts.control.fVelX = modelOut.controlModel.fVelX;
    simOpts.control.fVelY = modelOut.controlModel.fVelY;
    simOpts.control.rtSteps = opts.feedbackDelaySteps;
    if opts.fitNoiseModel
        simOpts.noiseMatrix = generateNoiseFromModel( 100000, modelOut.noiseModel );
    else 
        simOpts.noiseMatrix = randn(100000, size(opts.decoded_u,2));
    end
    simOpts.noise.sdnX = modelOut.sdnX;
    simOpts.noise.sdnY = modelOut.sdnY;
    modelOut.simOpts = simOpts;
end

%Utility function that takes a matrix of trial start and end indices and
%returns a vector of all indices inside those trials. 
function [ idx ] = expandEpochIdx( epochs )
    idx = [];
    for e=1:size(epochs,1)
        idx = [idx, epochs(e,1):epochs(e,2)];
    end
end
