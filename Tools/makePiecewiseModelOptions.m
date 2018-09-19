function [ opts ] = makePiecewiseModelOptions( )
    %This function makes a default arguments struct for the
    %fitPiecewiseModel() function that fits the PLM model to a given
    %dataset.
    
    %cursor state, target position, and decoded control vectors
    opts.pos = randn(10000,2);
    opts.vel = randn(10000,2);
    opts.targPos = randn(10000,2);
    opts.decoded_u = randn(10000,2);

    %model fitting options
    opts.modelOpts.noVel = false; %whether to fit fVel (if true, doesn't fit fVel)
    opts.modelOpts.nKnots = 12; %number of breakpoints in the piecewise linear functions
    opts.modelOpts.noNegativeFTarg = false; %if true, constrains fTarg to be positive

    %alpha and beta values of the decoder that was being used in this
    %dataset (used for to set up the user's forward model)
    opts.filtAlpha = 0.96; 
    opts.filtBeta = 1;

    opts.reachEpochsToFit = [1 100; 101 200; 201 300]; %will fit only on data within these epochs
    opts.feedbackDelaySteps = 10; %number of feedback delay steps for the user's forward model
    opts.timeStep = 0.02; %length of the time step, in seconds
    opts.fitNoiseModel = true; %will return an autoregressive noise model if true
    opts.fitSDN = true; %will fit a signal-dependent noise function if true
end
