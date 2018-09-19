# Overview

bciSim is a set of functions for simulating, modeling and parameterizing linear velocity decoders with exponential smoothing (including steady-state velocity Kalman filters). This folder contains the functions, example scripts, and an example dataset.

# Example scripts

- To run the example scripts in the "Examples" folder, the "Sample Dataset" and "Tools" folders must be on the MATLAB path. 

- Examples\exampleFittingAndSimulations.m gives an example of how to fit the control policy model first introduced in Willett et. al 2016 ("Feedback control policies employed by people using intracortical brain-computer interfaces") to one block of closed-loop data. It also shows how to use that model to predict performance under different gain and smoothing conditions. It uses the example dataset T8.2015.03.24 (a "high gain" session where both high and low decoder gains were tested).

- Examples\exampleSimulations.m shows how to use the simulator to rapidly simulate cursor movements. The simulator can be used to predict which gain and smoothing parameters will be optimal for online performance with a specific user and on a specific task. The simulator is written in C and has a mex interface. It will have to be compiled for your system (which can be accomplished with compileSimBci.m).

- Examples\exampleReparameterization.m shows how to reparameterize a setady-state velocity Kalman filter into the (alpha, beta, D) parameterization used in the paper “Feedback control policies employed by people using intracortical brain-computer interfaces” (and others). This parameterization allows easy reporting and understanding of the Kalman filter's gain and smoothing properties.

# Key Functions

- Tools\reparamKalman.m converts a steady-state velocity Kalman filter to the (alpha, beta, D) parameterization.

- Tools\simBci.mex is the mex interface to the simulator. It is called to simulate a single trajectory. It requires the simulation options to be specified with an options struct that can be created with makeBciSimOptions.m

- Tools\simBatch.m can be used to simulate a batch of cursor trajectories. 

- Tools\fitPiecewiseModel.m can be used to fit a control policy model (and a corresponding noise model) to closed-loop cursor control data. It requires an options struct that can be created with makePiecewiseModelOptions.m

# Sample Dataset T8.2015.03.24

- Contains the cursor kinematics and decoder output (but not neural data) from one of the datasets analyzed in Willett et. al 2016 ("Feedback control policies employed by people using intracortical brain-computer interfaces"). It is used in the script exampleFittingAndSimulations.m.

- This sample dataset contains 15 different blocks. Each block is collected under 1 of 4 possible decoder gains. The “alpha”, “beta” and “dwellTimes” variables give the alpha (smoothing), beta (gain), and dwell time (in seconds) used on any given block. The “secondsSinceSystemBoot”, “conditionNumber”, “cursorPos”, “cursorVel”, “targetPos” and “decodedControlVector” variables contain time series data where each sample corresponds to one row. “cursorVel” is units of cm/sec. The “trialEpochs” variable contains the trial start and end times (in # of samples) for each of the 1869 trials.
