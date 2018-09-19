#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "mex.h"

#define MAX_PWL_KNOTS 100
#define MAX_DIM 100

//the MATLAB function makeBciSimOptions defines many of these variables

struct simTrial {
    double targetPos[MAX_DIM];    
    double dwellTime;
    double maxTrialTime;
    double targRad;
    int continuousHoldRule;
};

struct simPlant {
    double alpha;
    double beta;
    
    int nDim;
    
    int nonlinType;
    double n1;
    double n2;
    
    int nfStatic;
    double fStaticX[MAX_PWL_KNOTS];
    double fStaticY[MAX_PWL_KNOTS];
};

struct simForwardModel {
    int delaySteps;
    int forwardSteps;
};

struct simNoise {  
    double *noiseMatrix; 
    int noiseIdx;
    int nColsForNoiseMatrix;
    
    double sdnX[MAX_PWL_KNOTS];
    double sdnY[MAX_PWL_KNOTS];
    int nsdn; 
};

struct simController {
    double fTargX[MAX_PWL_KNOTS];
    double fTargY[MAX_PWL_KNOTS];
    int nfTarg;
    
    double fVelX[MAX_PWL_KNOTS];
    double fVelY[MAX_PWL_KNOTS];
    int nfVel;
    
    double targetDeadzone;
    int rtSteps;
};

struct simulator {
    int loopIdx;
    int maxLoops;
    double loopTime;
    
    //These are (nDim x maxLoops) or (2*nDim x maxLoops) column major matrices that
    //store the cursor state (xMatrix), decoded control vector (uMatrix), control vector (cMatrix),
    //and the internal model estiamte of the cursor state (xHatMatrix)
    double *xMatrix;
    double *uMatrix;
    double *cMatrix;
    double *xHatMatrix;
    
    struct simTrial trial;
    struct simForwardModel forwardModel;
    struct simPlant plant;
    struct simNoise noise;
    struct simController control;
};

void simulate(struct simulator *sim);

#endif