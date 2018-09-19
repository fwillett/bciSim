//This function is for input checking and initializing the simulator.
//The actual simulation is done by simulator.c

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "simulator.h"
#include "pwl_interp_1d.h"

struct simulator sim;
int initialized = 0;

//Various input checking and utility functions
void checkFields(const mxArray *m, char fields [][20] , int numFields, char *structName);
void checkVectorLen(mxArray *vector, int len, char *errMsg);
void checkMatRows(mxArray *mat, int rows, char *errMsg);
void checkMatrixSizeEquality(mxArray *m1, mxArray *m2, char *errMsg);
void copyPwlFunction(double *x, double *y, int *nKnots, mxArray *x_src, mxArray *y_src);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    //These define the expected fields for the input struct. This function will throw an error if the input does not contain these fields.
    char *funcString;
    char fieldsOpts[][20] = {"trial","plant","forwardModel","noise","control","loopTime"}; 
    char fieldsTrial[][20] = {"dwellTime","maxTrialTime","continuousHoldRule"};
    char fieldsPlant[][20] = {"alpha","beta","nonlinType","n1","n2","fStaticX","fStaticY","nDim"};
    char fieldsForwardModel[][20] = {"delaySteps","forwardSteps"};
    char fieldsNoise[][20] = {"sdnX","sdnY"};
    char fieldsControl[][20] = {"fTargX","fTargY","fVelX","fVelY","rtSteps","targetDeadzone"};
    char fieldsRunOpts[][20] = {"noiseMatrix","noiseIdx","targetPos","initC","initX","targRad"}; 
    
    mxArray *trial;
	mxArray *plant;
	mxArray *forwardModel;
	mxArray *noise;
    mxArray *control;
    
    mxArray *initX;
    mxArray *initC;
    mxArray *noiseMatrix;
    mxArray *targetPos;
    const mxArray *opts;
    
    int nInitRows;
    int maxLoops;
    int x;
    int y;
    
    /* Check for proper number of arguments. */
    if ( nrhs != 2 ) {
        mexErrMsgTxt("This function takes two input arguments: a struct of data and a string specifying the sub-function to call.");
    } 
    
    funcString = mxArrayToString(prhs[1]);
    opts = prhs[0];
    
    //The function will 'initialize' or 'run' based on the second input.
    //'initialize' sets struct fields to prepare to run the simulator. 
    //'run' simulates a single movement. 
    if (funcString==NULL)
    {
        mexErrMsgTxt("Second input must be a string.");
    }
    else if (strcmp(funcString,"init")==0)
    {
        //Populates the simulator struct with all specifications.
        //Does a lot of input checking.
        if( nlhs != 0)
            mexErrMsgTxt("When calling 'init', must have zero outputs.");
    
        checkFields(opts,fieldsOpts,6,"opts");
        
        sim.loopTime = mxGetScalar(mxGetField(opts,0,"loopTime"));
        trial = mxGetField(opts,0,"trial");
        plant = mxGetField(opts,0,"plant");
        forwardModel = mxGetField(opts,0,"forwardModel");
        noise = mxGetField(opts,0,"noise");
        control = mxGetField(opts,0,"control");
        
        checkFields(trial,fieldsTrial,3,"trial");
        checkFields(plant,fieldsPlant,8,"plant");
        checkFields(forwardModel,fieldsForwardModel,2,"forwardModel");
        checkFields(noise,fieldsNoise,2,"noise");
        checkFields(control,fieldsControl,6,"control");
        	        
        sim.plant.nDim = (int)mxGetScalar(mxGetField(plant,0,"nDim"));
        sim.plant.alpha = mxGetScalar(mxGetField(plant,0,"alpha"));
        sim.plant.beta = mxGetScalar(mxGetField(plant,0,"beta"));
        sim.plant.n1 = mxGetScalar(mxGetField(plant,0,"n1"));
        sim.plant.n2 = mxGetScalar(mxGetField(plant,0,"n2"));
        sim.plant.nonlinType = (int)mxGetScalar(mxGetField(plant,0,"nonlinType"));
                
        sim.trial.dwellTime = mxGetScalar(mxGetField(trial,0,"dwellTime"));
        sim.trial.maxTrialTime = mxGetScalar(mxGetField(trial,0,"maxTrialTime"));
        sim.trial.continuousHoldRule = (int)mxGetScalar(mxGetField(trial,0,"continuousHoldRule"));
        
        sim.forwardModel.delaySteps = (int)mxGetScalar(mxGetField(forwardModel,0,"delaySteps"));
        sim.forwardModel.forwardSteps = (int)mxGetScalar(mxGetField(forwardModel,0,"forwardSteps"));
        
        sim.control.targetDeadzone = mxGetScalar(mxGetField(control,0,"targetDeadzone"));
        sim.control.rtSteps = (int)mxGetScalar(mxGetField(control,0,"rtSteps"));
        
        checkMatrixSizeEquality(mxGetField(noise,0,"sdnX"), mxGetField(noise,0,"sdnY"), "Dimensions of opts.noise.sdnX and opts.noise.sdnY should be equal");
        checkMatrixSizeEquality(mxGetField(control,0,"fTargX"), mxGetField(control,0,"fTargY"), "Dimensions of opts.control.fTargX and opts.control.fTargY should be equal");
        checkMatrixSizeEquality(mxGetField(control,0,"fVelX"), mxGetField(control,0,"fVelY"), "Dimensions of opts.control.fVelX and opts.control.fVelY should be equal");
        checkMatrixSizeEquality(mxGetField(plant,0,"fStaticX"), mxGetField(plant,0,"fStaticY"), "Dimensions of opts.plant.fStaticX and opts.plant.fStaticY should be equal");
                
        copyPwlFunction(sim.noise.sdnX, sim.noise.sdnY, &sim.noise.nsdn, mxGetField(noise,0,"sdnX"), mxGetField(noise,0,"sdnY"));
        copyPwlFunction(sim.control.fTargX, sim.control.fTargY, &sim.control.nfTarg, mxGetField(control,0,"fTargX"), mxGetField(control,0,"fTargY"));
        copyPwlFunction(sim.control.fVelX, sim.control.fVelY, &sim.control.nfVel, mxGetField(control,0,"fVelX"), mxGetField(control,0,"fVelY"));
        copyPwlFunction(sim.plant.fStaticX, sim.plant.fStaticY, &sim.plant.nfStatic, mxGetField(plant,0,"fStaticX"), mxGetField(plant,0,"fStaticY"));

        //keep track of whether we have successfully made it all the way through an initialization, in which case we can assume data is safe
        initialized = 1;
    }
    else if (strcmp(funcString,"run")==0)
    {
        //update target position and initial plant/control history values and then simulate a movement
        if( nlhs != 5)
            mexErrMsgTxt("When calling 'run', must have five outputs.");
        
        if(!initialized)
            mexErrMsgTxt("Initialize the model first by calling 'init'.");
        
        checkFields(opts,fieldsRunOpts,6,"opts");
        
        targetPos = mxGetField(opts,0,"targetPos");
        noiseMatrix = mxGetField(opts,0,"noiseMatrix");
        initX = mxGetField(opts,0,"initX");
        initC = mxGetField(opts,0,"initC");
        
        sim.trial.targRad = mxGetScalar(mxGetField(opts,0,"targRad"));
        
        checkVectorLen(targetPos, sim.plant.nDim, "Dimensions of opts.targetPos sohuld match opts.plant.nDim");
        memcpy(sim.trial.targetPos, mxGetPr(targetPos), sim.plant.nDim*sizeof(double));
        
        checkMatRows(noiseMatrix, sim.plant.nDim, "Number of rows in opts.noiseMatrix should be equal to opts.plant.nDim.");
        sim.noise.noiseMatrix = mxGetPr(noiseMatrix);
        sim.noise.noiseIdx = ((int)mxGetScalar(mxGetField(opts,0,"noiseIdx")))-1;
        sim.noise.nColsForNoiseMatrix = mxGetN(noiseMatrix);
        if(sim.noise.noiseIdx >= sim.noise.nColsForNoiseMatrix)
            mexErrMsgTxt("opts.noiseIdx is greater than the number of columns of opts.noiseMatrix.");
                
        checkMatRows(initC, sim.plant.nDim, "Number of rows in opts.initC should be equal to opts.plant.nDim.");
        checkMatRows(initX, 2*sim.plant.nDim, "Number of rows in opts.initX should be equal to 2*opts.plant.nDim.");
        
        if(mxGetN(initC)!=mxGetN(initX))
            mexErrMsgTxt("opts.initX and opts.initC should have the same number of columns.");
        
        nInitRows = mxGetN(initC);
        if(nInitRows==0){
            mexErrMsgTxt("opts.initX and opts.initC must have at least one column in order to initialize the cursor state.");
        }
        else if(nInitRows<sim.forwardModel.delaySteps){
            mexErrMsgTxt("opts.initX and opts.initC must have at least as many columns as opts.forwardModel.delaySteps.");
        }
                
        //allocate memory for state and control vector matrices
        sim.maxLoops = nInitRows + ((int)ceil(sim.trial.maxTrialTime / sim.loopTime));
                
        sim.xMatrix = mxCalloc(2 * sim.plant.nDim * sim.maxLoops, sizeof(double));
        sim.xHatMatrix = mxCalloc(2 * sim.plant.nDim * sim.maxLoops, sizeof(double));
        
        sim.uMatrix = mxCalloc(sim.plant.nDim * sim.maxLoops, sizeof(double));
        sim.cMatrix = mxCalloc(sim.plant.nDim * sim.maxLoops, sizeof(double));
              
        //initialize history
        memcpy(sim.xMatrix, mxGetPr(initX), (2 * sim.plant.nDim * nInitRows)*sizeof(double));
        memcpy(sim.cMatrix, mxGetPr(initC), (sim.plant.nDim * nInitRows)*sizeof(double));
        sim.loopIdx = nInitRows;
        
        //simulate
        simulate(&sim);
        
        //return results matrices
        plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        plhs[3] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        plhs[4] = mxCreateDoubleScalar(sim.loopIdx);
        
        mxSetPr(plhs[0], sim.xMatrix);
        mxSetM(plhs[0], 2*sim.plant.nDim);
        mxSetN(plhs[0], sim.maxLoops);
        
        mxSetPr(plhs[1], sim.xHatMatrix);
        mxSetM(plhs[1], 2*sim.plant.nDim);
        mxSetN(plhs[1], sim.maxLoops);
        
        mxSetPr(plhs[2], sim.uMatrix);
        mxSetM(plhs[2], sim.plant.nDim);
        mxSetN(plhs[2], sim.maxLoops);
        
        mxSetPr(plhs[3], sim.cMatrix);
        mxSetM(plhs[3], sim.plant.nDim);
        mxSetN(plhs[3], sim.maxLoops);
    }
    else
    {
        mexErrMsgTxt("The second input must equal \"init\" or \"run\".");
    }

    mxFree(funcString);
    return;
}

//--sub-functions for input checking--
void checkFields(const mxArray *m, char fields [][20], int numFields, char *structName)
{
	int f;
    char errMsg[200];
	for(f=0; f<numFields; f++)
	{
		if(mxGetField(m,0,fields[f])==NULL)
		{
            sprintf(errMsg,"ERROR: Field \"%s\" not found in %s \n",fields[f],structName);
            mexErrMsgTxt(errMsg);
		}
	}
}

void checkVectorLen(mxArray *vector, int len, char *errMsg)
{
	if(!((mxGetM(vector)==len) && (mxGetN(vector)==1)) &&
		!((mxGetM(vector)==1) && (mxGetN(vector)==len)))
	{
		mexErrMsgTxt(errMsg);
	}
}

void checkMatRows(mxArray *mat, int rows, char *errMsg)
{
	if(mxGetM(mat)!=rows)
	{
		mexErrMsgTxt(errMsg);
	}
}

void checkMatrixSizeEquality(mxArray *m1, mxArray *m2, char *errMsg)
{
    if((mxGetN(m1)!=mxGetN(m2)) || (mxGetM(m1)!=mxGetM(m2)))
    {
        mexErrMsgTxt(errMsg);
    }
}

void copyPwlFunction(double *x, double *y, int *nKnots, mxArray *x_src, mxArray *y_src)
{
    if( mxGetM(x_src) > mxGetN(x_src) )
        *nKnots = (int)mxGetM(x_src);
    else
        *nKnots = (int)mxGetN(x_src);

    memcpy(x, mxGetPr(x_src), (*nKnots)*sizeof(double));
    memcpy(y, mxGetPr(y_src), (*nKnots)*sizeof(double));
}

