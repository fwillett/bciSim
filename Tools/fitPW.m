function [ model, predVals ] = fitPW( cVec, posHat, velHat, targPos, fitOpts )
    %Fits a piecewise model of the user's control policy. 
    targDist = matVecMag(targPos - posHat,2);
    speed = matVecMag(velHat,2);
    nDim = size(cVec, 2);

    %find the knots (breakpoints) based on percentiles of the data
    maxDist = prctile(targDist, 99);
    maxSpeed = prctile(speed, 98);
    targDist(targDist>maxDist) = maxDist;
    speed(speed>maxSpeed) = maxSpeed;
    
    distEdges = prctile(targDist,linspace(0,100,fitOpts.nKnots));
    distEdges = unique(distEdges);
    nKnots = length(distEdges);
    speedEdges = prctile(speed,linspace(0,100,nKnots));
    
    %make the design matrices for piecewise linear fitting
    [ distExp ] = cpwlDesignMatrix( distEdges, targDist );
    [ velExp ] = cpwlDesignMatrix( speedEdges, speed );
    
    %expand weights features into vectors
    toTargVec = bsxfun(@times, targPos - posHat, 1./matVecMag(targPos - posHat,2));
    toTargVec(isnan(toTargVec))=0;
    velUnitVec = bsxfun(@times, velHat, 1./matVecMag(velHat,2));
    velUnitVec(isnan(velUnitVec))=0;
    
    distVec = [];
    totalVelVec = [];
    for t=1:nKnots
        distVec = [distVec, bsxfun(@times, toTargVec, distExp(:,t))];
        totalVelVec = [totalVelVec, bsxfun(@times, velUnitVec, velExp(:,t))];
    end
    if fitOpts.noVel
        goodIdxOneDim = ~all(distVec==0,2);
    else
        goodIdxOneDim = ~all(distVec==0,2) & ~all(totalVelVec==0,2);
    end
    goodIdxStack = repmat(goodIdxOneDim, nDim, 1);
    
    %stack dimensions
    distVecStack = [];
    totalVelVecStack = [];
    for n=1:nDim
        distVecStack = [distVecStack; distVec(:,n:nDim:end)];
        totalVelVecStack = [totalVelVecStack; totalVelVec(:,n:nDim:end)];
    end
    
    %special cases of the model
    if fitOpts.noVel
        predictors = distVecStack;
        LB = -inf(1,nKnots + nDim);
        UB = inf(1,nKnots + nDim);
    else
        predictors = [distVecStack, totalVelVecStack];
        LB = -inf(1, nKnots*2 + nDim);
        UB = [inf(1,nKnots), zeros(1,nKnots), inf(1,nDim)];
    end
    
    %option to clamp fTarg at 0
    if fitOpts.noNegativeFTarg
        LB(1:nKnots) = 0;
    end
    
    %bias offsets
    biasCols = zeros(size(cVec,1)*nDim, nDim);
    biasCols(1:size(cVec,1), 1) = 1;
    biasCols((1:size(cVec,1))+size(cVec,1), 2) = 1;
    predictors = [predictors, biasCols];
    
    %stack columns
    response = cVec(:); 
    
    %remove idx
    predictors = predictors(goodIdxStack,:);
    response = response(goodIdxStack,:); 
    
    %quad prog, constrain fVel to be negative
    A = predictors'*predictors;
    q = zeros(size(predictors,2),1);
    for n=1:size(predictors,2)
        q(n) = -sum(predictors(:,n).*response);
    end
    quadOpts = optimoptions('quadprog', 'display', 'off');
    coefFinal = quadprog(A,q,[],[],[],[],LB,UB,[],quadOpts);
        
    if fitOpts.noVel
        model.fTargX = distEdges;
        model.fTargY = coefFinal(1:nKnots)';
        model.fVelX = [];
        model.fVelY = [];
        model.bias = coefFinal((end-nDim+1):end);
    else
        model.fTargX = distEdges;
        model.fTargY = coefFinal(1:nKnots)';
        model.fVelX = speedEdges;
        model.fVelY = coefFinal((nKnots+1):(2*nKnots))';
        model.bias = coefFinal((end-nDim+1):end);
    end
    
    predVals = applyPiecewiseModel( model, posHat, velHat, targPos );
end

function [ mag ] = matVecMag( mat, dim )
    if length(size(mat))>2
        error('Cannot use on matrices bigger than 2 dimensions');
    end
    mag = sqrt(sum(mat.^2,dim));
end
    