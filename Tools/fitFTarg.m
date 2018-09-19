function [fTarg, weightVec] = fitFTarg(posErr, decVectors, maxDist, nKnots, noNegative)
    
    if nargin<5
        noNegative = false;
    end
    
    targDist = sqrt(sum(posErr.^2,2));
    goodIdx = targDist<=maxDist;
    nDim = size(posErr, 2);
    
    distNorm = targDist/maxDist;
    distEdges = prctile(distNorm(goodIdx),linspace(0,100,nKnots));
    distEdges = unique(distEdges);
    nKnots = length(distEdges);
    distExp = cpwlDesignMatrix( distEdges, distNorm );

    %expand weights features into vectors
    toTargVec = bsxfun(@times, posErr, 1./targDist);
    toTargVec(isnan(toTargVec))=0;
    
    distVec = [];
    for t=1:nKnots
        distVec = [distVec, bsxfun(@times, toTargVec, distExp(:,t))];
    end
    goodIdx = ~all(distVec==0,2) & goodIdx;
    goodIdxStack = repmat(goodIdx, nDim, 1);
        
    %stack dimensions
    distVecStack = [];
    for n=1:nDim
        distVecStack = [distVecStack; distVec(:,n:nDim:end)];
    end
    response = decVectors(:); %stack columns
    
    %quad prog, constrain fVel to be negative
    if noNegative
        predictors = distVecStack(goodIdxStack,:);
        response = response(goodIdxStack,:);
        
        A = predictors'*predictors;
        q = zeros(size(predictors,2),1);
        for n=1:size(predictors,2)
            q(n) = -sum(predictors(:,n).*response);
        end
        quadOpts = optimoptions('quadprog', 'display', 'off');
          
        LB = zeros(1,nKnots);
        UB = inf(1,nKnots);  
        coefFinal = quadprog(A,q,[],[],[],[],LB,UB,[],quadOpts);
    else
        coefFinal = distVecStack(goodIdxStack,:) \ response(goodIdxStack,:);
    end
    
    distEdges = distEdges*maxDist;
    fTarg = [distEdges', coefFinal];
    weightVec = bsxfun(@times, toTargVec, interp1(fTarg(:,1), fTarg(:,2), targDist));
end