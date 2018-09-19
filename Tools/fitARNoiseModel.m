function [ arModel ] = fitARNoiseModel( noise, fitEpochs, maxLags )
    %Fits the multi-variate autoregressive noise model using least squares.
    %This is called in fitPiecewiseModel(). This function sweeps over
    %different numbers of lags and chooses a number of time lags to use
    %based on when R2 stops increasing. Fitting is complicated here by the
    %fact that fitEpochs can be discontinuous. 
    rIdx = [];
    for r=1:size(fitEpochs,1)
        rIdx = [rIdx, (fitEpochs(r,1)+maxLags):fitEpochs(r,2)];
    end
    
    nDim = size(noise,2);
    designMat = zeros(length(rIdx), maxLags * nDim);
    for r=1:length(rIdx)
        loopIdx = (rIdx(r)-1):-1:(rIdx(r)-maxLags);
        tmp = noise(loopIdx,:)';
        designMat(r,:) = tmp(:);
    end
    
    nFolds = 10;
    R2_vals = zeros(maxLags,nFolds,nDim);
    coef = cell(maxLags,nDim);
    nIdxPerFold = floor(size(designMat,1)/nFolds);
    allIdx = 1:(nIdxPerFold*nFolds);
    
    for noiseDim = 1:nDim
        for nLags=1:maxLags
            testIdx = 1:nIdxPerFold;
            for f=1:nFolds
                trainIdx = setdiff(allIdx, testIdx);
                coefInner = designMat(trainIdx, 1:(nLags*nDim)) \ noise(rIdx(trainIdx), noiseDim); 
                predNoise = designMat(testIdx, 1:(nLags*nDim)) * coefInner;
                
                truth = noise(rIdx(testIdx), noiseDim);
                SSERR = sum((predNoise - truth).^2);
                SSTOT = sum((truth-repmat(mean(truth),size(truth,1),1)).^2);
                R2_vals(nLags,f,noiseDim) = 1 - SSERR./SSTOT;

                testIdx = testIdx + nIdxPerFold;
            end
            coef{nLags, noiseDim} = designMat(:, 1:(nLags*nDim)) \ noise(rIdx, noiseDim); 
        end
    end
    
    %select the number of lags to use based on when R2 stops increasing
    lagsPerDim = zeros(nDim,1);
    meanR2 = squeeze(mean(R2_vals,2));
    for d=1:nDim
        for n=1:(maxLags-1)
            if mean(mean(R2_vals(n,:,:),3)) > 0.9*mean(mean(R2_vals(n+1,:,:),3))
                break;
            end
        end
        lagsPerDim(d) = n;
    end
    
    for d=1:nDim
        if mean(R2_vals(n,:,d),2)<0
            lagsPerDim(d) = 0;
        end
    end
    
    %Format the final model.
    %Remove outliers before getting the final covariance.
    tmpNoise = noise(rIdx,:);
    for d=1:nDim
        tmpNoise(abs(tmpNoise(:,d))>6*std(tmpNoise(:,d)),d)=0;
    end
    
    nLagsToUse = max(lagsPerDim);
    arModel.nLags = nLagsToUse;
    residualNoise = [];
    
    if arModel.nLags == 0
        arModel.coef = [];
    else
        arModel.coef = cell(nDim,1);
        for d=1:nDim
            arModel.coef{d} = coef{nLagsToUse, d};
            if lagsPerDim(d)>0
                residue = tmpNoise(:,d) - designMat(:, 1:(nLagsToUse*nDim)) * coef{nLagsToUse, d};
                residualNoise = [residualNoise, residue];
            else
                residualNoise = [residualNoise, tmpNoise(:,d)];
            end
        end
        arModel.covEps = cov(residualNoise);
    end
    arModel.covNoise = cov(tmpNoise);
    arModel.maxLags = maxLags;
    arModel.meanR2 = meanR2;
end

