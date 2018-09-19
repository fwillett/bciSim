function [ ts ] = generateNoiseFromModel( nPoints, arModel )
    %Applies the multi-variate autoregressive model to generate a simulated time series history
    %consistent with the model.
    
    %More recent lags are at the top of arCoeff.
    if arModel.nLags==0
        ts = mvnrnd(zeros(length(arModel.covNoise),1),arModel.covNoise,nPoints);
    else
        nDim = length(arModel.coef);
        currentHistory = mvnrnd(zeros(nDim,1), arModel.covNoise, arModel.nLags)';
        
        ts = zeros(nPoints, nDim);
        indNoiseSeries = mvnrnd(zeros(nDim,1),arModel.covEps,nPoints);
        for n=1:nPoints
            for d=1:nDim
                ts(n,d) = sum(sum(arModel.coef{d} .* currentHistory(:))) + indNoiseSeries(n,d);
            end
            currentHistory = [ts(n,:)', currentHistory(:,1:(end-1))];
        end
    end
end

