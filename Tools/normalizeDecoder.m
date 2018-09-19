function [ normFactor ] = normalizeDecoder(posErr, decVectors, farDistInterval)
    %normFactor = normalizeDecoder(posErr, decVectors, farDistInterval)
    %returns a normFactor which is what you need to multiply your decoding matrix by in
    %order to normalize it so that the mean of its projected output is 1. 
    
    %posErr is a T x D matrix of position error vector samples (T = number
    %of samples, D = number of dimensions).
    
    %decVectors is a T x D matrix of decoded vectors
    
    %farDistInterval is a distance interval that is considered to be "far"
    %from the target. This function will only normalize based on data where
    %the cursor is "far" from the target.
    
    targDist = sqrt(sum(posErr.^2,2));
    unitErr = bsxfun(@times, posErr, 1./targDist);
    projDec = sum(unitErr.*decVectors,2);
    
    farIdx = (targDist > farDistInterval(1)) & (targDist < farDistInterval(2));
    normFactor = 1/mean(projDec(farIdx));
end
