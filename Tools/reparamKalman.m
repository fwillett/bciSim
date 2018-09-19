function [ alpha, beta, D ] = reparamKalman( K, A, H, posErr, featureVectors, farDistInterval)
    %reparamKalman( K, A, H, posErr, featureVectors, farDistInterval) converts from Kalman filter matrices to the (alpha, beta, D) parameterization.
    %
    %K is a D x N Kalman gain matrix (D = number of dimensions, N = number
    %of neural features).
    %
    %A is a D x D state transition matrix.
    %
    %H is an N x D observation model matrix.
    %
    %posErr is a T x D matrix of position errors (T = number of time samples), where each row is a
    %position error vector (TargetPos-CursorPos).
    %
    %featureVectors is a T x N matrix feature vectors corresponding to the
    %position errors in posErr.
    %
    %farDistInterval is a 1 x 2 vector describing a distance window that is
    %considered "far" from the target, that is used for decoder
    %normalization. The first entry is the beginning of the window and the
    %second entry is the end of the window. 
    aMat = (eye(length(A))-K*H)*A;
    alpha = mean(diag(aMat));

    decVectors = featureVectors * K';
    normFactor = normalizeDecoder(posErr, decVectors/(1-alpha), farDistInterval);
    D = normFactor*K'/(1-alpha);
    beta = 1 / normFactor;
end

