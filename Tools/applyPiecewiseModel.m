function [ cVec ] = applyPiecewiseModel( model, posHat, velHat, targPos )
    %cVec = applyPiecewiseModel( model, posHat, velHat, targPos ) uses a
    %model that was fit with fitPiecewiseModel() and applies it to a given
    %dataset to simulate the user's control vector.
    %
    %model is a struct from fitPiecewiseModel().
    %
    %posHat is a T x D matirx of internal model cursor positions (where T =
    %number of samples and D = number of dimensions)
    %
    %velHat is a T x D matrix of internal model cursor velocities
    %
    %targPos is a T x D matrix of target positions
    
    dist = matVecMag(targPos - posHat, 2);
    speed = matVecMag(velHat, 2);
    
    toTargVec = bsxfun(@times, (targPos - posHat), 1./dist);
    velVec = bsxfun(@times, velHat, 1./speed);
    toTargVec(dist==0,:) = 0;
    velVec(speed==0,:)=0;

    distWeight = interp1([-0.01 model.fTargX 10*model.fTargX(end)], [model.fTargY(1) model.fTargY model.fTargY(end)], dist,'linear','extrap');
    distWeight(distWeight<0) = 0;
    
    cVec = zeros(size(posHat));
    cVec = cVec + bsxfun(@times, toTargVec, distWeight);
    cVec = bsxfun(@plus, cVec, model.bias(1,:));
    
    if ~isempty(model.fVelX)
        speedWeight = interp1([-0.01 model.fVelX 10*model.fVelX(end)], [model.fVelY(1) model.fVelY model.fVelY(end)], speed,'linear','extrap');
        speedWeight(speedWeight>0) = 0;
        cVec = cVec + bsxfun(@times, velVec, speedWeight);
    end
end

function [ mag ] = matVecMag( mat, dim )
    if length(size(mat))>2
        error('Cannot use on matrices bigger than 2 dimensions');
    end
    mag = sqrt(sum(mat.^2,dim));
end
