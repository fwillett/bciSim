function [ features ] = cpwlDesignMatrix( breaks, xData )
    %make a design matrix for a least squares regression fit of a
    %continuous piecewise linear function
    [~,binIdx]=histc(xData, breaks);
    features = zeros(length(xData), length(breaks));
    for segIdx=1:(length(breaks)-1)
        dataIdx = (binIdx==segIdx);
        features(dataIdx,segIdx)=(-xData(dataIdx)+breaks(segIdx+1))/(breaks(segIdx+1)-breaks(segIdx));
        features(dataIdx,segIdx+1)=(xData(dataIdx)-breaks(segIdx))/(breaks(segIdx+1)-breaks(segIdx));
    end
end

