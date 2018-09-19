function [optAlpha, optBeta, simOut] = alphaBetaSweep(simOpts, alpha, beta, targDist, targRad, dwellTime)
    %[optAlpha, optBeta, simOut] = alphaBetaSweep(simOpts, alpha, beta, targDist, targRad, dwellTime)    
    %finds the optimal gain and smoothing value (alpha & beta) for the specified simulation settings. 
    %
    %simOpts is a simulation options struct (makeBciSimOptions defines the
    %default values and fields).
    %
    %alpha is a vector of alpha parameters to test
    %
    %beta is a vector of beta parameters to test
    %
    %targDist is a scalar value describing the target distance
    %
    %targRad is a scalar value describing the target radius
    %
    %dwellTime is a scalar value describing the dwell time (in seconds)
    
    simOpts.trial.targRad = targRad;
    simOpts.trial.dwellTime = dwellTime;
    
    %At each alpha and beta value, we do a sweep of fVel slopes and pick
    %the best-performing one to report (this makes the assumption that the user
    %will adapt their fVel). 
    velSlopes = linspace(0,-2,10);

    out = cell(length(alpha),length(beta),length(velSlopes));
    timeMat = zeros(length(alpha),length(beta));
    nTrials = 50;

    for a=1:length(alpha)
        disp([num2str(a) ' / ' num2str(length(alpha))]);
        for b=1:length(beta)
            tmpOut = cell(length(velSlopes),1);
            for v=1:length(velSlopes)
                %specify the plant and control policy
                simOpts.plant.alpha = alpha(a);
                simOpts.plant.beta = beta(b);
                simOpts.control.fVelX = [0 100];
                simOpts.control.fVelY = [0 velSlopes(v)*100];

                %simulate a batch of movements
                startPos = repmat([0 0], nTrials, 1);
                targPos = repmat([targDist 0], nTrials, 1);
                [ tmpOut{v} ] = simBatch( simOpts, targPos, startPos );
                timeMat(a,b,v) = mean(tmpOut{v}.movTime);
            end
            [~, minIdx] = min(timeMat(a,b,:));
            out{a,b} = tmpOut{minIdx};
        end
    end

    %%
    %return and plot results
    optTimes = squeeze(min(timeMat,[],3));
    [~,minIdx] = min(optTimes(:));
    [alphaIdx, betaIdx] = ind2sub(size(optTimes), minIdx);
    optAlpha = alpha(alphaIdx);
    optBeta = beta(betaIdx);
    simOut = out{alphaIdx, betaIdx};
    
    bLabels = cell(length(beta),1);
    for b=1:length(bLabels)
        bLabels{b} = num2str(beta(b)/targDist,2);
    end

    aLabels = cell(length(alpha),1);
    for b=1:length(aLabels)
        aLabels{b} = num2str(alpha(b),2);
    end

    figure('Position',[36         108        1104         331]);
    subplot(1,2,1);
    hold on
    imagesc(optTimes,[0 10]); 
    plot(betaIdx, alphaIdx, 'wx', 'LineWidth', 2, 'MarkerSize', 10);
    xlim([1 length(beta)]);
    ylim([1 length(alpha)]);
    colormap(jet);
    set(gca,'YDir','normal');
    set(gca,'XTick',1:3:length(beta),'XTickLabel',bLabels(1:3:end));
    set(gca,'YTick',1:3:length(alpha),'YTickLabel',aLabels(1:3:end));
    xlabel('Beta (TD/s)');
    ylabel('Alpha');
    title('Average Movement Time (s)');
    colorbar;
    
    subplot(1,2,2);
    hold on
    for t=1:size(simOut.reachEpochs,1)
        loopIdx = simOut.reachEpochs(t,1):simOut.reachEpochs(t,2);
        plot(simOut.pos(loopIdx,1), simOut.pos(loopIdx,2));
    end
    halfRad = simOpts.trial.targRad / 2;
    rectangle('Position', [simOut.targPos(1,1)-halfRad, simOut.targPos(1,2)-halfRad, halfRad*2, halfRad*2], ...
        'Curvature', [1 1], 'LineWidth', 2);
    axis equal;
    title('Example Simulated Trajectories');
end
