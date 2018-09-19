%%
%This example script demonstrates how to use the simulator itself. The .mex file
%should be compiled first (compileSimBci.m).

%%
%This section runs an alpha/beta parameters sweep to find optimal gain and
%smoothing parameters.

%This function fills a simulator options struct with defaults.
opts = makeBciSimOptions( );

%For this example, set some of the parameters to different values than the default.
opts.trial.dwellTime = 1.0;
opts.noiseMatrix = randn(10000,2)*1.5;
opts.forwardModel.delaySteps = 10;
opts.forwardModel.forwardSteps = 10;

%save original options
originalOpts = opts;

%Specify the alpha and beta values to test.
alpha = fliplr(1-logspace(log10(0.005),log10(0.8),15));
beta = logspace(log10(0.3),log10(6.25),20);

%At each alpha and beta value, we do a sweep of fVel slopes and pick
%the best-performing one to report (this makes the assumption that the user
%will adapt their fVel).
velSlopes = linspace(0,-2,10);

timeMat = zeros(length(alpha),length(beta),length(velSlopes));
nTrials = 50;
allResults = cell(length(alpha), length(beta));

for a=1:length(alpha)
    disp([num2str(a) ' / ' num2str(length(alpha))]);
    for b=1:length(beta)
        tmpResults = cell(length(velSlopes),1);
        for v=1:length(velSlopes)
            %specify the plant and control policy
            opts.plant.alpha = alpha(a);
            opts.plant.beta = beta(b);
            opts.control.fVelX = [0 1];
            opts.control.fVelY = [0 velSlopes(v)];
            
            %simulate a batch of movements
            startPos = repmat([0 0], nTrials, 1);
            targPos = repmat([1 0], nTrials, 1);
            tmpResults{v} = simBatch( opts, targPos, startPos );
            timeMat(a,b,v) = mean(tmpResults{v}.movTime);
        end
        
        [~,bestIdx] = min(squeeze(timeMat(a,b,:)));
        allResults{a,b} = tmpResults{bestIdx};
    end
end

%For these values,  we will plot example trajectories.
aToPlot = length(alpha):-3:1;
bToPlot = 1:3:length(beta);

%First, summarize the results with a heatmap figure.
bLabels = cell(length(beta),1);
for b=1:length(bLabels)
    bLabels{b} = num2str(beta(b),2);
end

aLabels = cell(length(alpha),1);
for b=1:length(aLabels)
    aLabels{b} = num2str(alpha(b),2);
end

figure
hold on
imagesc(squeeze(min(timeMat,[],3)),[0 10]); 
for aIdx = 1:length(aToPlot)
    for bIdx = 1:length(bToPlot)
        plot(bToPlot(bIdx), aToPlot(aIdx), 'ko');
    end
end
colormap(jet);
set(gca,'YDir','normal');
set(gca,'XTick',1:3:length(beta),'XTickLabel',bLabels(1:3:end));
set(gca,'YTick',1:3:length(alpha),'YTickLabel',aLabels(1:3:end));
xlabel('Beta (TD/s)');
ylabel('Alpha');
title('Average Movement Time (s)');
colorbar;
xlim([1 length(beta)]);
ylim([1 length(alpha)]);

%Then, plot example trajectories for selected alpha/beta settings (marked
%with black circles on the heat map).
figure('Position',[75         224        1152         881]);
for aIdx = 1:length(aToPlot)
    for bIdx = 1:length(bToPlot)
        subplot(length(aToPlot), length(bToPlot), (aIdx-1)*length(bToPlot) + bIdx);
        hold on
        resultsToPlot = allResults{aToPlot(aIdx), bToPlot(bIdx)};
        for t=1:length(resultsToPlot.movTime)
            loopIdx = resultsToPlot.reachEpochs(t,1):resultsToPlot.reachEpochs(t,2);
            plot(resultsToPlot.pos(loopIdx,1), resultsToPlot.pos(loopIdx,2), 'LineWidth',1);
        end
        rectangle('Position',[1-opts.trial.targRad, 0-opts.trial.targRad, opts.trial.targRad*2, opts.trial.targRad*2],...
            'Curvature',[1 1], 'LineWidth', 1, 'EdgeColor','k');
        axis equal;
        xlim([-0.2 1.4]);
        ylim([-0.8 0.8]);
        title(['\alpha=' num2str(alpha(aToPlot(aIdx)),2) ', \beta=' num2str(beta(bToPlot(bIdx)),2)]);
        set(gca,'XTick',[],'YTick',[]);
        axis off;
    end
end

%This function can also be used to perform an alpha & beta parameter sweep:
%alphaBetaSweep(originalOpts, alpha, beta, 1, 0.2, 1);
%%
%This section simulates the average movement time as a function of dwell time, target
%radius, noise STD, and feedback delay.
nTrials = 250;
startPos = repmat([0 0], nTrials, 1);
targPos = repmat([1 0], nTrials, 1);

meanMovTime = zeros(10,1);            
dwellTime = linspace(0.1,5,10);
targRad = linspace(0.05,0.25,10);
noiseStd = linspace(0.2,3,10);
visDelay = round(linspace(0,20,10));

defaultOpts = makeBciSimOptions( );
defaultOpts.trial.maxTrialTime = 15;
defaultOpts.trial.dwellTime = 1;
defaultOpts.trial.continuousHoldRule = 1;
defaultOpts.noiseMatrix = randn(1000000,2)*1.5;

figure
subplot(2,2,1);
opts = defaultOpts;
for x = 1:length(dwellTime)
    opts.trial.dwellTime = dwellTime(x);
    out = simBatch( opts, targPos, startPos );
    meanMovTime(x) = mean(out.movTime);
end
plot(dwellTime, meanMovTime, '-o');
xlabel('Dwell Time (s)');
ylabel('Mean Movement Time (s)');

subplot(2,2,2);
opts = defaultOpts;
for x = 1:length(targRad)
    opts.trial.targRad = targRad(x);
    out = simBatch( opts, targPos, startPos );
    meanMovTime(x) = mean(out.movTime);
end
plot(targRad, meanMovTime, '-o');
xlabel('Target Radius');
ylabel('Mean Movement Time (s)');

subplot(2,2,3);
opts = defaultOpts;
for x = 1:length(noiseStd)
    opts.noiseMatrix = randn(10000,2)*noiseStd(x);
    out = simBatch( opts, targPos, startPos );
    meanMovTime(x) = mean(out.movTime);
end
plot(noiseStd, meanMovTime, '-o');
xlabel('Noise STD');
ylabel('Mean Movement Time (s)');

subplot(2,2,4);
opts = defaultOpts;
for x = 1:length(visDelay)
    opts.forwardModel.delaySteps = visDelay(x);
    opts.forwardModel.forwardSteps = visDelay(x);
    out = simBatch( opts, targPos, startPos );
    meanMovTime(x) = mean(out.movTime);
end
plot(visDelay, meanMovTime, '-o');
xlabel('Feedback Delay (# of steps)');
ylabel('Mean Movement Time (s)');

