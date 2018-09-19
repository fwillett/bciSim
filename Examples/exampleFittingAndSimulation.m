%This script fits the piecewise linear model to the slowest gain condition,
%and then uses it to predict what performance will be like on all other
%conditions.

%The example dataset is from participant T8's "high gain" session, and this
%script produces a figure similar to Figs 3 and 4. In this session,
%participant T8 completed 15 blocks of closed-loop 2D cursor control on a
%center-out-back task under 1 of 4 possible decoder gains. 

%If your operating system is not 64 bit Windows, you will have to compile
%the mex simulation function (you can use the script compileSimBci.m).

%The README.txt file explains the format for the sample dataset.

%%
%Load the sample dataset and specify fitting parameters.
fileName = 'T8.2015.03.24.mat';
feedbackDelayTimeSteps = round(0.4/0.02);
reactionTimeSteps = round(0.4/0.02);
farDistance = 14;
data = load(fileName);

%Prepare to fit model on the slowest condition by getting the "reachEpochs"
%(list of start and end indices for each trial that occurred during a slow
%block). 
[~,conFit] = min(data.beta); 
conditionNumberByTrial = data.conditionNumber(data.trialEpochs(:,1));
reachEpochs = data.trialEpochs(conditionNumberByTrial==conFit,:);
reachEpochs(:,1) = reachEpochs(:,1) + reactionTimeSteps;

%Create a data fitting options struct to pass into fitPiecewiseModel. See
%makePiecewiseModelOptions() for explanations of each parameter.
opts.pos = data.cursorPos;
opts.vel = data.cursorVel;
opts.targPos = data.targetPos;
opts.decoded_u = data.decodedControlVector;

opts.modelOpts.noVel = false;
opts.modelOpts.nKnots = 12;
opts.modelOpts.noNegativeFTarg = true;

opts.filtAlpha = data.alpha(conFit);
opts.filtBeta = data.beta(conFit);

opts.reachEpochsToFit = reachEpochs;
opts.feedbackDelaySteps = feedbackDelayTimeSteps;
opts.timeStep = 0.02;
opts.fitNoiseModel = true;
opts.fitSDN = true;

%Fit the PLM model.
disp(['Fitting model to condition ' num2str(conFit) ' (may take a minute)']);
model = fitPiecewiseModel( opts );

%Now use the model to simulate performance on all conditions.
perfStats = cell(2,4,length(data.alpha));
simResults = cell(length(data.alpha),1);

velCoef = linspace(0,-1,15);
for c=1:length(data.alpha)
    disp(['Simulating condition ' num2str(c)]);
    
    %Prepare simulation input.
    simOpts = model.simOpts;
    simOpts.plant.alpha = data.alpha(c);
    simOpts.plant.beta = data.beta(c);
    simOpts.trial.dwellTime = data.dwellTimes(c);
    simOpts.trial.targRad = data.cursorRadius + data.targetRadius;

    targAngles = linspace(0,2*pi,9);
    targAngles = targAngles(1:8);
    targPos = [cos(targAngles)', sin(targAngles)']*farDistance;
    targPos = repmat(targPos, 62, 1);
    startPos = repmat([0 0], size(targPos,1), 1);

    %This part of the code implements a simple model of user adaptation by
    %iterating over possibilities for how the user will adapt the
    %velocity damping component of their control policy. We assume the user
    %will choose whatever velocity damping coefficient works best. 
    avgMovTime = zeros(length(velCoef),1);
    for v=1:length(velCoef)
        simOpts.control.fVelX = [0 farDistance];
        simOpts.control.fVelY = [0, velCoef(v)];
        tmp = simBatch( simOpts, targPos, startPos );
        avgMovTime(v) = mean(tmp.movTime);
    end

    [~,minIdx] = min(avgMovTime);
    simOpts.control.fVelX = [0 farDistance];
    simOpts.control.fVelY = [0, velCoef(minIdx)];
    
    %Simulate using the best velocity damping coefficient.
    simResults{c} = simBatch( simOpts, targPos, startPos );

    %Compute performance summary statistics
    trlIdx = conditionNumberByTrial==c;
    [perfStats{1,1,c}, perfStats{1,2,c}, perfStats{1,3,c}, perfStats{1,4,c}] = ...
        trajectoryPerformance(data.cursorPos, data.targetPos, data.cursorRadius + data.targetRadius, data.trialEpochs(trlIdx,:), 0.02);
    [perfStats{2,1,c}, perfStats{2,2,c}, perfStats{2,3,c}, perfStats{2,4,c}] = ...
        trajectoryPerformance(simResults{c}.pos, simResults{c}.targPos, data.cursorRadius + data.targetRadius, simResults{c}.reachEpochs, 0.02);
end

%%
%Plot actual performance vs. predicted performance.
[conList, ~, conIdx] = unique(round(data.beta/0.1)*0.1);
yLimits = [0 1; 
    0 2.2;
    0 3;
    0 1.2];
yLabels = {'Dial-in Time (s)','Translation Time (s)','Total Movement Time (s)','Path Efficiency'};

figure('Position',[49         174        1313         343]);
for p=1:4
    realLine = zeros(length(conList),3);
    for c=1:length(conList)
        blockIdx = find(conIdx==c);
        tmp = [];
        for b=1:length(blockIdx)
            tmp = [tmp; perfStats{1,p,blockIdx(b)}];
        end
        
        notNan = ~isnan(tmp);
        [realLine(c,1),~,realLine(c,2:3)] = normfit(tmp(notNan));
    end

    simLine = zeros(length(conList),3);
    for c=1:length(conList)
        blockIdx = find(conIdx==c);
        tmp = [];
        for b=1:length(blockIdx)
            tmp = [tmp; perfStats{2,p,blockIdx(b)}];
        end
        
        notNan = ~isnan(tmp);
        [simLine(c,1),~,simLine(c,2:3)] = normfit(tmp(notNan));
    end
    
    subplot(1,4,p);
    hold on;
    plot(conList/farDistance, realLine(:,1), '-ko', 'LineWidth', 2);
    plot(conList/farDistance, simLine(:,1), 'r', 'LineWidth', 2);
    errorbar(conList/farDistance, realLine(:,1), realLine(:,1)-realLine(:,2), realLine(:,3)-realLine(:,2), 'ko', 'LineWidth', 2);
    
    ylim(yLimits(p,:));
    xlabel('Gain (TD/s)');
    ylabel(yLabels{p});
    
    if p==1
        legend({'Observed','Simulated'},'Location','NorthWest');
    end
end

%%
%Plot actual vs. simulated trajectories.
conToPlot = [4 5 6 7];
targAngles = targAngles(1:8);
targList = [cos(targAngles)', sin(targAngles)']*farDistance;
targList = [targList; 0 0];

simTargNum = repmat((1:8)', 62, 1);
realTargNum = zeros(size(data.trialEpochs,1),1);
for t=1:size(data.trialEpochs,1)
    tmp = repmat(data.targetPos(data.trialEpochs(t,1),:)-[0, 40.5],9,1) - targList;
    tmp = sqrt(sum(tmp.^2,2));
    [~,realTargNum(t)] = min(tmp);
end
isOuter = realTargNum~=9;
targColors = hsv(8)*0.8;
targRad = data.cursorRadius + data.targetRadius;

figure('Position',[ 31          88        1382         562]);
for c=1:length(conToPlot)
    subplot(2,4,c);
    hold on
    trlIdx = find(conditionNumberByTrial==conToPlot(c));
    for t=1:length(trlIdx)
        if isOuter(t)
            loopIdx = data.trialEpochs(trlIdx(t),1):data.trialEpochs(trlIdx(t),2);
            plot(data.cursorPos(loopIdx,1), data.cursorPos(loopIdx,2), 'LineWidth',1,'Color',targColors(realTargNum(trlIdx(t)),:));
        end
    end
    for t=1:size(targList,1)
        rectangle('Position',[targList(t,1)-targRad, 40.5+targList(t,2)-targRad, targRad*2, targRad*2],...
            'Curvature',[1 1], 'LineWidth', 1, 'EdgeColor','k');
    end
    xlim([-18 18]);
    ylim([-18 18]+40.5);
    axis equal;
    title(['Gain = ' num2str(data.beta(conToPlot(c))/farDistance,3)]);
    if c==1
        ylabel('Observed');
    end
    set(gca,'XTick',[],'YTick',[]);
    
    subplot(2,4,c+4);
    hold on
    for t=1:round(length(trlIdx)/2)
        loopIdx = simResults{conToPlot(c)}.reachEpochs(t,1):simResults{conToPlot(c)}.reachEpochs(t,2);
        plot(simResults{conToPlot(c)}.pos(loopIdx,1), simResults{conToPlot(c)}.pos(loopIdx,2), 'LineWidth',1,'Color',targColors(simTargNum(t),:));
    end
    for t=1:size(targList,1)
        rectangle('Position',[targList(t,1)-targRad, targList(t,2)-targRad, targRad*2, targRad*2],...
            'Curvature',[1 1], 'LineWidth', 1, 'EdgeColor','k');
    end
    xlim([-18 18]);
    ylim([-18 18]);
    axis equal;
    if c==1
        ylabel('Simulated');
    end
    set(gca,'XTick',[],'YTick',[]);
end
