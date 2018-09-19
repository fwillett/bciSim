%%
%This script gives an example of how to re-parameterize a steady-state Kalman filter to
%the (alpha, beta, D) parameterization. A Kalman filter is calibrated on
%simulated open-loop data, then reparameterized. 

%%
%Parameters
nNeurons = 40; %number of "neurons" to include
neuronSNR = 0.2; %higher values mean the "neurons" have more movement information
movTime = 0.75; %in seconds, how long the simulated minimum-jerk trajectories last
nRep = 10; %how many blocks of 8 center-out trajectories to include

%%
%First, make some simulated opne-loop data. The simulated neurons will
%linearly encode a control vector during a set of minimum-jerk,
%center-out trajectories.

%generate 8 center-out targets
theta = linspace(0,2*pi,9)';
targList =  [cos(theta(1:8)), sin(theta(1:8))];

%generate minimum jerk trajectories
XY = [];
TargXY = [];
VelXY = [];

for s=1:size(targList,1)
    tAxis = linspace(0,movTime,round(movTime*50));
    minJerkPos = [targList(s,1)*(10*(tAxis/movTime).^3-15*(tAxis/movTime).^4+6*(tAxis/movTime).^5)
        targList(s,2)*(10*(tAxis/movTime).^3-15*(tAxis/movTime).^4+6*(tAxis/movTime).^5)]';
    
    XY = [XY; minJerkPos];
    TargXY = [TargXY; repmat(targList(s,:), length(minJerkPos), 1)];
    VelXY = [VelXY; [0 0; diff(minJerkPos)*50]];
end

%repeat these center-out trajectories nRep times
XY = repmat(XY, nRep, 1);
TargXY = repmat(TargXY, nRep, 1);
VelXY = repmat(VelXY, nRep, 1);

%Have the "neurons" encode a simple linear control policy (they will encode
%TargetPosition - CursorPosition). 
controlVector = TargXY - XY;

%Generate tuning coefficients randomly.
E = randn(2,40)*0.2;

%Simulate the neural activity with Gaussian noise.
simNeural = controlVector*E + randn(size(controlVector,1), 40)*2;

%%
%Make the Kalman filter with MATLAB's "kalman" function.

%first, fit the state transition model and observation model with least
%squares. A is the satate transition matrix, Q is the state transition
%covariance, C is the observation matrix and R is the observation
%covariance.
A = controlVector(1:(end-1),:) \ controlVector(2:end,:);
Q = cov(controlVector(2:end,:)-controlVector(1:(end-1),:)*A);

C = (controlVector \ simNeural)'; 
R = cov(simNeural-controlVector*C');

%Now get Kalman gain matrix.
sys = ss(A,[zeros(2) eye(2)],C,[],0.02);
[k,L,P,M,Z] = kalman(sys,Q,R,[]);
kalmanGain = M;            

%Apply the Kalman filter to make sure it works.
decVel = zeros(size(controlVector));
for t=2:size(controlVector,1)
    decVel(t,:) = A*decVel(t-1,:)' + kalmanGain*(simNeural(t,:)'-C*A*decVel(t-1,:)');
end

figure('Position',[624   243   615   735]);
ax1 = subplot(2,1,1);
hold on
plot(decVel(:,1));
plot(controlVector(:,1));
legend({'Kalman Filter Output (X Dim)','Encoded Control Vector (X Dim)'});
xlabel('Time Step');

ax2 = subplot(2,1,2);
hold on
plot(decVel(:,2));
plot(controlVector(:,2));
legend({'Kalman Filter Output (Y Dim)','Encoded Control Vector (Y Dim)'});
xlabel('Time Step');

linkaxes([ax1, ax2]);

%Compute correlation between decoded velocity and the encoded control
%vector.
disp(['X Velocity Correlation: ' num2str(corr(decVel(:,1), controlVector(:,1)))]); 
disp(['Y Velocity Correlation: ' num2str(corr(decVel(:,2), controlVector(:,2)))]); 

%%
%Reparameterize the Kalman filter into (alpha, beta, D) format so we can understand and
%report its gain (beta) and smoothing (alpha) properties.
[ alpha, beta, D ] = reparamKalman( kalmanGain, A, C, TargXY - XY, simNeural, [0.8 1]);

%Prove that this reparameterization is nearly equivalent (won't be exactly
%equivalent, since a matrix multiplication is approximated by a scalar
%multiplication by alpha).
decVelReparam = zeros(size(controlVector));
for t=2:size(controlVector,1)
    decVelReparam(t,:) = alpha*decVelReparam(t-1,:) + (1-alpha)*beta*(simNeural(t,:)*D);
end

%Plot reparameterized decoder time series vs. Kalman filter time series
figure('Position',[624   243   615   735]);
ax1 = subplot(2,1,1);
hold on
plot(decVel(:,1));
plot(decVelReparam(:,1));
legend({'Kalman Filter Output (X Dim)','Reparameterized Filter Output (X Dim)'});
xlabel('Time Step');

ax2 = subplot(2,1,2);
hold on
plot(decVel(:,2));
plot(decVelReparam(:,2));
legend({'Kalman Filter Output (Y Dim)','Reparameterized Filter Output (Y Dim)'});
xlabel('Time Step');

linkaxes([ax1, ax2]);

%Plot reparameterized decoder vs. Kalman filter scatter plot
figure;
ax1 = subplot(2,1,1);
hold on;
plot(decVel(:,1),decVelReparam(:,1),'.');
title('Decoder Outputs (X)');
xlabel('Kalman Output');
ylabel('Reparameterized Linear Output');
axis equal;
plot([-1,1],[-1,1],'--k');

ax2 = subplot(2,1,2);
hold on;
plot(decVel(:,2),decVelReparam(:,2),'.');
title('Decoder Outputs (Y)');
xlabel('Kalman Output');
ylabel('Reparameterized Linear Output');
axis equal;
plot([-1,1],[-1,1],'--k');

linkaxes([ax1, ax2]);

