%% Setup: run simulation
load simulationSettings
load orbitGeometry
load spacecraftProperties

% Set up initial conditions - arbitrary for testing KF
[r0, v0] = oe2mci(a,e,i,raan,argp,trueAnomaly0);   % initial position w.r.t. inertial frame
w0_principal = [0,0,0]';
q0 = [0, 0, 0, 1];
y0 = [r0', v0', q0, w0_principal'];               % initial state vector
disp("Initial conditions ready!")

% Run "ground truth" simulation
[t, y_out] = ode113(@trueOrbitAndAttitude, [0:tstep:t_period]' , y0, options);
pos = y_out(:,1:3); vel = y_out(:,4:6); q = y_out(:,7:10); w = y_out(:,11:13);
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
errorTheta = zeros(size(theta)); errorPhi = zeros(size(theta)); errorPsi = zeros(size(theta));
disp("Simulation done!");

%% Run Kalman Filter with coarse attitude sensing suite
% Initial conditions
variances = [9296.8361, 5622.6309, 21349.0534, 1504.4700, .000068567, .0000403, .0012959]*10^(-6);
P0 = diag(variances).^2; Q = (1/100) * P0;
x0 = [q0'; w0_principal];

x_est = zeros(7,length(t));
P_est = zeros(7,7,length(t));
z_pre_gyro = zeros(3,length(t));
z_post_gyro = zeros(3,length(t));
z_pre_horizon = zeros(3,length(t));
z_post_horizon = zeros(3,length(t));
z_pre_sun = zeros(3, length(t));
z_post_sun = zeros(3, length(t));

x_est(:,1) = x0;
P_est(:,:,1) = P0;
x = x0; P = P0;

% Run KF
disp("KF started")
for idx = 2:length(t)
    
    % Predict step
    [x,P] = KF_TimeUpdate(x,P,Q,0,t(idx)-t(idx-1));

    % Model measurements using estimated state and no noise
    wModel = gyroMeasurement(x(5:7)', false);               
    [nadirModel, nadirRef] = horizonSensorMeasurement(x(1:4)', pos(idx,:), false);
    [sunModel, sunRef] = sunSensorMeasurement(x(1:4)', t(idx), false);

    % Sequential update step - gyro
    wMeas = gyroMeasurement(w(idx,:), true); % actual measurement
    H = sensitivityMatrix("gyro", q(idx,:), 0);
    R = measurementCovarianceMatrix("gyro");
    [x, P, K] = KF_MeasurementUpdate(x,wMeas',wModel',P,H,R);

    % Sequential update step - horizon sensor
    [nadirMeas, nadirRef] = horizonSensorMeasurement(q(idx,:), pos(idx,:), true);  % actual measurement
    H = sensitivityMatrix("horizon", q(idx,:), nadirRef);
    R = measurementCovarianceMatrix("horizon");
    [x, P, K] = KF_MeasurementUpdate(x,nadirMeas, nadirModel, P,H,R);

    % Sequential update step - sun sensor
    if (eclipseCheck(pos(idx,:), t(idx)) == false)
        [sunMeas, sunRef] = sunSensorMeasurement(q(idx,:), t(idx), true); % actual measurement
        H = sensitivityMatrix("sun sensor", q(idx,:), sunRef);
        R = measurementCovarianceMatrix("sun sensor");
        z_pre_sun(:,idx) = sunModel - sunMeas;
        [x, P, K] = KF_MeasurementUpdate(x,sunMeas, sunModel,P,H,R);
        [sunModelPost, sunRefPost] = sunSensorMeasurement(x(1:4)', t(idx), false); 
        z_post_sun(:,idx) = sunModelPost - sunMeas;

    end

    x_est(:,idx) = x;
    P_est(:,:,idx) = P;

    % Pre- and post-fit residuals
    wModelPost = gyroMeasurement(x(5:7)', false);
    [nadirModelPost, nadirRefPost] = horizonSensorMeasurement(x(1:4)',pos(idx,:), false);
    z_pre_gyro(:,idx) = wModel - wMeas;
    z_pre_horizon(:,idx) = nadirModel - nadirMeas;
    z_post_gyro(:,idx) = wModelPost - wMeas;
    z_post_horizon(:,idx) = nadirModelPost - nadirMeas;
end

x_est = x_est(:,1:length(t)); P_est = P_est(:,:,1:length(t));

Pq1 = reshape(P_est(1,1,:), [1,length(t)]); 
Pq2 = reshape(P_est(2,2,:), [1,length(t)]);
Pq3 = reshape(P_est(3,3,:), [1,length(t)]);
Pq4 = reshape(P_est(4,4,:), [1,length(t)]);
Pwx = reshape(P_est(5,5,:), [1,length(t)]);
Pwy = reshape(P_est(6,6,:), [1,length(t)]);
Pwz = reshape(P_est(7,7,:), [1,length(t)]);

disp("KF Done (Coarse Sensing)")


pltTitle = sprintf("Kalman Filter: Quaternion Attitude (Coarse Sensing Suite)");
figure(); sgtitle(pltTitle);
subplot(4,3,1); hold on; title('q1');
plot(t, q(:,1)); plot(t, x_est(1,:));
grid on; hold off; legend('True', 'Estimated');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,2); hold on; title('q1 Error');
plot(t, x_est(1,:)'-q(:,1));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,3); hold on; title('q1 Covariance');
plot(t, sqrt(Pq1));
grid on; hold off; 
xlabel('Time (s)'); ylabel('\sigma');
subplot(4,3,4); hold on; title('q2');
plot(t, q(:,2)); plot(t, x_est(2,:));
grid on; hold off; legend('True', 'Estimated');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,5); hold on; title('q2 Error');
plot(t, x_est(2,:)'-q(:,2));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,6); hold on; title('q2 Covariance');
plot(t, sqrt(Pq2));
grid on; hold off; 
xlabel('Time (s)'); ylabel('\sigma');
subplot(4,3,7); hold on; title('q3');
plot(t, q(:,3)); plot(t, x_est(3,:));
grid on; hold off; legend('True', 'Estimated');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,8); hold on; title('q3 Error');
plot(t, x_est(3,:)'-q(:,3));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,9); hold on; title('q3 Covariance');
plot(t, sqrt(Pq3));
grid on; hold off; 
xlabel('Time (s)'); ylabel('\sigma');
subplot(4,3,10); hold on; title('q4');
plot(t, q(:,4)); plot(t, x_est(4,:));
grid on; hold off; legend('True', 'Estimated');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,11); hold on; title('q4 Error');
plot(t, x_est(4,:)'-q(:,4));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,12); hold on; title('q4 Covariance');
plot(t, sqrt(Pq4));
grid on; hold off; 
xlabel('Time (s)'); ylabel('\sigma');

pltTitle = sprintf("EKF Time Update\n Angular Velocity (Coarse Sensing Suite)");
figure(); sgtitle(pltTitle);
subplot(2,2,1); hold on; title('True Angular Velocity');
plot(t, w(:,1)); plot(t, w(:,2)); plot(t, w(:,3));
grid on; hold off;
xlabel('Time (s)'); ylabel('Rad/s'); legend('\omega_x', '\omega_y', '\omega_z');
subplot(2,2,2); hold on; title('Estimated Angular Velocity');
plot(t, x_est(5,:)); plot(t, x_est(6,:)); plot(t, x_est(7,:));
grid on; hold off;
xlabel('Time (s)'); ylabel('Rad/s'); legend('\omega_x', '\omega_y', '\omega_z');
subplot(2,2,3); hold on; title('Angular Velocity Estimation Error');
plot(t, x_est(5,:)'-w(:,1)); plot(t, x_est(6,:)'-w(:,2)); plot(t, x_est(7,:)'-w(:,3));
grid on; hold off;
xlabel('Time (s)'); ylabel('Rad/s'); legend('\omega_x', '\omega_y', '\omega_z');
subplot(2,2,4); hold on; title('Covariances');
plot(t, sqrt(Pwx)); plot(t, sqrt(Pwy)); plot(t, sqrt(Pwz));
grid on; hold off;
xlabel('Time (s)'); ylabel('\sigma'); legend('\omega_x', '\omega_y', '\omega_z');

figure();
subplot(3,2,1); title('Sun Sensor'); hold on; grid on;
plot(t, vecnorm(z_pre_sun)); plot(t,vecnorm(z_post_sun)); legend('Pre-Fit', 'Post-Fit');
xlabel('Time (s)'); ylabel('Residual Magnitude');
subplot(3,2,2); title('Sun Sensor'); hold on; grid on;
plot(t, vecnorm(z_pre_sun)-vecnorm(z_post_sun)); legend('z_{pre} - z_{post}');
xlabel('Time (s)'); ylabel('Residual Magnitude');
subplot(3,2,3); title('Horizon Sensor'); hold on; grid on;
plot(t, vecnorm(z_pre_horizon)); plot(t,vecnorm(z_post_horizon)); legend('Pre-Fit', 'Post-Fit');
xlabel('Time (s)'); ylabel('Residual Magnitude');
subplot(3,2,4); title('Horizon Sensor'); hold on; grid on;
plot(t, vecnorm(z_pre_horizon)-vecnorm(z_post_horizon)); legend('z_{pre} - z_{post}');
xlabel('Time (s)'); ylabel('Residual Magnitude');
subplot(3,2,5); title('Gyroscope'); hold on; grid on;
plot(t, vecnorm(z_pre_gyro)); plot(t,vecnorm(z_post_gyro)); legend('Pre-Fit', 'Post-Fit');
xlabel('Time (s)'); ylabel('Residual Magnitude');
subplot(3,2,6); title('Gyroscope'); hold on; grid on;
plot(t, vecnorm(z_pre_gyro)-vecnorm(z_post_gyro)); legend('z_{pre} - z_{post}');
xlabel('Time (s)'); ylabel('Residual Magnitude');

%% Run Kalman Filter with fine attitude sensing suite
% Initial conditions
variances = [450848.5655, 120673.3098, 80446.0581, 82931.6508, 0.3352, 1.91663, 0.06357]*10^(-6);
P0 = diag(variances).^2; Q = (1/100) * P0;
x0 = [q0'; w0_principal];

x_est = zeros(7,length(t));
P_est = zeros(7,7,length(t));
z_pre_gyro = zeros(3,length(t));
z_post_gyro = zeros(3,length(t));
z_pre_st = zeros(3,length(t));
z_post_st = zeros(3,length(t));

x_est(:,1) = x0;
P_est(:,:,1) = P0;
x = x0; P = P0;

% Run KF
disp("KF started")
idxIncrement = 1;
for idx = 2:length(t)
    
    % Predict step
    [x,P] = KF_TimeUpdate(x,P,Q,0,t(idx)-t(idx-1));

    % Model measurements
    wModel = gyroMeasurement(x(5:7)', false);               % modeled measurement
    [starModel, starRef, starIdxModel] = starTrackerMeasurement(x(1:4)', false);

    % Sequential update step - gyro
    wMeas = gyroMeasurement(w(idx,:), true); % actual measurement
    H = sensitivityMatrix("gyro", q(idx,:), 0);
    R = measurementCovarianceMatrix("gyro");
    [x, P, K] = KF_MeasurementUpdate(x,wMeas',wModel',P,H,R);

    % Sequential update step - star tracker
    [starMeas, starRef, starIdxMeas] = starTrackerMeasurement(q(idx,:),false);  % actual measurement
    [c, iMeas, iModel] = intersect(starIdxMeas, starIdxModel);  % which stars are in view in both the model and measurement?
    starModel = starModel(:,iModel); starMeas = starMeas(:,iMeas); starRef = starRef(iMeas,:);
    R = measurementCovarianceMatrix("star tracker");
    for k = 1:size(starModel,2)
        H = sensitivityMatrix("star tracker", q(idx,:), starRef(k,:));
        [x, P, K] = KF_MeasurementUpdate(x,starMeas(:,k),starModel(:,k),P,H,R);
    end

    x_est(:,idx) = x;
    P_est(:,:,idx) = P;

    % Pre- and post-fit residuals
    wModelPost = gyroMeasurement(x(5:7)', false);
end

x_est = x_est(:,1:length(t)); P_est = P_est(:,:,1:length(t));

Pq1 = reshape(P_est(1,1,:), [1,length(t)]); 
Pq2 = reshape(P_est(2,2,:), [1,length(t)]);
Pq3 = reshape(P_est(3,3,:), [1,length(t)]);
Pq4 = reshape(P_est(4,4,:), [1,length(t)]);
Pwx = reshape(P_est(5,5,:), [1,length(t)]);
Pwy = reshape(P_est(6,6,:), [1,length(t)]);
Pwz = reshape(P_est(7,7,:), [1,length(t)]);

disp("KF done (Fine Sensing)")

pltTitle = sprintf("Kalman Filter: Quaternion Attitude (Fine Sensing Suite)");
figure(); sgtitle(pltTitle);
subplot(4,3,1); hold on; title('q1');
plot(t, q(:,1)); plot(t, x_est(1,:));
grid on; hold off; legend('True', 'Estimated');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,2); hold on; title('q1 Error');
plot(t, x_est(1,:)'-q(:,1));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,3); hold on; title('q1 Covariance');
plot(t, sqrt(Pq1));
grid on; hold off; 
xlabel('Time (s)'); ylabel('\sigma');
subplot(4,3,4); hold on; title('q2');
plot(t, q(:,2)); plot(t, x_est(2,:));
grid on; hold off; legend('True', 'Estimated');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,5); hold on; title('q2 Error');
plot(t, x_est(2,:)'-q(:,2));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,6); hold on; title('q2 Covariance');
plot(t, sqrt(Pq2));
grid on; hold off; 
xlabel('Time (s)'); ylabel('\sigma');
subplot(4,3,7); hold on; title('q3');
plot(t, q(:,3)); plot(t, x_est(3,:));
grid on; hold off; legend('True', 'Estimated');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,8); hold on; title('q3 Error');
plot(t, x_est(3,:)'-q(:,3));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,9); hold on; title('q3 Covariance');
plot(t, sqrt(Pq3));
grid on; hold off; 
xlabel('Time (s)'); ylabel('\sigma');
subplot(4,3,10); hold on; title('q4');
plot(t, q(:,4)); plot(t, x_est(4,:));
grid on; hold off; legend('True', 'Estimated');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,11); hold on; title('q4 Error');
plot(t, x_est(4,:)'-q(:,4));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,12); hold on; title('q4 Covariance');
plot(t, sqrt(Pq4));
grid on; hold off; 
xlabel('Time (s)'); ylabel('\sigma');

pltTitle = sprintf("EKF Time Update\n Angular Velocity (Fine Sensing Suite)");
figure(); sgtitle(pltTitle);
subplot(2,2,1); hold on; title('True Angular Velocity');
plot(t, w(:,1)); plot(t, w(:,2)); plot(t, w(:,3));
grid on; hold off;
xlabel('Time (s)'); ylabel('Rad/s'); legend('\omega_x', '\omega_y', '\omega_z');
subplot(2,2,2); hold on; title('Estimated Angular Velocity');
plot(t, x_est(5,:)); plot(t, x_est(6,:)); plot(t, x_est(7,:));
grid on; hold off;
xlabel('Time (s)'); ylabel('Rad/s'); legend('\omega_x', '\omega_y', '\omega_z');
subplot(2,2,3); hold on; title('Angular Velocity Estimation Error');
plot(t, x_est(5,:)'-w(:,1)); plot(t, x_est(6,:)'-w(:,2)); plot(t, x_est(7,:)'-w(:,3));
grid on; hold off;
xlabel('Time (s)'); ylabel('Rad/s'); legend('\omega_x', '\omega_y', '\omega_z');
subplot(2,2,4); hold on; title('Covariances');
plot(t, sqrt(Pwx)); plot(t, sqrt(Pwy)); plot(t, sqrt(Pwz));
grid on; hold off;
xlabel('Time (s)'); ylabel('\sigma'); legend('\omega_x', '\omega_y', '\omega_z');

%% Useful Functions: Kalman Filter
% KF/EKF Measurement Update
function [xk, Pk, K] = KF_MeasurementUpdate(x,y,z,P,H,R)

    % Measurement update step
    K = P * H' * inv(H*P*H' + R);
    xk = x + K*(y-z);
    Pk = (eye(7) - K*H)*P*(eye(7) - K*H)' + K*R*K';

    % State involves a quaternion, which must be normalized
    q = xk(1:4);
    xk(1:4) = q ./ norm(q);
    
end

% KF/EKF Time Update
% x0 --> 7x1, state vector at prev. timestep
% P0 --> 7x7, state error cov at prev. timestep
% u0 --> nx1, control input
function [xk, Pk] = KF_TimeUpdate(x0,P0,Q,u0,dt)

    stateTransitionMatrix = @(dt) [1 0 0 0 0.5*dt -0.5*dt 0.5*dt;
                                   0 1 0 0 0.5*dt 0.5*dt -0.5*dt;
                                   0 0 1 0 -0.5*dt 0.5*dt 0.5*dt;
                                   0 0 0 1 -0.5*dt -0.5*dt -0.5*dt;
                                   0 0 0 0 1 0 0;
                                   0 0 0 0 0 1 0;
                                   0 0 0 0 0 0 1];
    
    B = zeros(7,1);         % Control input matrix - We currently have no controller! - 7xn

    % Implement time update step
    xk = stateTransitionMatrix(dt)*x0 + B*u0;
    Pk = stateTransitionMatrix(dt)*P0*stateTransitionMatrix(dt)' + Q;

    % State involves a quaternion, which must be normalized
    q = xk(1:4);
    xk(1:4) = q ./ norm(q);  

end

%% Problem 1/2: Define sensitivity matrix H and measurement error covariance R
% sensor --> sensor name as string (gyro, horizon, star tracker, sun sensor)
% q_true --> true attitude
% s_mci --> reference vector, if applicable
% H for sun sensor AND star tracker
% depends on current attitude and true MCI reference vector to star (or nadir vector in
% MCI for horizon)
% H for gyro is just identity
function [H] = sensitivityMatrix(sensor, q_true, s_mci)
    if sensor == "gyro"
        H = [zeros(3,4), eye(3,3)];
    else
        % Get attitude components
        q1 = q_true(1); q2 = q_true(2); q3 = q_true(3); q4 = q_true(4);
    
        % Partial derivatives
        dsx_dq1 = 2*q1*s_mci(1) + 2*q2*s_mci(2) + 2*q3*s_mci(3);
        dsx_dq2 = -2*q2*s_mci(1) + 2*q1*s_mci(2) - 2*q4*s_mci(3);
        dsx_dq3 = -2*q3*s_mci(1) + 2*q4*s_mci(2) + 2*q1*s_mci(3);
        dsx_dq4 = 2*q4*s_mci(1) + 2*q3*s_mci(2) - 2*q2*s_mci(3);
    
        dsy_dq1 = 2*q2*s_mci(1) - 2*q1*s_mci(2) + 2*q4*s_mci(3);
        dsy_dq2 = 2*q1*s_mci(1) + 2*q2*s_mci(2) + 2*q3*s_mci(3);
        dsy_dq3 = -2*q4*s_mci(1) - 2*q3*s_mci(2) + 2*q2*s_mci(3);
        dsy_dq4 = -2*q3*s_mci(1) + 2*q4*s_mci(2) + 2*q1*s_mci(3);
    
        dsz_dq1 = 2*q3*s_mci(1) - 2*q4*s_mci(2) - 2*q1*s_mci(3);
        dsz_dq2 = 2*q4*s_mci(1) + 2*q3*s_mci(2) - 2*q2*s_mci(3);
        dsz_dq3 = 2*q1*s_mci(1) + 2*q2*s_mci(2) + 2*q3*s_mci(3);
        dsz_dq4 = 2*q2*s_mci(1) - 2*q1*s_mci(2) + 2*q4*s_mci(3);
    
        H = [dsx_dq1, dsx_dq2, dsx_dq3, dsx_dq4, 0, 0, 0;
                 dsy_dq1, dsy_dq2, dsy_dq3, dsy_dq4, 0, 0, 0;
                 dsz_dq1, dsz_dq2, dsz_dq3, dsz_dq4, 0, 0, 0];
    end
end

function R = measurementCovarianceMatrix(sensor)
    if sensor == "gyro"
        load sensorProperties sigmaGyro_rads
        R = (sigmaGyro_rads)^2 * eye(3);
    elseif sensor == "star tracker"
        load sensorProperties sigmaST
        R = (sigmaST)^2 * eye(3);
    elseif sensor =="sun sensor"
        load sensorProperties sigmaSS_rad
        R = (sigmaSS_rad)^2 * eye(3);
    elseif sensor == "horizon"
        load sensorProperties sigmaHorizon
        R = (sigmaHorizon)^2 * eye(3);
    end
end

%% Useful Functions: sensor modeling
% Define sensor offsets, noise, other properties
function setSensorProperties()

    % Gyroscope - see Wertz p. 200
    load sensorProperties
    clear sigmaGyro
    A_gyro = 1;         % m^2
    sigmaGyro_rads = (0.003 * (2*pi/180) / 3600);                 % rad/s
    sigmaGyro_dt = sigmaGyro_rads * 4 * A_gyro / (3e8 * 3e8);       % measured time
    biasGyro = sigmaGyro_dt * 0.001;
    gyroOffset = biasGyro * (2*rand([1 3]) - 1)

    % Horizon sensor - determines nadir direction
    sigmaHorizon = 0.1 * pi / 180;      % rad
    biasHorizon = sigmaHorizon * 0.001;       % rad
    offsetHorizon = biasHorizon * (2*rand([1 3]) - 1);
    
    % Star tracker - same as before, see Wertz p. 190
    sigmaST = 1.0908e-5;
    starTrackerBias = sigmaST * 0.001;
    starTrackerOffset = biasST * (2*rand([3 1]) - 1);

    % Sun sensor - use apparent size of Sun from Mars, then scale by current
    ssMaxCurrent = 0.1;     % A
    sigmaSS_rad = 0.35 * pi / 180;  % rad
    sigmaSS_current = ssMaxCurrent * (0.35 * pi / 180);      % A
    biasSS = sigmaSS_rad * 0.001;
    sunSensorOffset = biasSS * (2*rand([3 1]) - 1);
    
    save sensorProperties
end

% Horizon Sensor
function [nadirMeasPrincipal,nadirRefMCI] = horizonSensorMeasurement(q_true, r, noise)

    load sensorProperties
    noise = normrnd(0, sigmaHorizon.^2, [3 1]);
    nadirRefMCI = (-1 * r ./ norm(r))';
    dcm_true = quaternion2DCM(q_true(1), q_true(2), q_true(3), q_true(4));
    nadirRefPrincipal = dcm_true * nadirRefMCI;
    nadirMeasPrincipal = nadirRefPrincipal;
    if (noise == true)
        nadirMeasPrincipal = nadirMeasPrincipal + noise + offsetHorizon';
    end
    nadirMeasPrincipal = nadirMeasPrincipal ./ norm(nadirMeasPrincipal);

end

% Ring-Laser Gyroscope
function wMeas = gyroMeasurement(w_true, noise)

    load sensorProperties sigmaGyro_dt gyroOffset A_gyro

    % Calculate expected measurement
    c = 3e8;        % speed of light [m/s]
    dt_expected = 4*A_gyro.*w_true ./ (c*c);     % expected elapsed time

    % Add noise & bias
    if (noise == true)
        gyroNoise = normrnd(0, sigmaGyro_dt^2, [1 3]);
        dt_meas = dt_expected + gyroOffset + gyroNoise;
        wMeas = dt_meas*c*c/(4*A_gyro);
    else
        wMeas = w_true;
    end
end

% Cosine Law Sun Sensor
function [sunMeasPrincipal, sunRefMCI] = sunSensorMeasurement(q_true, t, noise)
    load inertiaTensors rotationBodyToPrincipal
    load sensorProperties sigmaSS_current sunSensorOffset
    
    % True attitude
    dcmMCI2Principal = quaternion2DCM(q_true(1), q_true(2), q_true(3), q_true(4));
    
    % True sun position
    sunRefMCI = calculateSunVectorMCI(t);
    sunRefPrincipal = dcmMCI2Principal * sunRefMCI;
    sunBody = rotationBodyToPrincipal' * sunRefPrincipal;

    maxCurrent = 0.1; % [mA], Wertz 157
    angles = zeros(size(sunBody));

    % Measure along each axis
    for axis = 1:3
        for direction = [-1,1]
            sensorNormalVector = zeros(1,3);
            sensorNormalVector(axis) = direction;
            current = maxCurrent * dot(sensorNormalVector, sunBody);
            if (current > 0)    % sun is in sensor FOV

                % add noise and bias, then limit for saturation
                if (noise == true)
                    current = current + normrnd(0, sigmaSS_current^2, [1 1]) + sunSensorOffset(axis);
                    current = min(current, maxCurrent); current = max(current, -1*maxCurrent);
                end

                % calculate angle
                angle = acosd(current ./ maxCurrent);
                if (direction == -1)
                    angle = 180 - angle;
                end
                angles(axis) = angle;
            end
        end
    end

    % Rotate back to principal axis frame
    sunMeasuredBody = cosd(angles);
    sunMeasPrincipal = rotationBodyToPrincipal * sunMeasuredBody;
end

% Star Tracker
function [starMeasPrincipal, starRefMCI, starIndexes] = starTrackerMeasurement(q_true, noise)
    load sensorProperties sigmaST starTrackerOffset starCatalog
    dcm_true = quaternion2DCM(q_true(1), q_true(2), q_true(3), q_true(4));  
    starTrackerNoise = normrnd(0, sigmaST^2, [3 1]); 
    starTrackerAlignment = inv(dcm_true) * [1,0,0]';     % angle of star tracker in MCI frame, assuming it points along s/c x-axis
    starsInView = (starCatalog*starTrackerAlignment > 0);
    starIndexes = find(starsInView);
    starRefMCI = starCatalog(starIndexes,:);
    starRefPrincipal = (dcm_true*starRefMCI');
    if (noise == true)
        starMeasPrincipal = starRefPrincipal + starTrackerNoise + starTrackerOffset;
    else
        starMeasPrincipal = starRefPrincipal;
    end
    starMeasPrincipal = starMeasPrincipal ./ vecnorm(starMeasPrincipal')';  % normalize
end

%% Useful Functions: Integration
% Integrate orbit and attitude with ALL environmental torques
function [statedot] = trueOrbitAndAttitude(t, state)
    % Constants
    load inertiaTensors Ixx Iyy Izz
    load orbitGeometry mu_mars startMJD secondsPerEarthDay 

    MJD = startMJD + (t/secondsPerEarthDay);
  

    % Normalize quaternion
    state(7:10) = state(7:10) ./ norm(state(7:10));

    % State is [rx ry rz vx vy vz q1 q2 q3 q4 wx wy wz]
    r = state (1:3);
    v = state (4:6);
    q1 = state(7); q2 = state(8); q3 = state(9); q4 = state(10);
    wx = state(11); wy = state(12); wz = state(13);
    statedot = zeros (13,1);

    % Acceleration due to 1/ r ^2 gravity
    normr = norm(r);
    rhat = r/normr;
    acc = - (mu_mars/(normr * normr)) * rhat ;
    statedot (1:3) = v;
    statedot (4:6) = acc;

    % Calculate rate of change of attitude, expressed as quaternion
    statedot(7:10) = quaternionEOM([q1 q2 q3 q4], wx, wy, wz);

    % Calculate angular acceleration with torque effects
    gg = gravityGradientTorque(r, v, [q1 q2 q3 q4]);
    drag = aerodynamicTorque(r, v, [q1 q2 q3 q4]);
    mag = magneticTorque(r, [q1 q2 q3 q4], MJD);
    srp = srpTorque(t, r, [q1 q2 q3 q4]);
    torque = gg + drag + mag + srp;
    statedot(11) = (torque(1)-(Izz - Iyy)*(wy*wz))/Ixx;
    statedot(12) = (torque(2)-(Ixx - Izz)*(wz*wx))/Iyy;
    statedot(13) = (torque(3)-(Iyy - Ixx)*(wx*wy))/Izz;
end

%% Useful Functions: Disturbance Torques

% Magnetic Field Torque
% R_sat --> S/C position [MCI frame]
% q --> S/C attitude [MCI frame]
% m_body --> dipole direction [MCI frame]
% B0_body --> planet magnetic field strength
function torque = magneticTorque(R_sat, q, MJD)

    % Constants (note that we model Mars' magnetic field like Earth's, but
    % with reduced magnitude)
    load spacecraftProperties m_sc  % s/c magnetic field
    load orbitGeometry R_mars marsRotationRate B0 % [km] [deg/day] [T]

    inertial2Principal = quaternion2DCM(q(1),q(2),q(3),q(4));       % DCM: Principal Axis rel. to inertial
    referenceMJD = 44237;           % reference MJD (Dec 31, 1979)
    greenwichRA = 98.8279;          % RA of Greenwich meridian on Dec 31, 1979
    dipoleLongitude = 109.3759;     % [deg]c
    dipoleCoelevation = 168.6;      % [deg]

    % Estimate planet magnetic field in INERTIAL frame
    R = norm(R_sat);        % s/c position vector
    Rhat_sat = R_sat/R;
    dipoleRA = mod(greenwichRA + (MJD - referenceMJD)*marsRotationRate + dipoleLongitude,360);
    m_body = [sind(dipoleCoelevation)*cosd(dipoleRA); sind(dipoleCoelevation)*sind(dipoleRA); cosd(dipoleCoelevation)];    % Wertz, H-23
    B_MCI = -B0 * (R_mars/R)^3 * (3*dot(m_body, Rhat_sat).*Rhat_sat - m_body);

    % Convert planet magnetic field to principal axis frame
    B_principal = inertial2Principal * B_MCI;

    % Calculate torque
    torque = cross(m_sc, B_principal);
end

% SRP Torque
% For now, we assume Mars' orbital speed is negligible
% r --> s/c position in MCI frame
% q --> s/c attitude in MCI frame
function [torque] = srpTorque(t, r, q)
    
    % Constants
    torque = 0;
    load orbitGeometry R_mars
    load spacecraftProperties exposedSurfaceAreas exposedUnitVectors exposedCentroids Cs Cd

    P = 586.2/(3e8);                            % Solar pressure at Mars
    sunVectorMCI = calculateSunVectorMCI(t);

    % check if in eclipse
    
    if (eclipseCheck(r,t) == false)
        inertial2Principal = quaternion2DCM(q(1),q(2),q(3),q(4));
        sunVectorPrincipal = inertial2Principal * sunVectorMCI;

        % Calculate contribution from each surface
        for j = 1:length(exposedSurfaceAreas)
    
            % Check that surface is unshaded
            if dot(sunVectorPrincipal, exposedUnitVectors(j, 1:3)) > 0
                theta = acos(dot(sunVectorPrincipal,exposedUnitVectors(j, 1:3)));   % Angle between sun vector and unit normal
                
                % Radiation force on surface in principal axis frame
                df = (-1*P) * ((1-Cs)*sunVectorPrincipal + 2*(Cs*cos(theta) + (1/3)*Cd)*exposedUnitVectors(j, 1:3)')*cos(theta)*exposedSurfaceAreas(j);
                dM = cross(exposedCentroids(j,1:3), df);
                torque = torque + dM;
            end
        end
    end
end

% Atmospheric Drag Torque
% r_mci, v_mci --> s/c position and velocity in MCI
% q --> current attitude
function torque = aerodynamicTorque(r_mci, v_mci, q)

    % Constants
    load spacecraftProperties exposedSurfaceAreas exposedUnitVectors exposedCentroids
    load orbitGeometry R_mars
    inertial2Principal = quaternion2DCM(q(1),q(2),q(3),q(4));      % DCM: Principal Axis rel. to inertial
    v_principal = inertial2Principal * v_mci;                      % S/C velocity expressed in principal axis frame
    altitude = 1000 * (norm(r_mci) - R_mars);                      % altitude [m]
    rho = (0.699*exp(-0.00009*altitude)) / (0.1921 * (249.7 - 0.00222*altitude));   % atmospheric density [kg/m^3]
    Cdrag = 1.28;                     % unitless drag coefficient - approximate as flat plate normal to flow
    vnorm = norm(v_principal); vhat = v_principal / vnorm;
    % Total torque
    torque = 0;
    for j = 1:length(exposedSurfaceAreas)
        df = (-1/2) * Cdrag * rho * vnorm*vnorm * dot(vhat, exposedUnitVectors(j,1:3)) * vhat * exposedSurfaceAreas(j);
        dM = cross(exposedCentroids(j,1:3), df);
        torque = torque + dM;
    end
end

% Gravity gradient torque in princpal axes (lec7, slide 9) 
function torque = gravityGradientTorque(r_mci, v_mci, q)

    % constants
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    load inertiaTensors Ixx Iyy Izz
    load orbitGeometry mu_mars
    rvec = [r_mci(1), r_mci(2), r_mci(3)]; vvec = [v_mci(1) v_mci(2) v_mci(3)];
    R = norm(rvec);

    % Calculate graivty gradient torque in principal axes
    inertial2Principal = quaternion2DCM(q1,q2,q3,q4);       % DCM: Principal Axis rel. to inertial
    inertial2RTNDCM = inertialToRTN(rvec, vvec);            % DCM: RTN rel. to inertial
    RTN2Principal = inertial2Principal * inertial2RTNDCM';  % DCM: principal axis rel. to RTN
    r_rtn = [R,0,0];                                        % position vector in RTN frame
    r_xyz = RTN2Principal * r_rtn';                         % position vector in XYZ frame
    cx = r_xyz(1)/R; cy = r_xyz(2)/R; cz = r_xyz(3)/R;
    torque = (3*mu_mars/(R^3))*[(Izz-Iyy)*cy*cz; (Ixx-Izz)*cz*cx; (Iyy-Ixx)*cx*cy];
end

%% Useful Functions: Orbit Geometry

% Calculate sun position
function r_sun_mci = calculateSunVectorMCI(t)
    load sunOrbitElements
    load orbitGeometry marsYear
    e_mars = 0.0934;                    % Mars eccentricity w.r.t. Sun
    apparentMotion = deg2rad(360 / marsYear);    % apparent motion of sun, in rad/s
    M_current = sunInitialMeanAnomaly + apparentMotion*t;   % current mean anomaly of sun
    E_current = calculateEccentricAnomaly(M_current,e_mars,1e-8);   % current eccentric anomaly of sun
    currentTrueAnomaly = 2*atand(sqrt((1+e_mars)/(1-e_mars))*tan(E_current/2));
    [sunVectorMCI, sunVelocityMCI] = oe2mci(as,es,is,raans,argps,currentTrueAnomaly);
    r_sun_mci = sunVectorMCI ./ norm(sunVectorMCI);
end

% Check if in eclipse
% r --> s/c position in MCI
% t --> elapsed simulation time
function eclipse = eclipseCheck(r, t)
      load orbitGeometry R_mars
      sunVectorMCI = calculateSunVectorMCI(t);
      r_perp = r - (sunVectorMCI * dot(r, sunVectorMCI));
      if ((norm(r_perp) > R_mars || dot(r, sunVectorMCI) > 0))
          eclipse = false;
      else
          eclipse = true;
      end
end

% Calculate eccentric anomaly
% Params: M --> mean anomaly (rad), e --> eccentricity, tolerance (rad)
% Returns: E --> eccentric anomaly (rad)
function E = calculateEccentricAnomaly(M,e,tolerance)
    E = 0;
    delta = Inf;
    while (abs(delta) >= tolerance)
        delta = (E - e*sin(E) - M) / (1 - e*cos(E));
        E = E - delta;
    end
end

% Converts orbit elements to Mars-centered inertial frame
function [r_mci, v_mci] = oe2mci(a,e,i,raan,argp,trueAnomaly)

    mu_mars = 4.28284e4;    % gravitational parameter (km^3/s^2)
    E = 2*atand(sqrt((1-e)/(1+e))*tand(trueAnomaly/2));   % eccentric anomaly (deg)
    n = sqrt(mu_mars/(a^3));     % mean motion (1/s)
    % perifocal position (km)
    r_pqw = [a*(cosd(E) - e), a*sqrt(1-e^2)*sind(E), 0];     
    % perifocal velocity (km/s)
    v_pqw = (a*n)/(1-e*cosd(E)) .* [-1*sind(E), sqrt(1-e^2)*cosd(E), 0];
    % Rotations
    R1 = [cosd(-raan) sind(-raan) 0;
          -sind(-raan) cosd(-raan) 0;
          0 0 1];
    R2 = [1 0 0;
          0 cosd(-i) sind(-i);
          0 -sind(-i) cosd(-i)];
    R3 = [cosd(-argp) sind(-argp) 0;
          -sind(-argp) cosd(-argp) 0;
          0 0 1];
    R = R1*R2*R3;
    r_mci = R*r_pqw';    % MCI position(km)
    v_mci = R*v_pqw';    % MCI velocity (km)
end

% Calculate rotation matrix from inertial to RTN
function [dcm] = inertialToRTN(rvec,vvec)
    
    Rhat = rvec ./ norm(rvec);
    Nhat = cross(rvec,vvec) ./ norm(cross(rvec,vvec));
    That = cross(Nhat, Rhat);
    % Initial conditions
    dcm = [Rhat; That; Nhat];
end

%% Useful Functions: Attitude Representations
% Convert DCM to Euler (312)
function [theta, phi, psi] = dcm2Euler312(A)
    theta = asin(A(3,2));
    phi = atan2(A(1,2), A(2,2));
    psi = atan2(A(3,1),A(3,3));
end

% Convert DCM to Euler (313)
function [theta, phi, psi] = dcm2Euler313(A)
    theta = acos(A(3,3));
    phi = -atan2(A(3,1),A(3,2));
    psi = atan2(A(1,3),A(2,3));
end

% Convert Quaternion to DCM - Q4 is scalar component
function [dcm] = quaternion2DCM(Q1,Q2,Q3,Q4) 
    magnitude = Q1^2 + Q2^2 + Q3^2 + Q4^2;
    dcm = [Q4^2+Q1^2-Q2^2-Q3^2,   2*(Q1*Q2+Q3*Q4),        2*(Q1*Q3-Q2*Q4);
         2*(Q1*Q2-Q3*Q4),       Q4^2-Q1^2+Q2^2-Q3^2,    2*(Q2*Q3+Q1*Q4);
         2*(Q1*Q3+Q2*Q4),       2*(Q2*Q3-Q1*Q4),        Q4^2-Q1^2-Q2^2+Q3^2];
    dcm = dcm ./ magnitude;
end

% Convert DCM to quaternion
function [q1,q2,q3,q4] = dcm2Quaternion(A)
    q4 = (1/2) * (1 + A(1,1) + A(2,2) + A(3,3))^(1/2);
    q1 = (A(2,3) - A(3,2)) / (4*q4);
    q2 = (A(3,1) - A(1,3)) / (4*q4);
    q3 = (A(1,2) - A(2,1)) / (4*q4);
end

% Kinematic EOM for quaternion
function [qdot] = quaternionEOM(q, wx, wy, wz)
    omega = [0 wz -wy wx; -wz 0 wx wy; wy -wx 0 wz; -wx -wy -wz 0];
    q = q ./ norm(q);            % Normalize
    qdot = (1/2) * omega * q';   % Calculate derivative
end