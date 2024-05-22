%% Problem 1: Set Up Sensor Properties (basic settings)
clear; clc;

% Gyroscope
biasGyro = 0.01*pi/180;                      % lecture 1, spinning gyro
sigmaGyro = 0.003 * (pi/180) * (1/3600);     % deg/hr to rad/s (Wertz p.200) 
gyroOffset = biasGyro * (2*rand([1 3]) - 1); % random offset from -bias to + bias

% Sun sensor
biasSS = 5 * pi/180;                     % lecture 1, digital SS
sigmaSS = 0.5 * pi/180;                  % deg to rad
sunSensorOffset = biasSS * (2*rand([3 1]) - 1); % random offset from -bias to + bias

% Star Tracker
biasST = 0.01 * pi/180;                             % lecture 1
sigmaST = 2.25 * pi / (180*3600);                   % arcsec to rad/s (Wertz 190)
starTrackerOffset = biasST * (2*rand([3 1]) - 1);   % random offset from -bias to + bias
starCatalog = [1,0,0; 0,1,0; 0,0,1; -sqrt(2)/2, 0, sqrt(2)/2; 0, sqrt(3)/2, 1/2];       % Known star vectors in MCI frame

% save sensorProperties

%% Problem 4: Set Sensor Properties for Detailed Hardware Models
clear; clc;

% Ring Laser Gyroscope
biasGyro = 1e-25;
sigmaGyro = 1e-27;  
gyroOffset = biasGyro * (2*rand([1 3]) - 1); % random offset from -bias to + bias

% Sun sensor cosine law
biasSS = 1e-5;
sigmaSS = 5e-4;
sunSensorOffset = biasSS * (2*rand([3 1]) - 1); % random offset from -bias to + bias

% Star Tracker
n_stars = 50;
biasST = 0.01 * pi/180;                             % lecture 1
sigmaST = 2.25 * pi / (180*3600);                   % arcsec to rad/s (Wertz 190)
starTrackerOffset = biasST * (2*rand([3 1]) - 1);   % random offset from -bias to + bias
starCatalog = 2 * rand([n_stars 3]) - 1;            % generate random stars
starCatalog = starCatalog ./ vecnorm(starCatalog')';            % normalize

% save sensorProperties

%% Problem 2: Attitude Estimation Error - Run Simulation
load simulationSettings
load orbitGeometry
load spacecraftProperties

% Set up initial position
[r0, v0] = oe2mci(a,e,i,raan,argp,trueAnomaly0);        % initial position w.r.t. inertial frame

% Set up initial angular velocity
w0_body = [marsMeanMotion, 0, 0];                       % Initial angular velocity should be equal to Mars' mean motion
w0_principal = rotationBodyToPrincipal * w0_body';      % Convert to principal axis frame (most calculations will be in this frame for simplicity)

% Set up initial attitude
zb_target = calculateSunVectorMCI(0);
xb_target = cross(r0, zb_target) ./ norm(cross(r0, zb_target));
yb_target = cross(zb_target, xb_target);
dcm_initial_body = [xb_target'; yb_target'; zb_target'];
dcm_initial_principal = rotationBodyToPrincipal * dcm_initial_body;
[q10, q20, q30, q40] = dcm2Quaternion(dcm_initial_principal); % initial attitude w.r.t. inertial frame
q0 = [q10, q20, q30, q40];
y0 = [r0', v0', q0, w0_principal'];               % initial state vector
disp("Initial conditions ready!")

% Run "ground truth" simulation
[t, y_out] = ode113(@trueOrbitAndAttitude, [0:tstep:t_period]' , y0, options);
pos = y_out(:,1:3); vel = y_out(:,4:6); q = y_out(:,7:10); w = y_out(:,11:13);
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
errorTheta = zeros(size(theta)); errorPhi = zeros(size(theta)); errorPsi = zeros(size(theta));
disp("Simulation done!");

%% Problem 2: Attitude Estimation Error - Plot Results
% Measure attitude throughout simulation
q_meas = zeros(length(t),4);
w_meas = zeros(length(t),3);
for idx = 1:length(t)
    [qmt, wmt] = measureAttitude(q(idx,:), w(idx,:), pos(idx,:), t(idx));
    q_meas(idx,:) = qmt;
    w_meas(idx,:) = wmt;
end
disp("Measurements modeled");

% Extract Euler angles and attitude estimation errors
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
theta_meas = zeros(size(t)); phi_meas = zeros(size(t)); psi_meas = zeros(size(t));
for k = 1:size(t)
    dcm_true = quaternion2DCM(q(k,1),q(k,2),q(k,3),q(k,4));
    dcm_meas = quaternion2DCM(q_meas(k,1),q_meas(k,2),q_meas(k,3),q_meas(k,4));
    [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm_true);
    [theta_meas(k), phi_meas(k), psi_meas(k)] = dcm2Euler312(dcm_meas);
end
disp("Euler Angles done!")

pltTitle = sprintf("MAVEN Orbit\n Measured vs. True Attitude About Principal Axes");
figure(); sgtitle(pltTitle);
subplot(2,3,1); hold on; title('Angular Velocity (Measured)');
plot(t, w_meas(:,1)); plot(t, w_meas(:,2)); plot(t, w_meas(:,3));
grid on; hold off; ylim([-5e-4, 5e-4])
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); legend('X', 'Y', 'Z');

subplot(2,3,2); hold on; title('Angular Velocity (Actual)');
plot(t, w(:,1)); plot(t, w(:,2)); plot(t, w(:,3));
grid on; hold off; ylim([-5e-4, 5e-4]);
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); legend('X', 'Y', 'Z');

subplot(2,3,3); hold on; title('Angular Velocity Error');
plot(t, w_meas(:,1)-w(:,1)); plot(t,  w_meas(:,2)-w(:,2)); plot(t,  w_meas(:,3)-w(:,3));
grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); legend('X', 'Y', 'Z');

subplot(2,3,4); hold on; title('Euler Angles (Measured)'); hold on;
plot(t, phi_meas); plot(t, theta_meas); plot(t, psi_meas); grid on;
xlabel('Time (s)'); ylabel('Radians'); legend('\phi', '\theta', '\psi');

subplot(2,3,5); hold on; title('Euler Angles (Actual)'); hold on;
plot(t, phi); plot(t, theta); plot(t, psi); grid on;
xlabel('Time (s)'); ylabel('Radians'); legend('\phi', '\theta', '\psi');

subplot(2,3,6); hold on; title('Euler Angle Error'); hold on;
plot(t, phi_meas-phi); plot(t, theta_meas-theta); plot(t, psi_meas-psi); grid on;
xlabel('Time (s)'); ylabel('Radians'); legend('\phi', '\theta', '\psi');

%% Problem 3: Error DCM
errorQ = q_meas(:,1:4) - q(:,1:4);
errorDCM = zeros(3,3,length(t));
for idx = 1:length(t)
    errorDCM(:,:,idx) = quaternion2DCM(errorQ(idx,1),errorQ(idx,2),errorQ(idx,3),errorQ(idx,4));
end

xError = reshape(errorDCM(2,3,:), [1 length(t)]); 
yError = reshape(errorDCM(3,1,:), [1 length(t)]); 
zError = reshape(errorDCM(1,2,:), [1 length(t)]); 
figure(); hold on; grid on;
plot(t,xError); plot(t,yError); plot(t,zError);
legend('X', 'Y', 'Z');
xlabel('Time (s)'); ylabel('Error (Radians)');
pltTitle = sprintf('Attitude Measurement Error About Principal Axes\n Small Euler Angle Approximation');
title(pltTitle); hold off;

%% Problem 5: Test EKF
% Define initial states
variances = abs(mean([q-q_meas, w-w_meas]));
P0 = diag(variances);
Q = (1/100)*P0;
x0 = [q0'; w0_principal];

% Run EKF
x_est = zeros(7,length(t));
P_est = zeros(7,7,length(t));
x_est(:,1) = x0;
P_est(:,:,1) = P0;
x = x0; P = P0;
for idx = 2:length(t)
    [x,P] = EKF_TimeUpdate(x,P0,Q,0,t(idx)-t(idx-1));
    x_est(:,idx) = x; P_est(:,:,idx) = P;
end
Pq1 = reshape(P_est(1,1,:), [1,length(t)]); 
Pq2 = reshape(P_est(2,2,:), [1,length(t)]);
Pq3 = reshape(P_est(3,3,:), [1,length(t)]);
Pq4 = reshape(P_est(4,4,:), [1,length(t)]);
Pwx = reshape(P_est(5,5,:), [1,length(t)]);
Pwy = reshape(P_est(6,6,:), [1,length(t)]);
Pwz = reshape(P_est(7,7,:), [1,length(t)]);

disp("EKF done")

%% Problem 5: Compare to rigorous propagation
[t, y_out] = ode113(@torqueFreeOrbitAndAttitude, [0:tstep:t_period]' , y0, options);
tf_pos = y_out(:,1:3); tf_vel = y_out(:,4:6); tf_q = y_out(:,7:10); tf_w = y_out(:,11:13);

pltTitle = sprintf("State Propagation vs. Torque-Free Numerical Simulation");
figure(); sgtitle(pltTitle);

subplot(1,3,1); hold on; title('Propagated Attitude');
plot(t, x_est(1,:)); plot(t, x_est(2,:)); plot(t, x_est(3,:)); plot(t, x_est(4,:));
grid on; hold off;
xlabel('Time (s)'); ylabel('Quaternion'); legend('q1', 'q2', 'q3', 'q4');

subplot(1,3,2); hold on; title('Integrated Attitude');
plot(t, tf_q(:,1)); plot(t, tf_q(:,2)); plot(t, tf_q(:,3)); plot(t, tf_q(:,4));
grid on; hold off;
xlabel('Time (s)'); ylabel('Attitude'); legend('q1', 'q2', 'q3', 'q4');

subplot(1,3,3); hold on; title('Error');
plot(t, x_est(1,:)'-tf_q(:,1)); plot(t, x_est(2,:)'-tf_q(:,2)); plot(t, x_est(3,:)'-tf_q(:,3)); plot(t, x_est(4,:)'-tf_q(:,4));
grid on; ylim([-1e-3, 1e-3]); hold off;
xlabel('Time (s)'); ylabel('Attitude'); legend('q1', 'q2', 'q3', 'q4');

%% Problem 5: Plot errors
pltTitle = sprintf("EKF Time Update\n Quaternion Attitude");
figure(); sgtitle(pltTitle);

subplot(4,3,1); hold on; title('q1');
plot(t, q(:,1)); plot(t, x_est(1,:));
grid on; hold off; legend('True', 'Measured');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,2); hold on; title('q1 Error');
plot(t, x_est(1,:)'-q(:,1));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,3); hold on; title('q1 Covariance');
plot(t, Pq1);
grid on; hold off; 
xlabel('Time (s)'); ylabel('Covariance');

subplot(4,3,4); hold on; title('q2');
plot(t, q(:,2)); plot(t, x_est(2,:));
grid on; hold off; legend('True', 'Measured');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,5); hold on; title('q2 Error');
plot(t, x_est(2,:)'-q(:,2));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,6); hold on; title('q2 Covariance');
plot(t, Pq2);
grid on; hold off; 
xlabel('Time (s)'); ylabel('Covariance');

subplot(4,3,7); hold on; title('q3');
plot(t, q(:,3)); plot(t, x_est(3,:));
grid on; hold off; legend('True', 'Measured');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,8); hold on; title('q3 Error');
plot(t, x_est(3,:)'-q(:,3));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,9); hold on; title('q3 Covariance');
plot(t, Pq3);
grid on; hold off; 
xlabel('Time (s)'); ylabel('Covariance');


subplot(4,3,10); hold on; title('q4');
plot(t, q(:,4)); plot(t, x_est(4,:));
grid on; hold off; legend('True', 'Measured');
xlabel('Time (s)'); ylabel('Attitude');
subplot(4,3,11); hold on; title('q4 Error');
plot(t, x_est(4,:)'-q(:,4));
grid on; hold off; 
xlabel('Time (s)'); ylabel('Error');
subplot(4,3,12); hold on; title('q4 Covariance');
plot(t, Pq4);
grid on; hold off; 
xlabel('Time (s)'); ylabel('Covariance');

pltTitle = sprintf("EKF Time Update\n Angular Velocity");
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
plot(t, Pwx); plot(t, Pwy); plot(t, Pwz);
grid on; hold off;
xlabel('Time (s)'); ylabel('\sigma ^2'); legend('\omega_x', '\omega_y', '\omega_z');

%% Problem 4: sensor modeling
% Angular velocity vector obtained from gyroscope about PRINCIPAL AXES
% Modeled as a ring-laser gyroscope
function [w_meas] = gyroMeasurement(q_true, w_true, t)

    % Calculate expected measurement
    c = 3e8;        % speed of light [m/s]
    A = 1;          % RLG area [m^2]
    dt_expected = 4*A.*w_true ./ (c*c);     % expected elapsed time

    % Add noise & bias
    load sensorProperties sigmaGyro gyroOffset
    gyroNoise = normrnd(0, sigmaGyro, [1 3]);
    dt_meas = dt_expected + gyroOffset + gyroNoise;
    w_meas = dt_meas*c*c/(4*A);
end

% Sun position vector in PRINCIPAL AXIS FRAME obtained from sun sensors
function [sunMeasuredPrincipal, sunReferencePrincipal] = sunSensorMeasurement(q_true, w_true, t)
    load inertiaTensors rotationBodyToPrincipal
    load sensorProperties sigmaSS sunSensorOffset
    
    % True attitude
    dcmMCI2Principal = quaternion2DCM(q_true(1), q_true(2), q_true(3), q_true(4));
    
    % True sun position
    sunMCI = calculateSunVectorMCI(t);
    sunReferencePrincipal = dcmMCI2Principal * sunMCI;
    sunBody = rotationBodyToPrincipal' * sunReferencePrincipal;

    maxCurrent = 0.1; % [mA], Wertz 157
    angles = zeros(size(sunBody));
    sunMeasuredBody = zeros(size(sunBody));

    % Measure along each axis
    for axis = 1:3
        for direction = [-1,1]
            sensorNormalVector = zeros(1,3);
            sensorNormalVector(axis) = direction;
            current = maxCurrent * dot(sensorNormalVector, sunBody);
            if (current > 0)    % sun is in sensor FOV

                % add noise and bias, then limit for saturation
                current = current + normrnd(0, sigmaSS, [1 1]) + sunSensorOffset(axis);
                current = min(current, maxCurrent); current = max(current, -1*maxCurrent);

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
    sunMeasuredPrincipal = rotationBodyToPrincipal * sunMeasuredBody;
end

% Return measured star vectors (principal axis frame) and corresponding reference
% star vectors (inertial frame)
function [starVectorMeasured, starVectorReference] = starTrackerMeasurement(q_true, w_true, t)
    load sensorProperties sigmaST starTrackerOffset starCatalog
    dcm_true = quaternion2DCM(q_true(1), q_true(2), q_true(3), q_true(4));  
    starTrackerNoise = normrnd(0, sigmaST, [3 1]); 
    starTrackerAlignment = inv(dcm_true) * [1,0,0]';     % angle of star tracker in MCI frame, assuming it points along s/c x-axis
    starVectorReference = starCatalog(starCatalog*starTrackerAlignment > 0,:)';
    starVectorMeasured = (dcm_true*starVectorReference) + starTrackerNoise + starTrackerOffset;
    starVectorMeasured = starVectorMeasured ./ vecnorm(starVectorMeasured')';  % normalize
end

%% Problem 5: KF/EKF Time Update
% x0 --> 7x1, state vector at prev. timestep
% P0 --> 7x7, state error cov at prev. timestep
% u0 --> nx1, control input
function [xk, Pk] = EKF_TimeUpdate(x0,P0,Q,u0,dt)

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
end

%% Useful Functions: attitude determination
% Measure attitude and apply q-method
function [q_meas, w_meas] = measureAttitude(q_true, w_true, r, t)

    % Set up sensors, true attitude, eclipse
    load sensorProperties
    dcm_true = quaternion2DCM(q_true(1), q_true(2), q_true(3), q_true(4));                  % DCM from inertial to principal frame
    eclipse = eclipseCheck(r,t);

    % Obtain sensor measurements
    w_meas = gyroMeasurement(q_true, w_true, t);
    [starVectorMeasured, starVectorReference] = starTrackerMeasurement(q_true, w_true, t);
    if (eclipse == false)
        [sunSensorMeasured, sunSensorReference] = sunSensorMeasurement(q_true, w_true, t);
    end

    % Reconstruct attitude using Q-method
    if (eclipse == false)
        w_star = 0.8 * ones(1,length(starVectorMeasured));
        w_sun = 0.2;
        w = [w_star, w_sun];     % We trust the star tracker more than the sun sensor
        m = [starVectorMeasured, sunSensorMeasured];
        v = [starVectorReference, sunSensorReference];
    else
        w_star = 0.8 * ones(1,length(starVectorMeasured));
        w = [w_star];     % We trust the star tracker more than the sun sensor
        m = [starVectorMeasured];
        v = [starVectorReference];
    end

    q_meas = qMethod(w, m, v);
end

% Q-Method (Statistical Attitude Determination)
% w --> sensor weights row vector (star tracker, sun sensor)
% m --> sensor measurements as column vectors in principal axis frame (star tracker, sun sensor)
% v --> reference vectors as column vectors in MCI (star tracker, sun sensor)
function q = qMethod(w, m, v)
    W = m .* sqrt(w);       % weighted measurements
    U = v .* sqrt(w);       % weighted reference
    B = W * U';
    S = B + B';
    Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
    sig = trace(B);
    K = [S-eye(length(S))*sig, Z;
        Z', sig];
    [VV, DD] = eig(K);
    [mx,mxidx] = max(diag(DD));
    q = VV(:,mxidx);        % Last component is scalar
end

%% Useful Functions: Integration
% Integrate orbit and attitude with ALL environmental torques
function [statedot] = trueOrbitAndAttitude(t, state)
    % Constants
    load inertiaTensors Ixx Iyy Izz
    load orbitGeometry mu_mars startMJD secondsPerEarthDay 

    MJD = startMJD + (t/secondsPerEarthDay);
  
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

% Integrate orbit and attitude with NO environmental torques
function [statedot] = torqueFreeOrbitAndAttitude(t, state)
    % Constants
    load inertiaTensors Ixx Iyy Izz
    load orbitGeometry mu_mars startMJD secondsPerEarthDay 

    MJD = startMJD + (t/secondsPerEarthDay);
  
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


%% Useful Functions: Attitude Reference Vectors

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


%% Useful Functions: Orbit Geometry
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