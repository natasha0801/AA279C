%%  Initial Conditions
load simulationSettings
load orbitGeometry
load spacecraftProperties
load rwProperties I_RW maxRWSpeed


% Set up initial condition
[r0, v0] = oe2mci(a,e,i,raan,argp,trueAnomaly0);   % initial position w.r.t. inertial frame
w0_des = rotationBodyToPrincipal * [marsMeanMotion, 0, 0]';
w0_principal = w0_des;

zb_target = calculateSunVectorMCI(0);
xb_target = cross(r0, zb_target) ./ norm(cross(r0, zb_target));
yb_target = cross(zb_target, xb_target);
dcm_des = rotationBodyToPrincipal * [xb_target'; yb_target'; zb_target'];
[q01, q02, q03, q04] = dcm2Quaternion(dcm_des);       
q0_des = [q01, q02, q03, q04]';
q0 = q0_des;
q0 = q0 ./ norm(q0);

% True state
tstep = 0.1;
t = 0:tstep:600;
     
q_true = zeros(4,length(t)); 
q_true(:,1) = q0;
w_true = zeros(3,length(t)); 
w_true(:,1) = w0_principal;
pos_true = zeros(3,length(t)); 
pos_true(:,1) = r0;
vel_true = zeros(3,length(t)); 
vel_true(:,1) = v0;

% Environmental torques
torqueGG = zeros(3,length(t));
torqueMag = zeros(3,length(t));
torqueDrag = zeros(3,length(t));
torqueSRP = zeros(3,length(t));
torqueTotal = zeros(3,length(t));

% Estimated state
variances = [9296.8361, 5622.6309, 21349.0534, 1504.4700, .000068567, .0000403, .0012959]*10^(-6);
P0 = diag(variances).^2; Q = (1/10) * P0;
x_est = zeros(7,length(t));
P_est = zeros(7,7,length(t));
x_est(:,1) = [q0; w0_principal];
P_est(:,:,1) = P0;

% Attitude control
q_des = zeros(4,length(t));
w_des = zeros(3,length(t));
rwSpeeds = zeros(3,length(t));
rwTorques = zeros(3,length(t));
thrust = zeros(4,length(t));
thrusterTorques = zeros(3,length(t));
controlTorquesRW = zeros(3,length(t));
controlTorquesThrust = zeros(3,length(t));
q_des(:,1) = q0_des; w_des(:,1) = w0_des;
kp = [0.9356, 1.0143, 1.0831];
kd = [69.4156, 75.2541, 80.3552];

% Other things (modes, etc.)
desat = false;
desatThreshold = 1e-3;
eclipse = zeros(1,length(t));

load inertiaTensors Ixx Iyy Izz
load orbitGeometry mu_mars startMJD secondsPerEarthDay
disp("Initial conditions ready!")

% Run simulation
for j = 1:length(t) 

    %%%%% GROUND TRUTH %%%%%
    q_true(:,j) = q_true(:,j) ./ norm(q_true(:,j));
    r = pos_true(:,j); v = vel_true(:,j);
    MJD = startMJD + (t(j)/secondsPerEarthDay);
    eclipse(j) = eclipseCheck(r,t(j));

    % Acceleration due to 1/ r ^2 gravity
    r = pos_true(:, j);
    normr = norm(r);
    rhat = r/normr;
    acc = - (mu_mars/(normr * normr)) * rhat ;
    pos_true(:, j+1) = pos_true(:, j) + vel_true(:, j)*tstep;
    vel_true(:, j+1) = vel_true(:, j) + acc*tstep;

    % Calculate rate of change of attitude, expressed as quaternion
    q1 = q_true(1,j); q2 = q_true(2,j); q3 = q_true(3,j); q4 = q_true(4,j);
    wx = w_true(1,j); wy = w_true(2,j); wz = w_true(3,j);
    qdot = quaternionEOM([q1 q2 q3 q4], wx, wy, wz);
    q_true(:,j+1) = q_true(:,j) + qdot*tstep;

    % Calculate  torque effects
    gg = gravityGradientTorque(r, v, [q1 q2 q3 q4]);
    drag = aerodynamicTorque(r, v, [q1 q2 q3 q4]);
    mag = magneticTorque(r, [q1 q2 q3 q4], MJD);
    srp = srpTorque(t(j), r, [q1 q2 q3 q4]);
    envTorque = gg + drag + mag + srp;
    torqueGG(:,j) = gg;
    torqueDrag(:,j) = drag;
    torqueMag(:,j) = mag;
    torqueSRP(:,j) = srp;
    torqueTotal(:,j) = envTorque;

    %%%% ATTITUDE ESTIMATION %%%%%
    
    % Predict step
    [x,P] = KF_TimeUpdate(x_est(:,j),P_est(:,:,j),Q, rwTorques(:,j)+thrusterTorques(:,j),tstep);

    % Update Step
    if (eclipse(j) == false)    
        % Sun sensor can be used in daylight
        [x, P] = coarseMeasurementUpdate(x, P, q_true(:,j), w_true(:,j), pos_true(:,j), t(j));
    else
        % Star tracker must be used in eclipse
        [x, P] = fineMeasurementUpdate(x, P, q_true(:,j), w_true(:,j), pos_true(:,j), t(j));
    end

    % Save
    x_est(:,j+1) = x;
    P_est(:,:,j+1) = P;

    %%%% CALCULATE TARGET ATTITUDE %%%%
    % In daylight: point toward sun
    if (eclipse(j) == false)
        zb_target = calculateSunVectorMCI(t(j));
        xb_target = cross(pos_true(:,j), zb_target) ./ norm(cross(pos_true(:,j), zb_target));
        yb_target = cross(zb_target, xb_target);
        dcm_des = rotationBodyToPrincipal * [xb_target'; yb_target'; zb_target'];
        [qt1, qt2, qt3, qt4] = dcm2Quaternion(dcm_des);   
        q_des(:,j) = [qt1, qt2, qt3, qt4]';
        q_des(:,j) = q_des(:,j) ./ norm(q_des(:,j));
        w_des(:,j) = rotationBodyToPrincipal * [marsMeanMotion, 0, 0]';
    % In darkness: inertial pointing
    else  
        q_des(:,j) = q_des(:,j-1);
        w_des(:,j) = [0, 0, 0]';
    end


    %%%% ACTUATION %%%%
    % Check if we're done with desaturation
    if norm(x_est(5:7,j+1)-w_des(:,j)) <= desatThreshold
        desat = false;
    end

    % Reaction wheel controller
    if desat == false && max(rwSpeeds(:,j)) < maxRWSpeed && min(rwSpeeds(:,j)) > -1*maxRWSpeed
    
        controlTorquesRW(:,j) = ControlLaw(q_des(:,j), x_est(1:4,j), w_des(:,j), x_est(5:7,j), kp, kd);
        [rwt, rws] = ReactionWheels(controlTorquesRW(:,j), rwSpeeds(:,j), tstep);
        rwTorques(:,j+1) = rwt;
        rwSpeeds(:,j+1) = rws;

    % Desaturation with thrusters
    else

        % (1) turn off RW and "dump" momentum into s/c
        load rwProperties A_RW
        desat = true;
        momentumDump = A_RW * (I_RW * rwSpeeds(:,j));       % momentum dumped into s/c
        rwSpeeds(:,j+1) = [0;0;0];
        rwTorques(1,j) = momentumDump(1) / (Ixx * tstep);
        rwTorques(2,j) = momentumDump(2) / (Iyy * tstep);
        rwTorques(3,j) = momentumDump(3) / (Izz * tstep);

        % (2) turn on thrusters until rotation rate is slowed down
        controlTorquesThrust(:,j) = ControlLaw(q_des(:,j), x_est(1:4,j), w_des(:,j), x_est(5:7,j), kp, kd);
        [torqueT, thrustT] = Thrusters(controlTorquesThrust(:,j));
        thrusterTorques(:,j+1) = torqueT;
        thrust(:,j+1) = thrustT;
    end

    wdotx = (envTorque(1) + rwTorques(1,j) + thrusterTorques(1,j) - (Izz - Iyy)*(wy*wz))/Ixx;
    wdoty = (envTorque(2) + rwTorques(2,j) + thrusterTorques(2,j) - (Ixx - Izz)*(wz*wx))/Iyy;
    wdotz = (envTorque(3) + rwTorques(3,j) + thrusterTorques(3,j) - (Iyy - Ixx)*(wx*wy))/Izz;
    
    w_true(1,j+1) = w_true(1,j) + wdotx*tstep;
    w_true(2,j+1) = w_true(2,j) + wdoty*tstep;
    w_true(3,j+1) = w_true(3,j) + wdotz*tstep;
end
 


%% Attitude Estimation Results
% For simulations < 1 orbit, there will be one eclipse period:
eclipseStart = min(t(eclipse==true));
eclipseEnd = max(t(eclipse==true));

% Plot attitude estimation (ground truth, estimated, error)
figure(); sgtitle('Quaternion Attitude Estimation');
q_true = q_true(:,1:length(t)); w_true = w_true(:,1:length(t));
q_est = x_est(1:4,1:length(t)); w_est = x_est(5:7,1:length(t));

for q = 1:4
    subplot(2,4,(2*q)-1); hold on; grid on;
    stitle = sprintf("q%i",q);
    title(stitle);
    plot(t, q_est(q,:), 'b'); plot(t,q_true(q,:), 'r--');
    legend('Estimated', 'True');
    xlabel('Time (s)'); ylabel('Quaternion Attitude');

    hold off;

    subplot(2,4,2*q); hold on; grid on;
    stitle = sprintf("q%i Estimation Error",q);
    title(stitle);
    plot(t, q_est(q,:) - q_true(q,:))
    xlabel('Time (s)'); ylabel('Attitude Estimation Error'); 
    hold off;
end

figure(); sgtitle('Angular Velocity Estimation');
for ax = 1:3
    subplot(3,2,(2*ax)-1); hold on; grid on;
    stitle = sprintf("Axis %i", ax);
    title(stitle);
    plot(t, w_est(ax,:), 'b'); plot(t,w_true(ax,:), 'r--');
    legend('Estimated', 'True');
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); hold off;

    subplot(3,2,2*ax); hold on; grid on;
    stitle = sprintf("Axis %i Estimation Error", ax);
    title(stitle);
    plot(t, w_est(ax,:) - w_true(ax,:));
    xlabel('Time (s)'); ylabel('Angular Velocity Error(rad/s)'); hold off;
end


%% Attitude Control Results
% Attitude control (ground truth, target, error)
figure(); sgtitle('Quaternion Attitude Control');
for q = 1:4
    subplot(2,4,(2*q)-1); hold on; grid on;
    stitle = sprintf("q%i",q);
    title(stitle);
    plot(t, q_true(q,:), 'b'); plot(t,q_des(q,:), 'r--');
    legend('True', 'Target');
    xlabel('Time (s)'); ylabel('Quaternion Attitude'); hold off;

    subplot(2,4,2*q); hold on; grid on;
    stitle = sprintf("q%i Control Error",q);
    title(stitle);
    plot(t, q_est(q,:) - q_des(q,:))
    xlabel('Time (s)'); ylabel('Attitude Control Error'); hold off;
end

figure(); sgtitle('Angular Velocity Control');
for q = 1:3
    subplot(3,2,(2*q)-1); hold on; grid on;
    stitle = sprintf("Axis %i",q);
    title(stitle);
    plot(t, w_true(q,:), 'b'); plot(t,w_des(q,:), 'r--');
    legend('True', 'Target');
    xlabel('Time (s)'); ylabel('Rad/s'); hold off;

    subplot(3,2,2*q); hold on; grid on;
    stitle = sprintf("Angular Velocity Control Error",q);
    title(stitle);
    plot(t, q_est(q,:) - q_des(q,:))
    xlabel('Time (s)'); ylabel('Rad/s'); hold off;
end

%% Control Action Results
rwTorques = rwTorques(:,1:length(t));
rwSpeeds = rwSpeeds(:,1:length(t));
thrust = thrust(:,1:length(t));
thrusterTorques = thrusterTorques(:,1:length(t));

% Control Torques from RW
figure(); sgtitle('Reaction Wheel Control');
for ax = 1:3
    subplot(1,4,ax); hold on; grid on;
    plot(t, rwTorques(ax,:)); plot(t, controlTorquesRW(ax,:));
    legend('Actual', 'Command');
    xlabel('Time (s)'); ylabel('Torque (Nm)');
    stitle = sprintf("Axis %i Torque", ax);
    title(stitle);
end

subplot(1,4,4); title('RW Speeds'); hold on;
plot(t, rwSpeeds(1,:)); plot(t,rwSpeeds(2,:)); plot(t, rwSpeeds(3,:));
xlabel('Time (s)'); ylabel('Speed (rad/s)');
legend('X', 'Y', 'Z');
grid on;

% Thruster Torques
figure(); sgtitle('Thruster Control');
for ax = 1:3
    subplot(1,4,ax); hold on; grid on;
    plot(t, thrusterTorques(ax,:)); plot(t, controlTorquesThrust(ax,:));
    legend('Actual', 'Command');
    xlabel('Time (s)'); ylabel('Torque (Nm)');
    stitle = sprintf("Axis %i Torque", ax);
    title(stitle);
end

subplot(1,4,4); title('Thrust'); hold on;
plot(t, thrust(1,:)); plot(t,thrust(2,:)); plot(t, thrust(3,:)); plot(t, thrust(4,:));
xlabel('Time (s)'); ylabel('Thrust (N)');
legend('T1', 'T2', 'T3', 'T4');
grid on;

%% Environmental Torque Results
% Plot environmental torques
figure(); sgtitle('Environmental Torques');
subplot(3,2,1);
plot(t, vecnorm(torqueTotal)); xlabel('Time (s)'); ylabel('Torque (N m)');
title('Total Environmental Torques'); grid on;
subplot(3,2,2); 
plot(t, vecnorm(torqueMag));
xlabel('Time (s)'); ylabel('Torque (N m)');
title('Magnetic Torque'); grid on;
subplot(3,2,3); 
plot(t, vecnorm(torqueGG));
xlabel('Time (s)'); ylabel('Torque (N m)');
title('Gravity Gradient Torque'); grid on;
subplot(3,2,4); 
plot(t, vecnorm(torqueDrag));
xlabel('Time (s)'); ylabel('Torque (N m)');
title('Aerodynamic Torque'); grid on;
subplot(3,2,5); 
plot(t, vecnorm(torqueSRP));
xlabel('Time (s)'); ylabel('Torque (N m)');
title('Solar Radiation Pressure Torque'); grid on;
subplot(3,2,6);
plot(t, vecnorm(pos_true(:,1:length(t)))); xlabel('Time (s)'); ylabel('Orbit Radius (m)');
title('Orbit Radius'); grid on;

%% Useful Functions: Actuation & Control
% rwTorque, torqueCommand --> resulting & desired torque about principal axes (3x1)
% rwSpeed, rwSpeed0 --> individual rw speed before and after command (3x1)
% dt --> time interval (s)
function [rwTorque, rwSpeed] = ReactionWheels(torqueCommand, rwSpeed0, dt)
    load rwProperties

    ldot_rw = pinv(A_RW) * torqueCommand;
    ldot_rw = max(min(ldot_rw, maxRWTorque), -1*maxRWTorque);

    wdot_rw = ldot_rw / I_RW;
    rwSpeed = rwSpeed0 + (wdot_rw.*dt.*(1+normrnd(0,0.01, [3 1])));
    rwSpeed = max(min(rwSpeed, maxRWSpeed), -1*maxRWSpeed);

    rwTorque = A_RW * (I_RW * (rwSpeed - rwSpeed0)/dt);              
end

% Mc --> command torque (3x1)
% Torque --> output torque (3x1)
% Thrust --> individual thrusters (4x1)
function [torque, thrust] = Thrusters(Mc)
    load thrusterGeometry
    
    maxThrust = 1; minThrust = 0;       % N
    thrust = pinv(A_thrusters) * Mc;
    thrust = thrust .* (1 + 0.001 * randn(4,1));    % uncertainty

    % Limit for saturation
    if max(thrust) > maxThrust
        thrust = thrust + maxThrust - max(thrust);
    end
    if min(thrust) < minThrust
        thrust = thrust + minThrust - min(thrust);
    end

    % Check both ends (maxima + minima)
    thrust = max(min(thrust,maxThrust), minThrust);

    torque = A_thrusters * thrust;

end

% Control law for large tracking angles
function [Mc] = ControlLaw(q_des, q_est, w_des, w_est, kp, kd)

    Mc = [0;0;0];

    % Error DCM
    A_est = quaternion2DCM(q_est(1),q_est(2),q_est(3),q_est(4));
    A_des = quaternion2DCM(q_des(1),q_des(2),q_des(3),q_des(4));
    A_err = A_est*A_des';

    % Control torque about each axis
    for i = 1:3
        j = mod(i,3) + 1;
        k = mod(j,3) + 1;
        Mc(i) = -kp(i) * 0.5 * (A_err(j,k) - A_err(k,j)) - kd(i)*(w_est(i));
    end

end

%% Useful Functions: Kalman Filter

function [x, P] = coarseMeasurementUpdate(x, P, q_true, w_true, r, t)
    % Model measurements using estimated state and no noise
    wModel = gyroMeasurement(x(5:7), false);               
    [nadirModel, nadirRef] = horizonSensorMeasurement(x(1:4), r, false);
    [sunModel, sunRef] = sunSensorMeasurement(x(1:4), t, false);

    % Sequential update step - gyro
    wMeas = gyroMeasurement(w_true, true); % actual measurement
    H = sensitivityMatrix("gyro", q_true, 0);
    R = measurementCovarianceMatrix("gyro");
    [x, P, K] = KF_MeasurementUpdate(x,wMeas,wModel,P,H,R);

    % Sequential update step - horizon sensor
    [nadirMeas, nadirRef] = horizonSensorMeasurement(q_true, r, true);  % actual measurement
    H = sensitivityMatrix("horizon", q_true, nadirRef);
    R = measurementCovarianceMatrix("horizon");
    [x, P, K] = KF_MeasurementUpdate(x,nadirMeas, nadirModel, P,H,R);

    % Sequential update step - sun sensor
    if (eclipseCheck(r, t) == false)
        [sunMeas, sunRef] = sunSensorMeasurement(q_true, t, true); % actual measurement
        H = sensitivityMatrix("sun sensor", q_true, sunRef);
        R = measurementCovarianceMatrix("sun sensor");
        [x, P, K] = KF_MeasurementUpdate(x,sunMeas,sunModel,P,H,R);
    end
end


function [x, P] = fineMeasurementUpdate(x, P, q_true, w_true, r, t)
    % Model measurements
    wModel = gyroMeasurement(x(5:7), false);               % modeled measurement
    [starModel, starRef, starIdxModel] = starTrackerMeasurement(x(1:4), false);

    % Sequential update step - gyro
    wMeas = gyroMeasurement(w_true, true); % actual measurement
    H = sensitivityMatrix("gyro", x(1:4), 0);
    R = measurementCovarianceMatrix("gyro");
    [x, P, K] = KF_MeasurementUpdate(x,wMeas,wModel,P,H,R);

    % Sequential update step - star tracker
    [starMeas, starRef, starIdxMeas] = starTrackerMeasurement(q_true,false);  % actual measurement
    [c, iMeas, iModel] = intersect(starIdxMeas, starIdxModel);  % which stars are in view in both the model and measurement?
    starModel = starModel(:,iModel); starMeas = starMeas(:,iMeas); starRef = starRef(iMeas,:);
    R = measurementCovarianceMatrix("star tracker");
    for k = 1:size(starModel,2)
        H = sensitivityMatrix("star tracker", x(1:4), starRef(k,:));
        [x, P, K] = KF_MeasurementUpdate(x,starMeas(:,k),starModel(:,k),P,H,R);
    end
end

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
    load inertiaTensors Ixx Iyy Izz

    stateTransitionMatrix = @(dt) [1 0 0 0 0.5*dt -0.5*dt 0.5*dt;
                                   0 1 0 0 0.5*dt 0.5*dt -0.5*dt;
                                   0 0 1 0 -0.5*dt 0.5*dt 0.5*dt;
                                   0 0 0 1 -0.5*dt -0.5*dt -0.5*dt;
                                   0 0 0 0 1 0 0;
                                   0 0 0 0 0 1 0;
                                   0 0 0 0 0 0 1];
    
    % Control input matrix is defined as [Mx, My, Mz]
    % Neglecting coupling between axes, we define the time update as:
    B = zeros(7,3);
    B(5,1) = dt/Ixx; B(6,2) = dt/Iyy; B(7,3) = dt/Izz;

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
    gyroOffset = biasGyro * (2*rand([1 3]) - 1);

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
    nadirRefMCI = (-1 * r ./ norm(r));
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
        gyroNoise = normrnd(0, sigmaGyro_dt^2, [3 1]);
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
    torque = torque';
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
    torque = torque';
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
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     �