%% HW 6
% NOTE: start date corresponds MAVEN start date (sep. 21, 2014)
% NOTE: we use convention where q1,q2,q3 = complex and q4 = scalar

%% Set up target attitude
% Load settings
load simulationSettings
load orbitGeometry
load spacecraftProperties

% Set up initial position
[r0, v0] = oe2mci(a,e,i,raan,argp,trueAnomaly0);        % initial position w.r.t. inertial frame
% Initial angular velocity should be equal to Mars' mean motion
w0_body = [marsMeanMotion, 0, 0];
w0_principal = rotationBodyToPrincipal * w0_body';

% Desired attitude: Z_body is aligned with sun vector, arbitrary X/Y
zb_target = calculateSunVectorMCI(0);
xb_target = cross(r0, zb_target) ./ norm(cross(r0, zb_target));
yb_target = cross(zb_target, xb_target);
dcm_initial_body = [xb_target'; yb_target'; zb_target'];
dcm_initial_principal = rotationBodyToPrincipal * dcm_initial_body;
[q10, q20, q30, q40] = dcm2Quaternion(dcm_initial_principal); % initial attitude w.r.t. inertial frame
q0 = [q10, q20, q30, q40];

y0 = [r0', v0', q0, w0_principal'];               % initial state vector
disp("Initial conditions ready!")

%% Run simulation with all torques
[t, y_out] = ode113(@trueOrbitAndAttitude, [0:tstep:t_period]' , y0, options);
pos = y_out(:,1:3); vel = y_out(:,4:6); q = y_out(:,7:10); w = y_out(:,11:13);
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
errorTheta = zeros(size(theta)); errorPhi = zeros(size(theta)); errorPsi = zeros(size(theta));

disp("Simulation done!");

%% Extract Euler angles and errors
for k = 1:size(t)
    dcm = quaternion2DCM(q(k,1),q(k,2),q(k,3),q(k,4));
    [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
    zb_target = calculateSunVectorMCI(t(k));
    xb_target = cross(r0, zb_target) ./ norm(cross(r0, zb_target));
    yb_target = cross(zb_target, xb_target);
    dcm_target_principal = rotationBodyToPrincipal * [xb_target'; yb_target'; zb_target'];
    [q1t, q2t, q3t, q4t] = dcm2Quaternion(dcm_initial_principal); % target attitude w.r.t. inerital frame
    targetDCM = quaternion2DCM(q10,q20,q30,q40); % OH Question: why does this rotation matrix look different from dcm_initial_principal?
    [theta_target, phi_target, psi_target] = dcm2Euler312(targetDCM);
    errorTheta(k) = theta(k) - theta_target;
    errorPhi(k) = phi(k) - phi_target;
    errorPsi(k) = psi(k) - psi_target;
end
disp("Euler Angles done!")

%% Extract environmental torques
torque_gg = zeros(size(pos));
torque_atm = zeros(size(pos));
torque_srp = zeros(size(pos));
torque_mag = zeros(size(pos));
f_atm = zeros(size(pos));
f_srp = zeros(size(pos));
MJD = zeros(size(t));
for k = 1:size(t)
    MJD(k) = startMJD + t(k)/secondsPerEarthDay;
    torque_gg(k,1:3) = gravityGradientTorque(pos(k,1:3)',vel(k,1:3)',q(k,1:4));
    [tdrag, fdrag] = aerodynamicTorque(pos(k,1:3)', vel(k,1:3)', q(k,1:4));
    [tsrp, fsrp] = srpTorque(t(k), pos(k,1:3)', q(k,1:4));
    torque_atm(k,1:3) = tdrag;  f_atm(k,1:3) = fdrag;
    torque_srp(k,1:3) = tsrp;   f_srp(k,1:3) = fsrp;
    torque_mag(k,1:3) = magneticTorque(pos(k,1:3)', q(k,1:4), MJD(k));
end
disp("Torques extracted! Starting plots")

%% Plot attitude results
pltTitle = sprintf("MAVEN Orbit\n Target: Sun-Pointing");
figure(); sgtitle(pltTitle);
subplot(3,3,1); hold on; title('Angular Velocity (X)');
plot(t, w(:,1), 'LineWidth', 1); grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
subplot(3,3,2); hold on; title('Angular Velocity (Y)');
plot(t, w(:,2), 'LineWidth', 1); grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
subplot(3,3,3); hold on; title('Angular Velocity (Z)');
plot(t, w(:,3), 'LineWidth', 1); grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');

subplot(3,3,5); hold on; title('Euler Angles (\theta)');
plot(t, theta, 'b'); yline(theta_target, 'r'); grid on;
xlabel('Time (s)'); ylabel('Radians'); legend('Measured',  'Target');

subplot(3,3,4); hold on; title('Euler Angles (\phi)');
plot(t, phi, 'b'); yline(phi_target, 'r'); grid on;
xlabel('Time (s)'); ylabel('Radians'); legend('Measured', 'Target');

subplot(3,3,6); hold on; title('Euler Angles (\psi)');
plot(t, psi, 'b'); yline(psi_target, 'r'); grid on;
xlabel('Time (s)'); ylabel('Radians'); legend('Measured', 'Target');

subplot(3,3,8); hold on; title('Attitude Error (\theta)');
plot(t, errorTheta, 'r'); grid on; xlabel('Time (s)'); ylabel('Radians');

subplot(3,3,7); hold on; title('Attitude Error (\phi)');
plot(t, errorPhi, 'r'); grid on; xlabel('Time (s)'); ylabel('Radians');

subplot(3,3,9); hold on; title('Attitude Error(\psi)');
plot(t, errorPsi, 'r'); grid on; xlabel('Time (s)'); ylabel('Radians');

%% Plot torques
plotTorque(t, torque_gg, theta, phi, psi, w, pos, 'Gravity Gradient');
plotTorque(t, torque_mag, theta, phi, psi, w, pos, 'Magnetic Field');
plotTorque(t, torque_atm, theta, phi, psi, w, pos, 'Aerodynamic');
plotTorque(t, torque_srp, theta, phi, psi, w, pos, 'Solar Radiation Pressure');
totalTorque = torque_gg + torque_mag + torque_atm + torque_srp;
plotTorque(t, totalTorque, theta, phi, psi, w, pos, 'Total');

%% Model Sensors
% Calculate measured attitude
q_meas = zeros(length(t),4);
w_meas = zeros(length(t),3);
for idx = 1:length(t)
    [qmt, wmt] = measureAttitude(q(idx,:), w(idx,:), t(idx));
    q_meas(idx,:) = qmt;
    w_meas(idx,:) = wmt;
end
disp("Measurements modeled");

% Convert to Euler angles
theta_meas = zeros(size(t)); phi_meas = zeros(size(t)); psi_meas = zeros(size(t));
for k = 1:size(t)
    dcm = quaternion2DCM(q_meas(k,1),q_meas(k,2),q_meas(k,3),q_meas(k,4));
    [theta_meas(k), phi_meas(k), psi_meas(k)] = dcm2Euler312(dcm);
end
disp("Angles calculated")

%% Plot results
pltTitle = sprintf("MAVEN Orbit\n Measured vs. True Attitude");
figure(); sgtitle(pltTitle);
subplot(2,2,1); hold on; title('Angular Velocity (Measured)');
plot(t, w_meas(:,1)); plot(t, w_meas(:,2)); plot(t, w_meas(:,3));
grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); legend('X', 'Y', 'Z');

subplot(2,2,2); hold on; title('Angular Velocity (Actual)');
plot(t, w(:,1)); plot(t, w(:,2)); plot(t, w(:,3));
grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); legend('X', 'Y', 'Z');

subplot(2,2,3); hold on; title('Euler Angles (Measured)'); hold on;
plot(t, phi_meas); plot(t, theta_meas); plot(t, psi_meas); grid on;
xlabel('Time (s)'); ylabel('Radians'); legend('\phi', '\theta', '\psi');

subplot(2,2,4); hold on; title('Euler Angles (Actual)'); hold on;
plot(t, phi); plot(t, theta); plot(t, psi); grid on;
xlabel('Time (s)'); ylabel('Radians'); legend('\phi', '\theta', '\psi');

%% Useful Functions: sensor modeling and attitude determination

% Problem 6: Measure Attitude
function [q_meas, w_meas] = measureAttitude(q_true, w_true, t)
    % Sensor properties
    gyroNoise = [0,0,0]; gyroOffset = [0,0,0];
    sunSensorNoise = [0,0,0]'; sunSensorOffset = [0,0,0]';
    starTrackerNoise = [0,0,0]'; starTrackerOffset = [0,0,0]';
    dcm_true = quaternion2DCM(q_true(1), q_true(2), q_true(3), q_true(4));                  % DCM from inertial to principal frame

    % Gyroscope measurement
    w_meas = w_true + gyroNoise + gyroOffset;

    % Sun sensor measurement & reference
    sunSensorReference = calculateSunVectorMCI(t);
    sunSensorMeasured = (dcm_true*sunSensorReference) + sunSensorNoise + sunSensorOffset;

    % Star tracker measurement & reference
    starCatalog = [1,0,0; 0,1,0; 0,0,1; -sqrt(2)/2, 0, sqrt(2)/2; 0, sqrt(3)/2, 1/2];       % Known star vectors in MCI frame
    starTrackerAlignment = inv(dcm_true) * [1,0,0]';     % angle of star tracker in MCI frame, assuming it points along s/c x-axis
    starVectorReference = starCatalog(starCatalog*starTrackerAlignment >= 0,:)';
    starVectorMeasured = (dcm_true*starVectorReference) + starTrackerNoise + starTrackerOffset;

    % Reconstruct attitude
    w_star = 0.8 * ones(1,length(starVectorMeasured));
    w_sun = 0.2;
    w = [w_star, w_sun];     % We trust the star tracker more than the sun sensor
    m = [starVectorMeasured, sunSensorMeasured];
    v = [starVectorReference, sunSensorReference];

    q_meas = statisticalAttitude(w, m, v);
end


% Problem 5a: Deterministic Attitude Determination
% Triad method to calculate attitude (no error cancelling)
% p1,2 --> measured vectors from star tracker and sun sensor, respectively
% r1,2 --> reference vectors for star and sun position, respectively
function attitudeDCM = deterministicAttitude(p1, p2, r1, r2)
    % create orthogonal triads
    m1 = p1; v1 = r1;
    m2 = cross(p1,p2)/norm(cross(p1,p2));   v2 = cross(r1,r2)/norm(cross(r1,r2));
    m3 = cross(m1, m2); v3 = cross(v1, v2);

    % calculate attitude dcm
    M = [m1, m2, m3];   % columns correspond to measurements
    V = [v1, v2, v3];   % columns correspond to reference vectors
    attitudeDCM = M * inv(V);   % attitude DCM
end

% Triad method to calculate attitude (with error distribution)
% y1,2 --> measured vectors from star tracker and sun sensor, respectively
% r1,2 --> reference vectors for star and sun position, respectively
function attitudeDCM = deterministicAttitudeWithErrors(y1, y2, f1, f2)
    p1 = (1/2) * (y1 + y2); p2 = (1/2) * (y1 - y2);
    r1 = (1/2) * (f1 + f2); r2 = (1/2) * (f1 - f2);
    attitudeDCM = deterministicAttitude(p1,p2,r1,r2);
end

% Problem 5b: Q-Method (Statistical Attitude Determination)
% w --> sensor weights row vector (star tracker, sun sensor)
% m --> sensor measurements as column vectors in principal axis frame (star tracker, sun sensor)
% v --> reference vectors as column vectors in MCI (star tracker, sun sensor)
function q = statisticalAttitude(w, m, v)
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

%% Useful Functions: Plots
function torque = plotTorque(t, torque, theta, phi, psi, w, pos, name)
    pltTitle = sprintf("%s Torque for MAVEN Orbit", name);
    figure(); sgtitle(pltTitle)
    subplot(3,3,1); hold on; title('Angular Velocity (X)');
    plot(t, w(:,1), 'LineWidth', 1); grid on; hold off;
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    subplot(3,3,2); hold on; title('Angular Velocity (Y)');
    plot(t, w(:,2), 'LineWidth', 1); grid on; hold off;
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    subplot(3,3,3); hold on; title('Angular Velocity (Z)');
    plot(t, w(:,3), 'LineWidth', 1); grid on; hold off;
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    subplot(3,3,4); hold on; title('Torque (X)');
    plot(t, torque(:,1), 'LineWidth', 1); ylabel('Torque (N m)'); xlabel('Time (s)'); grid on; hold off;
    subplot(3,3,5); hold on; title('Torque (Y)');
    plot(t, torque(:,2), 'LineWidth', 1); ylabel('Torque (N m)'); xlabel('Time (s)'); grid on; hold off;
    subplot(3,3,6); hold on; title('Torque (Z)');
    plot(t, torque(:,3), 'LineWidth', 1); ylabel('Torque (N m)'); xlabel('Time (s)'); grid on; hold off;
    subplot(3,3,7); hold on; title('Euler Angles');
    plot(t, theta, 'LineWidth', 1); plot(t, phi, 'LineWidth', 1); plot(t, psi, 'LineWidth', 1); legend('Theta', 'Phi', 'Psi'); grid on;  hold off;
    xlabel('Time (s)'); ylabel('Attitude (rad)');
    subplot(3,3,8); hold on; title('Orbit Radius');
    plot(t, vecnorm(pos'), 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Orbit Radius (km)'); hold off;
    subplot(3,3,9); hold on; title('Torque Magnitude');
    plot(t, vecnorm(torque'), 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Torque (N m)'); hold off;
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
    [drag, fdrag] = aerodynamicTorque(r, v, [q1 q2 q3 q4]);
    mag = magneticTorque(r, [q1 q2 q3 q4], MJD);
    [srp, fsrp] = srpTorque(t, r, [q1 q2 q3 q4]);
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
% r_mci --> s/c position in MCI frame
% q --> s/c attitude in MCI frame
function [torque] = srpTorque(t, r, q)
    
    % Constants
    torque = 0;
    load orbitGeometry R_mars
    load spacecraftProperties exposedSurfaceAreas exposedUnitVectors exposedCentroids Cs Cd

    P = 586.2/(3e8);                            % Solar pressure at Mars

    % check if in eclipse
    sunVectorMCI = calculateSunVectorMCI(t);
    r_perp = r - (sunVectorMCI * dot(r, sunVectorMCI));
    if (norm(r_perp) > R_mars || dot(r, sunVectorMCI) > 0)
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
function [torque, force] = aerodynamicTorque(r_mci, v_mci, q)

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
    force = 0;
    for j = 1:length(exposedSurfaceAreas)
        df = (-1/2) * Cdrag * rho * vnorm*vnorm * dot(vhat, exposedUnitVectors(j,1:3)) * vhat * exposedSurfaceAreas(j);
        dM = cross(exposedCentroids(j,1:3), df);
        torque = torque + dM;
        force = force + df;
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