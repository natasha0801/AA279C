%% HW 5-2
%% Set up Sun parameters
clear; clc;
load orbitGeometry;
[r_sun, v_sun] = planetEphemeris(juliandate(2014,9,21), 'Mars', 'Sun');

% Make sure reference frame is consistent w/ the rest of this program!
[a, e, i, raan, argp, trueAnomaly] = mci2oe(r_sun, v_sun);
[sunVectorMCI, sunVelocityECI] = oe2mci(a,e,i,raan,argp,trueAnomaly);
sunVectorMCI = sunVectorMCI / norm(sunVectorMCI);       % unit vector
save orbitGeometry; clear; clc;


%% Set up surface parameters
load spacecraftProperties;
exposedSurfaces = surfaces(unshaded == 1);
exposedSurfaceAreas = surfaceAreas(unshaded == 1);
exposedCentroids = centroids(unshaded == 1, :);
exposedUnitVectors = unitVectors(unshaded == 1, :);
save spacecraftProperties; clear; clc;

%% Integrate Orbit
% Set up initial conditions
% NOTE: start date corresponds MAVEN start date (sep. 21, 2014)
load simulationSettings
load orbitGeometry
load spacecraftProperties

[r0, v0] = oe2mci(a,e,i,raan,argp,trueAnomaly0);        % initial position w.r.t. inertial frame
mci2rtn_initial = inertialToRTN(r0',v0');
[q10, q20, q30, q40] = dcm2Quaternion(mci2rtn_initial); % initial attitude w.r.t. inertial frame
q0 = [q10, q20, q30, q40];
w0_body = [deg2rad(5), 0, 0];
w0_principal = rotationBodyToPrincipal * w0_body';
y0 = [r0', v0', q0, w0_principal'];               % initial state vector
disp("Initial conditions ready!")

% Run simulation
[t, y_out] = ode113(@integrateOrbitAndAttitude, [0:tstep:t_period]' , y0, options);
pos = y_out(:,1:3); vel = y_out(:,4:6); q = y_out(:,7:10); w = y_out(:,11:13);
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
disp("Simulation done!");

% Extract Euler angles
for k = 1:size(t)
    dcm = quaternion2DCM(q(k,1),q(k,2),q(k,3),q(k,4));
    [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
end
disp("Euler Angles done!")

% Extract environmental torques
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
    [tsrp, fsrp] = srpTorque(pos(k,1:3)', q(k,1:4));
    torque_atm(k,1:3) = tdrag;  f_atm(k,1:3) = fdrag;
    torque_srp(k,1:3) = tsrp;   f_srp(k,1:3) = fsrp;
    torque_mag(k,1:3) = magneticTorque(pos(k,1:3)', q(k,1:4), MJD(k));
end
disp("Torques extracted! Starting plots")

% Calculate expected maximum torques
torque_gg_max = 1.5 * mu_mars * (Izz-Ixx) ./ (vecnorm(pos')).^3;
torque_mag_max = 2*norm(m_sc)*(R_mars^3)*B0 ./ (vecnorm(pos')).^3;

%% Plot results
plotTorqueWithLimits(t, torque_gg, theta, phi, psi, w, pos, 'Gravity Gradient', torque_gg_max);
plotTorqueWithLimits(t, torque_mag, theta, phi, psi, w, pos, 'Magnetic Field', torque_mag_max);
plotTorque(t, torque_atm, theta, phi, psi, w, pos, 'Aerodynamic');
plotTorque(t, torque_srp, theta, phi, psi, w, pos, 'Solar Radiation Pressure');
totalTorque = torque_gg + torque_mag + torque_atm + torque_srp;
plotTorque(t, totalTorque, theta, phi, psi, w, pos, 'Total');

%% Useful Functions: Plotting
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


function torque = plotTorqueWithLimits(t, torque, theta, phi, psi, w, pos, name, maxTorque)
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
    plot(t, maxTorque, 'LineWidth', 1.5);
    plot(t, -1*maxTorque, 'LineWidth', 1.5); 
    legend('Torque', 'Max Limit', 'Min Limit');
    xlabel('Time (s)'); ylabel('Torque (N m)'); hold off;
end


%% Useful Functions: Integration
% Integrate orbit and attitude with ALL environmental torques
function [statedot] = integrateOrbitAndAttitude(t, state)
    disp(t)
    % Constants
    load inertiaTensors Ixx Iyy Izz
    load orbitGeometry mu_mars startMJD secondsPerEarthDay 

    % Calendar
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
    [srp, fsrp] = srpTorque(r, [q1 q2 q3 q4]);
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
function [torque, force] = srpTorque(r, q)
    
    % Constants
    torque = 0;
    force = 0;
    load orbitGeometry R_mars sunVectorMCI
    load spacecraftProperties exposedSurfaceAreas exposedUnitVectors exposedCentroids Cs Cd

    P = 586.2/(3e8);                            % Solar pressure at Mars

    % check if in eclipse
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
                force = force + df;
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


%% Useful Functions: Orbit Geometry
% Convert position/velocity to orbit elements
function [a, e, i, raan, argp, trueAnomaly] = mci2oe(r_eci, v_eci)
    mu = 398600.435507;             % Earth gravitational parameter (km^3 / s^2)
    r = norm(r_eci);                % radius magnitude (km)
    v = norm(v_eci);                % velocity magnitude (km/s)
    h = cross(r_eci, v_eci);        % specific angular momentum
    W = h/norm(h);
    i = atan(sqrt(W(1)^2 + W(2)^2)/W(3));   % inclination (radians)
    if (W(3) < 0)
        i = i + pi;
    end
    raan = atan(W(1)/(-1*W(2)));            % RAAN (radians)
    if (W(2) > 0)
        raan = raan + pi;
    end
    p = (norm(h))^2 / mu;                   % semilatus rectum (km)
    a = ((2/r) - (v^2)/mu)^(-1);            % semimajor axis (km)
    n = sqrt(mu/(a^3));                     % mean motion (1/s)
    e = sqrt(1-p/a);                        % eccentricity
    E = atan((dot(r_eci,v_eci)/(a*a*n))/(1-r/a));   % eccentric anomaly
    trueAnomaly = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));   % true anomaly (radians)
    u = atan((r_eci(3)/sin(i))/(r_eci(1)*cos(raan) + r_eci(2)*sin(raan)));  % arg. of latitude
    argp = u - trueAnomaly;     % argument of perigee (radians)

    % Convert angles to degrees
    i = rad2deg(i); raan = rad2deg(raan); argp = rad2deg(argp);
    trueAnomaly = mod(rad2deg(trueAnomaly) + 360, 360);
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
    dcm = [Q4^2+Q1^2-Q2^2-Q3^2,   2*(Q1*Q2+Q3*Q4),        2*(Q1*Q3-Q2*Q4);
         2*(Q1*Q2-Q3*Q4),       Q4^2-Q1^2+Q2^2-Q3^2,    2*(Q2*Q3+Q1*Q4);
         2*(Q1*Q3+Q2*Q4),       2*(Q2*Q3-Q1*Q4),        Q4^2-Q1^2-Q2^2+Q3^2];
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