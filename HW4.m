%% AA279C HW3
% Natasha Evans
load simulationSettings
load orbitGeometry
load inertiaTensors
t_period = 600;
Ixx = principalInertiaTensor(1,1);
Iyy = principalInertiaTensor(2,2);
Izz = principalInertiaTensor(3,3);

%% Problem 1a
% Define initial conditions
q0 = [0, 0, 0, 1];          % initial attitude
w0 = [0, 0, deg2rad(25)];   % initial angular velocity
y0 = [q0,w0];               % initial state vector

% Run simulation
[t, y_out] = ode113(@integrateAttitude, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz);
q1 = y_out(:,1); q2 = y_out(:,2); q3 = y_out(:,3); q4 = y_out(:,4);     % quaternion attitude
wx = y_out(:,5); wy = y_out(:,6); wz = y_out(:,7); 

% Quaternion to Euler angles
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
for k = 1:size(t)
    dcm = quaternion2DCM(q1(k),q2(k),q3(k),q4(k));
    [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
end

figure(); sgtitle('(1) Principal Axes Aligned with Inertial Frame')
subplot(1,2,1); title('Angular Velocity vs. Time'); hold on; grid on;
plot(t,wx); plot(t,wy); plot(t,wz); legend('wx', 'wy', 'wz');
xlabel('t (s)'); ylabel('Rad/s'); hold off;

subplot(1,2,2); title('Euler Angle vs. Time'); hold on; grid on;
plot(t, theta); plot(t, phi); plot(t, psi); legend('Theta', 'Phi', 'Psi');
xlabel('t (s)'); ylabel('Radians'); hold off;

%% Problem 1b
% Calculate r and v in inertial frame
[r0, v0] = oe2mci(a,e,k,raan,argp,trueAnomaly0);
y0 = [r0;v0];
[t, y_out] = ode113(@orbitingbody , [0:tstep:t_period]' , y0 , options );

% Set up initial conditions
rvecs = y_out(:,1:3); vvecs = y_out(:,4:6);
rvec0 = rvecs(1,:); vvec0 = vvecs(1,:);
mci2rtn_initial = inertialToRTN(rvec0,vvec0);
[q10, q20, q30, q40] = dcm2Quaternion(mci2rtn_initial);
q0 = [q10,q20,q30,q40];          % initial attitude w.r.t. inertial frame
w0 = [0, 0, deg2rad(15)];        % initial angular velocity
y0 = [q0,w0];                    % initial state vector

% Run attitude simulation
[t, y_out] = ode113(@integrateAttitude, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz);
q1 = y_out(:,1); q2 = y_out(:,2); q3 = y_out(:,3); q4 = y_out(:,4);     % quaternion attitude in inertial frame
wx = y_out(:,5); wy = y_out(:,6); wz = y_out(:,7);                      % principal axes angular velocity in inertial frame
w_rtn = zeros(size([wx, wy, wz]));      % allocate space for angular rates in RTN frame

% Convert quaternion attitude in inertial frame to Euler Angles in RTN frame
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
for k = 1:size(t)
    rvec = rvecs(k,:); vvec = vvecs(k,:);                           % s/c position and velocity (inertial frame)
    inertial2Principal = quaternion2DCM(q1(k),q2(k),q3(k),q4(k));   % DCM: Principal Axis rel. to inertial
    inertial2RTN = inertialToRTN(rvec, vvec);                       % DCM: RTN rel. to inertial
    RTN2Principal = inertial2Principal * inertial2RTN';
    [theta(k), phi(k), psi(k)] = dcm2Euler312(RTN2Principal);
    w_rtn(k,:) = RTN2Principal' * [wx(k), wy(k), wz(k)]';
end

figure();  sgtitle('(1) Principal Axes Aligned with RTN Frame')
xlabel('t (s)'); ylabel('Rad/s'); hold off;

subplot(1,2,1); title('Angular Velocity vs. Time - RTN Frame'); hold on; grid on;
plot(t, w_rtn(:,1)); plot(t, w_rtn(:,2)); plot(t, w_rtn(:,3)); legend('wx', 'wy', 'wz');
xlabel('t (s)'); ylabel('Radians'); hold off;

subplot(1,2,2); title('Euler Angles - RTN Frame'); hold on; grid on;
plot(t, theta); plot(t, phi); plot(t, psi); legend('Theta', 'Phi', 'Psi');
xlabel('t (s)'); ylabel('Radians'); hold off;

%% Problem 2: stability tests
testAxes=['X','Y','Z'];
rotationRate = deg2rad(25);
perturbation = deg2rad(5);

for rotationAxis = 1:3
    % Initial condition - no perturbation
    q0 = [0,0,0,1];         % initial attitude w.r.t. inertial frame
    w0 = [0, 0, 0];         % initial angular velocity
    w0(rotationAxis) = rotationRate;    % nonzero angular velocity about a principal axis
    y0 = [q0,w0];           % initial state vector
    
    % Run attitude simulation - no perturbation
    [t, y_out] = ode113(@integrateAttitude, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz);
    q1 = y_out(:,1); q2 = y_out(:,2); q3 = y_out(:,3); q4 = y_out(:,4);     % quaternion attitude in inertial frame
    wx = y_out(:,5); wy = y_out(:,6); wz = y_out(:,7);                      % principal axes angular velocity in inertial frame

    % Quaternion to Euler angles - multiple sequence options in case of
    % singularity
    theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
    use313 = false;
    for k = 1:size(t)
        dcm = quaternion2DCM(q1(k),q2(k),q3(k),q4(k));
        [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
        if imag(theta(k)) ~= 0 || imag(phi(k)) ~= 0 || imag(psi(k)) ~= 0
            use313 = true;
        end
    end
    if (use313)
        for k = 1:size(t)
            dcm = quaternion2DCM(q1(k),q2(k),q3(k),q4(k));
            [theta(k), phi(k), psi(k)] = dcm2Euler313(dcm);
        end
    end

    % Initial condition - with perturbation
    w0 = w0 + perturbation;
    y0 = [q0,w0];   % initial state vector
    
    % Run attitude simulation - with perturbation
    [t, y_out] = ode113(@integrateAttitude, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz);
    q1p = y_out(:,1); q2p = y_out(:,2); q3p = y_out(:,3); q4p = y_out(:,4);     % quaternion attitude in inertial frame
    wxp = y_out(:,5); wyp = y_out(:,6); wzp = y_out(:,7);                      % principal axes angular velocity in inertial frame

    % Quaternion to Euler angles
    thetap = zeros(size(t)); phip = zeros(size(t)); psip = zeros(size(t));
    use313 = false;
    for k = 1:size(t)
        dcm = quaternion2DCM(q1p(k),q2p(k),q3p(k),q4p(k));
        [thetap(k), phip(k), psip(k)] = dcm2Euler312(dcm);
        if imag(thetap(k)) ~= 0 || imag(phip(k)) ~= 0 || imag(psip(k)) ~= 0
            use313 = true;
        end
    end
    if (use313)
        for k = 1:size(t)
            dcm = quaternion2DCM(q1p(k),q2p(k),q3p(k),q4p(k));
            [thetap(k), phip(k), psip(k)] = dcm2Euler313(dcm);
        end
    end
    
    plttitle = sprintf('(2) Single-Spin Stability Test (%s Axis)', testAxes(rotationAxis));
    figure(); sgtitle(plttitle);
    subplot(2,2,1); title('Angular Velocity (No Perturbation)'); hold on; grid on;
    plot(t,wx); plot(t,wy); plot(t,wz); legend('wx', 'wy', 'wz');
    xlabel('t (s)'); ylabel('Rad/s'); hold off;
    
    subplot(2,2,2); title('Euler Angles (No Perturbation)'); hold on; grid on;
    plot(t, theta); plot(t, phi); plot(t, psi); legend('Theta', 'Phi', 'Psi');
    xlabel('t (s)'); ylabel('Radians'); hold off;

    subplot(2,2,3); title('Angular Velocity (Perturbed)'); hold on; grid on;
    plot(t,wxp); plot(t,wyp); plot(t,wzp); legend('wx', 'wy', 'wz');
    xlabel('t (s)'); ylabel('Rad/s'); hold off;
    
    subplot(2,2,4); title('Euler Angles (Perturbed)'); hold on; grid on;
    plot(t, thetap); plot(t, phip); plot(t, psip); legend('Theta', 'Phi', 'Psi');
    xlabel('t (s)'); ylabel('Radians'); hold off;
end 

%% Problem (3a) Program Euler Equations and RW specs - see integrateAttitudeDualSpin
Ir = 8e-3;              % RW moment of inertia (kg m^2)
wr0 = rpm2Rads(6500);   % RW max speed
wz0 = deg2rad(15);
perturbation = deg2rad(2);

%% (3b) Integrate, verify integration accuracy
% Initial condition - no perturbation
q0 = [0,0,0,1];                % initial attitude w.r.t. inertial frame
w0 = [0, 0, wz0, wr0];         % initial angular velocity
y0 = [q0,w0];                  % initial state vector

% Run attitude simulation - no perturbation
[t, y_out] = ode113(@integrateAttitudeDualSpin, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz, Ir, 3);
q1 = y_out(:,1); q2 = y_out(:,2); q3 = y_out(:,3); q4 = y_out(:,4);     % quaternion attitude in inertial frame
wx = y_out(:,5); wy = y_out(:,6); wz = y_out(:,7); wr = y_out(:,8);     % principal axes angular velocity in inertial frame

% Quaternion to Euler angles - 312 should work for rotation about Z-axis
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
for k = 1:size(t)
    dcm = quaternion2DCM(q1(k),q2(k),q3(k),q4(k));
    [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
end

% Plot Results
figure();
subplot(2,1,1); title('Dual-Spin Angular Velocity (No Perturbation)'); hold on; grid on;
plot(t,wx); plot(t,wy); plot(t,wz); plot(t,wr); legend('wx', 'wy', 'wz', 'wr');
xlabel('t (s)'); ylabel('Rad/s'); hold off;

subplot(2,1, 2); title('Dual-Spin Euler Angles (No Perturbation)'); hold on; grid on;
plot(t, theta); plot(t, phi); plot(t, psi); legend('Theta', 'Phi', 'Psi');
xlabel('t (s)'); ylabel('Radians'); hold off;

% Verify accuracy by checking angular momentum & energy conservation
L_principal = [Ixx.*wx, Iyy.*wy, Izz.*wz + Ir.*wr]; 
w_principal = [wx, wy, wz];
L_inertial = zeros(size(L_principal));
w_inertial = zeros(size(w_principal));
for k = 1:length(L_principal)
    Q4 = q4(k); Q1 = q1(k); Q2 = q2(k); Q3 = q3(k);
    R = quaternion2DCM(Q1,Q2,Q3,Q4);
    L_inertial(k,:) = R' * L_principal(k,:)';
    w_inertial(k,:) = R' * w_principal(k,:)';
end

% Conserved angular momentum
figure(); sgtitle('Conservation Laws')
subplot(1,2,1); grid on; hold on;
plot(t, L_inertial(:,1), 'LineWidth', 1); 
plot(t, L_inertial(:,2), 'LineWidth', 1);
plot(t, L_inertial(:,3), 'LineWidth', 1);
legend('Lx', 'Ly', 'Lz');
xlabel('t'); ylabel('kg*m*m/s'); title('Angular Momentum Vector in Inertial Frame');
hold off;

% Conserved rotational KE
% Use 2T = w dot L = const (lecture 3 slide 6)
energy = (1/2) * dot(L_inertial',w_inertial');
subplot(1,2,2);
plot(t, energy, 'LineWidth', 1); grid on;
xlabel('Time (s)'); ylabel('Rotational KE');
ylim([0 800]);
title('Rotational Kinetic Energy = (1/2) * w * L)');

%% 3(c): Stability Test about Z axis
% Initial condition - with perturbation
q0 = [0,0,0,1];                % initial attitude w.r.t. inertial frame
w0 = [0, 0, wz0, wr0] + perturbation;   % initial angular velocity (perturbed)
y0 = [q0,w0];                   % initial state vector

% Run attitude simulation - with perturbation
[t, y_out] = ode113(@integrateAttitudeDualSpin, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz, Ir, 3);
q1p = y_out(:,1); q2p = y_out(:,2); q3p = y_out(:,3); q4p = y_out(:,4);     % quaternion attitude in inertial frame
wxp = y_out(:,5); wyp = y_out(:,6); wzp = y_out(:,7); wrp = y_out(:,8);     % principal axes angular velocity in inertial frame

% Quaternion to Euler angles
thetap = zeros(size(t)); phip = zeros(size(t)); psip = zeros(size(t));
for k = 1:size(t)
    dcm = quaternion2DCM(q1p(k),q2p(k),q3p(k),q4p(k));
    [thetap(k), phip(k), psip(k)] = dcm2Euler312(dcm);
end

% Plot Results
figure(); sgtitle('3(c) Perturbed S/C with Initial Spin about Z-Axis')
subplot(2,1,1); title('Angular Velocity (Perturbed)'); hold on; grid on;
plot(t,wxp); plot(t,wyp); plot(t,wzp); legend('wx', 'wy', 'wz');
xlabel('t (s)'); ylabel('Rad/s'); hold off;

subplot(2,1,2); title('Euler Angles (Perturbed)'); hold on; grid on;
plot(t, thetap); plot(t, phip); plot(t, psip); legend('Theta', 'Phi', 'Psi');
xlabel('t (s)'); ylabel('Radians'); hold off;

%% 3(d): Stability about intermediate axis
% Initial condition - with perturbation
Ir = 0.08; wr0 = rpm2Rads(12000);
q0 = [0,0,0,1];                % initial attitude w.r.t. inertial frame
w0 = [0, wz0, 0, wr0] + perturbation;   % initial angular velocity (perturbed)
y0 = [q0,w0];                   % initial state vector

% Run attitude simulation - with perturbation
[t, y_out] = ode113(@integrateAttitudeDualSpin, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz, Ir, 2);
q1p = y_out(:,1); q2p = y_out(:,2); q3p = y_out(:,3); q4p = y_out(:,4);     % quaternion attitude in inertial frame
wxp = y_out(:,5); wyp = y_out(:,6); wzp = y_out(:,7); wrp = y_out(:,8);     % principal axes angular velocity in inertial frame

% Quaternion to Euler angles
thetap = zeros(size(t)); phip = zeros(size(t)); psip = zeros(size(t));
for k = 1:size(t)
    dcm = quaternion2DCM(q1p(k),q2p(k),q3p(k),q4p(k));
    [thetap(k), phip(k), psip(k)] = dcm2Euler312(dcm);
end

% Plot Results
figure(); sgtitle('3(d) Perturbed S/C with Initial Spin about Y-Axis')
subplot(2,1,1); title('Angular Velocity (Perturbed)'); hold on; grid on;
plot(t,wxp); plot(t,wyp); plot(t,wzp); legend('wx', 'wy', 'wz');
xlabel('t (s)'); ylabel('Rad/s'); hold off;

subplot(2,1,2); title('Euler Angles (Perturbed)'); hold on; grid on;
plot(t, thetap); plot(t, phip); plot(t, psip); legend('Theta', 'Phi', 'Psi');
xlabel('t (s)'); ylabel('Radians'); hold off;

%% 3(e) Stability about X-Axis
Ir = 0.008; wr0 = rpm2Rads(6500);
q0 = [0,0,0,1];                % initial attitude w.r.t. inertial frame
w0 = [wz0, 0, 0, wr0] + perturbation;   % initial angular velocity (perturbed)
y0 = [q0,w0];                   % initial state vector

% Run attitude simulation - with perturbation
[t, y_out] = ode113(@integrateAttitudeDualSpin, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz, Ir, 1);
q1p = y_out(:,1); q2p = y_out(:,2); q3p = y_out(:,3); q4p = y_out(:,4);     % quaternion attitude in inertial frame
wxp = y_out(:,5); wyp = y_out(:,6); wzp = y_out(:,7); wrp = y_out(:,8);     % principal axes angular velocity in inertial frame

% Quaternion to Euler angles
thetap = zeros(size(t)); phip = zeros(size(t)); psip = zeros(size(t));
for k = 1:size(t)
    dcm = quaternion2DCM(q1p(k),q2p(k),q3p(k),q4p(k));
    [thetap(k), phip(k), psip(k)] = dcm2Euler312(dcm);
end

% Plot Results
figure(); sgtitle('3(e) Perturbed S/C with Initial Spin about X-Axis')
subplot(2,1,1); title('Angular Velocity (Perturbed)'); hold on; grid on;
plot(t,wxp); plot(t,wyp); plot(t,wzp); legend('wx', 'wy', 'wz');
xlabel('t (s)'); ylabel('Rad/s'); hold off;

subplot(2,1,2); title('Euler Angles (Perturbed)'); hold on; grid on;
plot(t, thetap); plot(t, phip); plot(t, psip); legend('Theta', 'Phi', 'Psi');
xlabel('t (s)'); ylabel('Radians'); hold off;

%% Problem 4: Gravity Gradient
% 4(a): Remove Rotor.
% 4(b): Program gravity gradient torque (see gravityGradientTorque,
% intOrbitAttitudGravityGradient
mu_mars = 4.28284e4;
Ixx = principalInertiaTensor(1,1); Iyy = principalInertiaTensor(2,2); Izz = principalInertiaTensor(3,3);

% 4(c): Verify accuracy of gravity gradient torque
disp("---- Problem 4(c): Testing Gravity Gradient ----")
[r_lo, v_lo] = oe2mci(a,e,i,raan,argp,0);               % Periapsis
[r_hi, v_hi] = oe2mci(a,e,i,raan,argp,180);             % Apoapsis
[r_mid, v_mid] = oe2mci(a,e,i,raan,argp,90);           % Between Periapsis and Apoapsis
mci2rtn_90 = inertialToRTN(r_mid', v_mid');             % RTN frame partway through orbit
q1 = 0.0289; q2 = 0.1076; q3 = 0.1925; q4 = 0.9749;     % Arbitrary initial attitude (not aligned with RTN)
[q1m, q2m,q3m,q4m] = dcm2Quaternion(mci2rtn_90);        % Attitude aligned with RTN partway through orbit
torque_lo = simplifiedGravityGradient(r_lo(1), r_lo(2), r_lo(3), v_lo(1), v_lo(2), v_lo(3), q1, q2, q3, q4, mu_mars,Ixx,Iyy,Izz);
torque_hi = simplifiedGravityGradient(r_hi(1), r_hi(2), r_hi(3), v_hi(1), v_hi(2), v_hi(3), q1, q2, q3, q4, mu_mars,Ixx,Iyy,Izz);
torque_mid = simplifiedGravityGradient(r_mid(1), r_mid(2), r_mid(3), v_mid(1), v_mid(2), v_mid(3), q1, q2, q3, q4, mu_mars,Ixx,Iyy,Izz);
torque_rtn = simplifiedGravityGradient(r_mid(1), r_mid(2), r_mid(3), v_mid(1), v_mid(2), v_mid(3), q1m, q2m, q3m, q4m, mu_mars,Ixx,Iyy,Izz);
fprintf("Torque at Periapsis (Low Altitude): %.4e N/m\n", norm(torque_lo));
fprintf("Torque at 90 degrees True Anomaly (Medium Altitude): %.4e N/m\n", norm(torque_mid));
fprintf("Torque at Apoapsis (High Altitude): %.4e N/n\n", norm(torque_hi));
fprintf("Torque Aligned with RTN Frame %.4e N/m\n", norm(torque_rtn));

% 4(d): Numerically integrate w/ body axes aligned with orbital frame for
% a circular orbit
e_circ = 0;
[r0, v0] = oe2mci(a,e_circ,i,raan,argp,trueAnomaly0);        % initial position w.r.t. inertial frame
meanMotion = sqrt(mu_mars/(a*a*a));
mci2rtn_initial = inertialToRTN(r0',v0');
[q10, q20, q30, q40] = dcm2Quaternion(mci2rtn_initial);
q0 = [q10,q20,q30,q40];                % initial attitude w.r.t. inertial frame
w0 = [0, 0, meanMotion];               % initial angular velocity
y0 = [r0', v0', q0, w0];               % initial state vector

% Run attitude simulation
[t, y_out] = ode113(@intOrbitAttitudeGravGradient, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz, mu_mars);
pos = y_out(:,1:3); vel = y_out(:,4:6); q = y_out(:,7:10); w = y_out(:,11:13);
torque = zeros(size(t));
for k = 1:size(t)
    torquevec = simplifiedGravityGradient(pos(k,1),pos(k,2),pos(k,3),vel(k,1),vel(k,2),vel(k,3),q(k,1),q(k,2),q(k,3),q(k,4),mu_mars,Ixx,Iyy,Izz);
    torque(k) = norm(torquevec);
end

% Quaternion to Euler angles
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
for k = 1:size(t)
    dcm = quaternion2DCM(q(k,1),q(k,2),q(k,3),q(k,4));
    [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
end

% Plot results
figure(); sgtitle('Gravity Gradient for Circular Orbit, Initially Aligned with RTN Frame')
subplot(2,2,1); hold on; title('Orbit Radius');
plot(t, vecnorm(pos'), 'LineWidth', 1.5); grid on; ylim([a-10, a+10]);
xlabel('Time (s)'); ylabel('Orbit Radius (km)'); hold off;
subplot(2,2,2); hold on; title('Gravity Gradient');
plot(t, torque, 'LineWidth', 1.5); grid on; ylim([0, 1e-4])
ylabel('Gravity Gradient Torque (N m)'); xlabel('Time (s)'); hold off;
subplot(2,2,3); hold on; title('Angular Velocity');
plot(t, w(:,1)); plot(t,w(:,2)); plot(t,w(:,3)); 
legend('X', 'Y', 'Z'); grid on;  hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
subplot(2,2,4); hold on; title('Euler Angles');
plot(t, theta); plot(t, phi); plot(t, psi); legend('Theta', 'Phi', 'Psi'); grid on;  hold off;
xlabel('Time (s)'); ylabel('Attitude (rad)');


% 4(e) Numerically integrate, relevant to project, over multiple orbits!
T_orbit = 2*pi/meanMotion; long_tstep = 3 * tstep;
[r0, v0] = oe2mci(a,e,i,raan,argp,trueAnomaly0);        % initial position w.r.t. inertial frame
mci2rtn_initial = inertialToRTN(r0',v0');
[q10, q20, q30, q40] = dcm2Quaternion(mci2rtn_initial); % initial attitude w.r.t. inertia frame

% Want slow rotation about x- body axis such that spacecraft points at sun
marsYear = 687 * 24 * 60 * 60;         % Mars year duration (s)
marsMeanMotion = 2 * pi / marsYear;    % Mars mean motion (rad/s)
q0 = [q10,q20,q30,q40];                % initial attitude w.r.t. inertial frame
w0_body = [marsMeanMotion, 0, 0];           % initial angular velocity
w0_principal = rotationBodyToPrincipal * w0_body';
y0 = [r0', v0', q0, w0_principal'];               % initial state vector

% Run attitude simulation
[t, y_out] = ode113(@intOrbitAttitudeGravGradient, [0:long_tstep:2*T_orbit]' , y0, options, Ixx, Iyy, Izz, mu_mars);
pos = y_out(:,1:3); vel = y_out(:,4:6); q = y_out(:,7:10); w = y_out(:,11:13);
torque = zeros(size(pos));
for k = 1:size(t)
    torquevec = simplifiedGravityGradient(pos(k,1),pos(k,2),pos(k,3),vel(k,1),vel(k,2),vel(k,3),q(k,1),q(k,2),q(k,3),q(k,4),mu_mars,Ixx,Iyy,Izz);
    torque(k,1) = torquevec(1);
    torque(k,2) = torquevec(2);
    torque(k,3) = torquevec(3);
end

% Quaternion to Euler angles
theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
for k = 1:size(t)
    dcm = quaternion2DCM(q(k,1),q(k,2),q(k,3),q(k,4));
    [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
end

% Plot results
figure(); sgtitle('Gravity Gradient for Sun-Pointing MAVEN Orbit')
subplot(3,3,1); hold on; title('Angular Velocity (X)');
plot(t, w(:,1), 'LineWidth', 1); grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
subplot(3,3,2); hold on; title('Angular Velocity (Y)');
plot(t, w(:,2), 'LineWidth', 1); grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
subplot(3,3,3); hold on; title('Angular Velocity (Z)');
plot(t, w(:,3), 'LineWidth', 1); grid on; hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
subplot(3,3,4); hold on; title('Gravity Gradient (X)');
plot(t, torque(:,1), 'LineWidth', 1); ylabel('Torque (N m)'); xlabel('Time (s)'); grid on; hold off;
subplot(3,3,5); hold on; title('Gravity Gradient (Y)');
plot(t, torque(:,2), 'LineWidth', 1); ylabel('Torque (N m)'); xlabel('Time (s)'); grid on; hold off;
subplot(3,3,6); hold on; title('Gravity Gradient (Z)');
plot(t, torque(:,3), 'LineWidth', 1); ylabel('Torque (N m)'); xlabel('Time (s)'); grid on; hold off;
subplot(3,3,7); hold on; title('Euler Angles');
plot(t, theta, 'LineWidth', 1); plot(t, phi, 'LineWidth', 1); plot(t, psi, 'LineWidth', 1); legend('Theta', 'Phi', 'Psi'); grid on;  hold off;
xlabel('Time (s)'); ylabel('Attitude (rad)');
subplot(3,3,8); hold on; title('Orbit Radius');
plot(t, vecnorm(pos'), 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Orbit Radius (km)'); hold off;
subplot(3,3,9); hold on; title('Gravity Gradient Magnitude');
plot(t, vecnorm(torque'), 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Torque (N m)'); hold off;

%% Useful Functions

% Return gravity gradient torque in princpal axes (lec7, slide 9) 
function torque = simplifiedGravityGradient(rx,ry,rz,vx,vy,vz,q1,q2,q3,q4,mu,Ix,Iy,Iz)
    rvec = [rx, ry, rz]; vvec=[vx,vy,vz];                   % pos/vel in inertial frame
    R = norm(rvec);
    inertial2Principal = quaternion2DCM(q1,q2,q3,q4);       % DCM: Principal Axis rel. to inertial
    inertial2RTNDCM = inertialToRTN(rvec, vvec);            % DCM: RTN rel. to inertial
    RTN2Principal = inertial2Principal * inertial2RTNDCM';  % DCM: principal axis rel. to RTN
    r_rtn = [R,0,0];                                        % position vector in RTN frame
    r_xyz = RTN2Principal * r_rtn';                         % position vector in XYZ frame
    cx = r_xyz(1)/R; cy = r_xyz(2)/R; cz = r_xyz(3)/R;
    torque = (3*mu/(R^3))*[(Iz-Iy)*cy*cz; (Ix-Iz)*cz*cx; (Iy-Ix)*cx*cy];
end

% Integrate orbit and attitude with gravity gradient
function [statedot] = intOrbitAttitudeGravGradient(t, state, Ixx, Iyy, Izz, mu_mars)

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

    % Calculate angular acceleration with gravity gradient effects
    torque = simplifiedGravityGradient(r(1),r(2),r(3),v(1),v(2),v(3),q1,q2,q3,q4,mu_mars,Ixx,Iyy,Izz);
    statedot(11) = (torque(1)-(Izz - Iyy)*(wy*wz))/Ixx;
    statedot(12) = (torque(2)-(Ixx - Izz)*(wz*wx))/Iyy;
    statedot(13) = (torque(3)-(Iyy - Ixx)*(wx*wy))/Izz;

end

% RPM to rad/s
function rads = rpm2Rads(rpm)
    rads = (rpm / 60) * 2 * pi;
end

% Euler equations with generic momentum wheel aligned with Z-axis
function [statedot] = integrateAttitudeDualSpin(t, state, Ixx, Iyy, Izz, Ir, ax)

    % State is [q1 q2 q3 q4 wx wy wz wr]
    q1 = state(1); q2 = state(2); q3 = state(3); q4 = state(4);
    wx = state(5); wy = state(6); wz = state(7); wr = state(8);
    statedot = zeros (8,1);

    % Calculate rate of change of attitude, expressed as quaternion
    statedot(1:4) = quaternionEOM([q1 q2 q3 q4], wx, wy, wz);

    % Calculate angular acceleration from Euler, in rad/s/s
    if ax == 3
        statedot(5) = (1/Ixx) * (-(Ir*wr*wy) - (Izz - Iyy)*wy*wz);
        statedot(6) = (1/Iyy) * (Ir*wr*wx - (Ixx - Izz)*wz*wx);
        statedot(7) = 0;
        statedot(8) = 0;
    elseif ax == 1
        statedot(5) = 0;
        statedot(6) = (1/Iyy) * (-Ir*wr*wz - (Ixx - Izz)*wz*wx);
        statedot(7) = (1/Izz) * (Ir*wr*wx - (Iyy - Ixx)*wx*wy);
        statedot(8) = 0;
    elseif ax == 2
        statedot(5) = (1/Ixx) * (Ir*wr*wz - (Izz - Iyy)*wy*wz);
        statedot(6) = 0;
        statedot(7) = (1/Izz) * (-Ir*wr*wx - (Iyy - Ixx)*wx*wy);
        statedot(8) = 0;
    end
end

% Calculate rotation matrix from inertial to RTN
function [dcm] = inertialToRTN(rvec,vvec)
    
    Rhat = rvec ./ norm(rvec);
    Nhat = cross(rvec,vvec) ./ norm(cross(rvec,vvec));
    That = cross(Nhat, Rhat);
    
    % Initial conditions
    dcm = [Rhat; That; Nhat];
end

% DCM to Euler
function [theta, phi, psi] = dcm2Euler312(A)
    theta = asin(A(3,2));
    phi = atan2(A(1,2), A(2,2));
    psi = atan2(A(3,1),A(3,3));
end

% DCM to Euler (313)
function [theta, phi, psi] = dcm2Euler313(A)
    theta = acos(A(3,3));
    phi = -atan2(A(3,1),A(3,2));
    psi = atan2(A(1,3),A(2,3));
end

% Quaternion to DCM - Q4 is scalar component
function [dcm] = quaternion2DCM(Q1,Q2,Q3,Q4) 
    dcm = [Q4^2+Q1^2-Q2^2-Q3^2,   2*(Q1*Q2+Q3*Q4),        2*(Q1*Q3-Q2*Q4);
         2*(Q1*Q2-Q3*Q4),       Q4^2-Q1^2+Q2^2-Q3^2,    2*(Q2*Q3+Q1*Q4);
         2*(Q1*Q3+Q2*Q4),       2*(Q2*Q3-Q1*Q4),        Q4^2-Q1^2-Q2^2+Q3^2];
end

% DCM to quaternion
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

% Integrate Euler and quaternion equations from arbitrary initial conditions
function [statedot] = integrateAttitude(t, state, Ixx, Iyy, Izz)

    % State is [q1 q2 q3 q4 wx wy wz]
    q1 = state(1); q2 = state(2); q3 = state(3); q4 = state(4);
    wx = state(5); wy = state(6); wz = state(7);
    statedot = zeros (7,1);

    % Calculate rate of change of attitude, expressed as quaternion
    statedot(1:4) = quaternionEOM([q1 q2 q3 q4], wx, wy, wz);

    % Calculate angular acceleration from Euler, in rad/s/s
    statedot(5) = -(Izz - Iyy)*(wy*wz)/Ixx;
    statedot(6) = -(Ixx - Izz)*(wz*wx)/Iyy;
    statedot(7) = -(Iyy - Ixx)*(wx*wy)/Izz;
end

% ODE integration function for orbit propagation (2-body problem)
function [statedot] = orbitingbody(t , state)
    % Returns derivative of state of orbiting body
    % State is [ rx ry rz vx vy vz ] ' in [km] and [km/s]
    mu_mars = 4.28284e4;    % gravitational parameter (km^3/s^2)
    r = state (1:3);
    v = state (4:6);
    normr = norm (r) ;
    rhat = r / normr ;
    % Acceleration due to 1/ r ^2 gravity
    acc = - (mu_mars/(normr * normr)) * rhat ;
    statedot = zeros (6 , 1) ;
    statedot (1:3) = v ;
    statedot (4:6) = acc ;
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