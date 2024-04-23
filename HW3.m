%% AA279C HW3
% Natasha Evans
load simulationSettings

%% Problem 1: Numerical simulation of axial-symmetric satellite
% Define initial conditions and output vector
wx0 = deg2rad(15); wy0 = deg2rad(-20); wz0 = deg2rad(-12);  % initial rotation rates (deg/s)
y0 = zeros(6,1);        % zero initial angle
y0(4) = wx0;   % initial rotation rate about x axis (rad/s)
y0(5) = wy0;   % initial rotation rate about y axis (rad/s)
y0(6) = wz0;   % initial rotation rate about z axis (rad/s)

% Simulation outputs
[t, y_out] = ode113(@axialSymmetricEuler, [0:tstep:t_period]' , y0, options);
thetax_sim = y_out(:,1); thetay_sim = y_out(:,2); thetaz_sim = y_out(:,3);      % angle (rad)
wx_sim = y_out(:,4); wy_sim = y_out(:,5); wz_sim = y_out(:,6);                  % rotation rate (rad/s)

%% Problem 2: Analytical solution for axial symmetric satellite
% moments of inertia about principal axes
Ixx = 3000;
Iyy = 3000;
Izz = 5000;
lambda = (Izz-Iyy)*wz0/Ixx;
sampleTimes = 0:tstep:t_period;
numSamples = length(sampleTimes);

% Z-axis: constant rotation rate
wz = ones(numSamples,1) * wz0;
wdotz = zeros(numSamples,1);

% X- and Y- axis are coupled
wdotx = ones(numSamples,1); wdoty = ones(numSamples,1); 
wx = ones(numSamples,1); wy = ones(numSamples,1);
wx(1) = wx0; wy(1) = wy0;
for i = 1:numSamples
    wdotx(i) = -lambda*wy(i); 
    wdoty(i) = lambda*wx(i);
    if (i < numSamples)
        wx(i + 1) = wx(i) + tstep*wdotx(i);
        wy(i + 1) = wy(i) + tstep*wdoty(i);
    end
end

%% Problem 3: Compare results
error_x = wx_sim - wx;
error_y = wy_sim - wy;
error_z = wz_sim - wz;

% Angular velocity components
figure();
subplot(3,1,1); hold on; grid on;
titleStr = sprintf('Analytical vs. Numerical Simulation of Axially Symmetric S/C\n Timestep=%.1e s, RelTol=%.1e, AbsTol=%.1e\n',tstep,reltol,abstol);
title(titleStr);
plot(t,wx); plot(t,wx_sim); plot(t, error_x);
legend('Analytical', 'ODE113', 'Error');
xlabel('t (s)'); ylabel('wx (rad/s)'); hold off;

subplot(3,1,2); hold on; grid on;
plot(t,wy); plot(t,wy_sim); plot(t, error_y);
legend('Analytical', 'ODE113', 'Error');
xlabel('t (s)'); ylabel('wy (rad/s)'); hold off;

subplot(3,1,3); hold on; grid on;
plot(t,wz); plot(t,wz_sim); plot(t, error_z);
legend('Analytical', 'ODE113', 'Error');
xlabel('t (s)'); ylabel('wz (rad/s)'); hold off;

% Check angular velocity vector
figure();
subplot(1,2,1); hold on; grid on;
title('Angular Velocity Vector');
plot(t, sqrt(wx_sim.^2 + wy_sim.^2)); plot(t, abs(wz_sim));
legend('||w_{xy}||', '|w_z|');
xlabel('t (s)'); ylabel('rad/s'); hold off;

% Check angular momentum vector
subplot(1,2,2); hold on; grid on;
title('Angular Momentum Vector');
plot(t, sqrt((Ixx.*wx_sim).^2 + (Iyy.*wy_sim).^2)); plot(t, abs(Izz.*wz_sim));
legend('L_{xy}', 'L_z');
xlabel('t (s)'); ylabel('rad/s'); hold off;


%% Problem 3: Kinematic EOM for quaternion
% See quaternionEOM under Useful Functions.

%% Problem 4: Kinematic EOM for DCM
% See dcmEOM under Unseful Functions.

%% Problem 5: Integrate kinematic and Euler eqns for spacecraft (quaternions)
% Define initial conditions and output vector
wx0 = deg2rad(25); wy0 = deg2rad(-20); wz0 = deg2rad(-12);  % initial rotation rates (deg/s)
y0q = zeros(7,1);
y0q(1) = 0; y0q(2) = 0; y0q(3) = 0; y0q(4) =    1;     % zero initial angle
y0q(5) = wx0;   % initial rotation rate about x axis (rad/s)
y0q(6) = wy0;   % initial rotation rate about y axis (rad/s)
y0q(7) = wz0;   % initial rotation rate about z axis (rad/s)

% Simulation outputs - quaternion
[t, y_out] = ode113(@eulerWithQuaternion, [0:tstep:t_period]' , y0q, options);
q1 = y_out(:,1); q2 = y_out(:,2); q3 = y_out(:,3); q4 = y_out(:,4);     % quaternion attitude
wx = y_out(:,5); wy = y_out(:,6); wz = y_out(:,7);                 % rotation rate (rad/s)

% Plot results
figure();
subplot(1,2,1); title('Angular Velocity vs. Time'); hold on; grid on;
plot(t,wx); plot(t,wy); plot(t,wz); legend('wx', 'wy', 'wz');
xlabel('t (s)'); ylabel('Rad/s'); hold off;

subplot(1,2,2); title('Attitude vs. Time'); hold on; grid on;
plot(t, q1); plot(t,q2); plot(t,q3); plot(t,q4); legend('q1', 'q2', 'q3', 'q4');
xlabel('t (s)'); ylabel('Quaternion Component'); hold off;


%% Problem 7ab: angular momentum and angular velocity vector in inertial coordinates
% Vectors in principal axis frame
load inertiaTensors
Ixx = principalInertiaTensor(1,1); Iyy = principalInertiaTensor(2,2); Izz = principalInertiaTensor(3,3);
L_principal = [Ixx.*wx, Iyy.*wy, Izz.*wz]; 
w_principal = [wx, wy, wz];

% Vectors in inertial frame
L_inertial = zeros(size(L_principal));
w_inertial = zeros(size(w_principal));
for i = 1:length(L_principal)
    Q4 = q4(i); Q1 = q1(i); Q2 = q2(i); Q3 = q3(i);
    R = [Q4^2+Q1^2-Q2^2-Q3^2,   2*(Q1*Q2+Q3*Q4),        2*(Q1*Q3-Q2*Q4);
         2*(Q1*Q2-Q3*Q4),       Q4^2-Q1^2+Q2^2-Q3^2,    2*(Q2*Q3+Q1*Q4);
         2*(Q1*Q3+Q2*Q4),       2*(Q2*Q3-Q1*Q4),        Q4^2-Q1^2-Q2^2+Q3^2];
    L_inertial(i,:) = R' * L_principal(i,:)';
    w_inertial(i,:) = R' * w_principal(i,:)';
end

% Plot angular velocity
figure(); grid on; hold on;
plot(t, L_inertial(:,1), 'LineWidth', 1); 
plot(t, L_inertial(:,2), 'LineWidth', 1);
plot(t, L_inertial(:,3), 'LineWidth', 1);
legend('Lx', 'Ly', 'Lz');
xlabel('t'); ylabel('kg*m*m/s'); title('Angular Momentum Vector in Inertial Frame');
hold off;

% Plot herpolhode
figure(); 
subplot(1,2,1); grid on; hold on;
plot3(w_inertial(:,1), w_inertial(:,2), w_inertial(:,3));
xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal;
title('Herpolhode'); hold off;

% Verify that herpolhode is perpendicular to angular momentum vector
% Use 2T = w dot L = const (lecture 3 slide 6)
energy = (1/2) * dot(L_inertial',w_inertial');
subplot(1,2,2);
plot(t, energy, 'LineWidth', 1); grid on;
xlabel('Time (s)'); ylabel('Rotational KE');
ylim([0 800]);
title('Rotational Kinetic Energy = (1/2) * w * L)');

%% Problem 7c: orbital frame, body axes, and principal axes as f(t)
% Plot inertial axes at t0
figure(); axisScale = 1.0; offset = 0.1; hold on;
quiver3(0,0,0,axisScale,0,0, 'Color', 'g', 'LineWidth', 2); text(axisScale + offset,offset, offset, 'X', 'Color', 'g', 'FontSize', 12);
quiver3(0,0,0,0,axisScale,0, 'Color', 'g', 'LineWidth', 2); text(offset, axisScale + offset, offset, 'Y', 'Color', 'g', 'FontSize', 12);
quiver3(0,0,0,0,0,axisScale, 'Color', 'g', 'LineWidth', 2); text(offset, offset, axisScale + offset, 'Z', 'Color', 'g', 'FontSize', 12);
xlabel('X'); ylabel('Y'); zlabel('Z');

% Plot principal axes at t0
startIndex = 1000;
rotationInertialToPrincipal = quatElementsToDcm(q1(startIndex), q2(startIndex), q3(startIndex), q4(startIndex))';
principalAxes = rotationInertialToPrincipal * diag([1,1,1]) * axisScale;
xp = principalAxes(1,:); yp = principalAxes(2,:); zp = principalAxes(3,:);
quiver3(0,0,0,xp(1),xp(2),xp(3), 'Color', [0.8 0.2 0.5], 'LineWidth', 2); text(xp(1) + offset,xp(2)+ offset,xp(3)+ offset, 'X', 'Color', [0.8 0.2 0.5], 'FontSize', 12);
quiver3(0,0,0,yp(1),yp(2),yp(3), 'Color', [0.8 0.2 0.5], 'LineWidth', 2); text(yp(1)+ offset,yp(2)+ offset,yp(3)+ offset, 'Y', 'Color', [0.8 0.2 0.5], 'FontSize', 12);
quiver3(0,0,0,zp(1),zp(2),zp(3), 'Color', [0.8 0.2 0.5], 'LineWidth', 2); text(zp(1)+ offset,zp(2)+ offset,zp(3)+ offset, 'Z', 'Color', [0.8 0.2 0.5], 'FontSize', 12);

% Plot body axes at t0
rotationPrincipalToBody = rotationBodyToPrincipal';     % DCM from principal axis frame to body axis frame - constant!
bodyAxes = rotationPrincipalToBody * principalAxes * axisScale;
xb = bodyAxes(1,:); yb = bodyAxes(2,:); zb = bodyAxes(3,:);
quiver3(0,0,0,xb(1),xb(2),xb(3), 'Color', 'b', 'LineWidth', 2); text(xb(1) + offset,xb(2)+ offset,xb(3)+ offset, 'X', 'Color', 'b', 'FontSize', 12);
quiver3(0,0,0,yb(1),yb(2),yb(3), 'Color', 'b', 'LineWidth', 2); text(yb(1)+ offset,yb(2)+ offset,yb(3)+ offset, 'Y', 'Color', 'b', 'FontSize', 12);
quiver3(0,0,0,zb(1),zb(2),zb(3), 'Color', 'b', 'LineWidth', 2); text(zb(1)+ offset,zb(2)+ offset,zb(3)+ offset, 'Z', 'Color', 'b', 'FontSize', 12);

% Find endpoints of vectors throughout simulation
numSamples=50000;
principalX = zeros(numSamples,3); principalY = zeros(numSamples,3); principalZ = zeros(numSamples,3);
bodyX = zeros(numSamples,3); bodyY = zeros(numSamples,3); bodyZ = zeros(numSamples,3);
for i = 1:numSamples
    j = i + startIndex - 1;
    rotationInertialToPrincipal = quatElementsToDcm(q1(j), q2(j), q3(j), q4(j))';
    principalAxes = rotationInertialToPrincipal * diag([1,1,1]) * axisScale;
    bodyAxes = rotationPrincipalToBody * principalAxes * axisScale;
    principalX(i,:) = principalAxes(1,:); principalY(i,:) = principalAxes(2,:); principalZ(i,:) = principalAxes(3,:);
    bodyX(i,:) = bodyAxes(1,:); bodyY(i,:) = bodyAxes(2,:); bodyZ(i,:) = bodyAxes(3,:);
end

plot3(principalX, principalY, principalZ, 'Color', [0.8 0.2 0.5], 'LineWidth', 0.1);
plot3(bodyX, bodyY, bodyZ, 'Color', 'b', 'LineWidth', 0.1);
title('Inertial (Green), Body (Blue), and Principal (Pink) Axis Frames')
grid on; axis equal; hold off;


%% Useful Functions
% Kinematic EOM for DCM
function [dcmdot] = dcmEOM(R, wx, wy, wz)
    wcross = [0 -wz wy; wz 0 -wx; -wy wx 0];
    dcmdot = -wcross * R;
end

% Kinematic EOM for quaternion
function [qdot] = quaternionEOM(q, wx, wy, wz)
    omega = [0 wz -wy wx; -wz 0 wx wy; wy -wx 0 wz; -wx -wy -wz 0];
    q = q ./ norm(q);            % Normalize
    qdot = (1/2) * omega * q';   % Calculate derivative
end

% Integrate Euler and quaternion equations from arbitrary initial conditions
function [statedot] = eulerWithQuaternion(t, state)
    
    % moments of inertia about principal axes
    Ixx = 3862.4648282509110686078201979399;
    Iyy = 4187.3328357004620556836016476154;
    Izz = 4471.171637778315016475971788168;

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

% Integrate Euler and dcm equations from arbitrary initial conditions
function [statedot] = eulerWithDCM(t, state)
    
    % moments of inertia about principal axes
    Ixx = 3862.4648282509110686078201979399;
    Iyy = 4187.3328357004620556836016476154;
    Izz = 4471.171637778315016475971788168;

    % State is [R11 R12 R13 R21 R22 R23 R31 R32 R33 wx wy wz]
    R = [state(1:3)'; state(4:6)'; state(7:9)'];
    wx = state(10); wy = state(11); wz = state(12);
    statedot = zeros (12,1);

    % Calculate rate of change of attitude, expressed as quaternion
    Rdot = dcmEOM(R, wx, wy, wz);
    statedot(1:3) = Rdot(1,:);
    statedot(4:6) = Rdot(2,:);
    statedot(7:9) = Rdot(3,:);

    % Calculate angular acceleration from Euler, in rad/s/s
    statedot(10) = -(Izz - Iyy)*(wy*wz)/Ixx;
    statedot(11) = -(Ixx - Izz)*(wz*wx)/Iyy;
    statedot(12) = -(Iyy - Ixx)*(wx*wy)/Izz;
end

% ODE integration for Euler equations (no applied forces), Ix=Iy!=Iz
function [statedot] = axialSymmetricEuler(t , state)
    % moments of inertia about principal axes
    Ixx = 3000;
    Iyy = 3000;
    Izz = 5000;

    % State is [thetax thetay thetaz wx wy wz] ' in [rad] and [rad/s]
    theta = state (1:3);
    wx = state(4); wy = state(5); wz = state(6);
    statedot = zeros (6,1);
    statedot(1) = wx; statedot(2) = wy; statedot(3) = wz;     % thetadot = angular velocity
    statedot(4) = -(Izz - Iyy)*(wy*wz)/Ixx; % From Euler equations with no applied moment
    statedot(5) = -(Ixx - Izz)*(wz*wx)/Iyy;
    statedot(6) = -(Iyy - Ixx)*(wx*wy)/Izz;
end

% Converts quaternion to DCM
function [dcm] = quatElementsToDcm(Q1,Q2,Q3,Q4) 
    dcm = [Q4^2+Q1^2-Q2^2-Q3^2,   2*(Q1*Q2+Q3*Q4),        2*(Q1*Q3-Q2*Q4);
         2*(Q1*Q2-Q3*Q4),       Q4^2-Q1^2+Q2^2-Q3^2,    2*(Q2*Q3+Q1*Q4);
         2*(Q1*Q3+Q2*Q4),       2*(Q2*Q3-Q1*Q4),        Q4^2-Q1^2-Q2^2+Q3^2];
end
