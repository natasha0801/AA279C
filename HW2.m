%% AA 279C HW 2
load massAndGeometryProperties

%% Problem 1: Numerical Simulation of Orbit
% Orbit elements & Mars parameters
r_mars = 3389.5;        % Mars radius (km)
mu_mars = 4.28284e4;    % gravitational parameter (km^3/s^2)
i = 75;                 % inclination (deg)
h_a = 6200;             % apoapsis altitude (km)
h_p = 150;              % periapsis altitude (km)
r_a = h_a + r_mars;     % apoapsis radius (km)
r_p = h_p + r_mars;     % periapsis radius (km)
a = 0.5 * (r_a + r_p);          % semimajor axis (km)
e = (r_a - r_p)/(r_a + r_p);    % eccentricity
trueAnomaly0 = 0;       % initial true anomaly (deg)
argp = 0;               % initial argument of periapsis, assumed (deg)
raan = 0;               % initial raan, assumed (deg)
[r0, v0] = oe2mci(a,e,i,raan,argp,trueAnomaly0);
y0 = [r0;v0];

% Run simulation
options = odeset ( 'RelTol' , 1e-3 , 'AbsTol' , 1e-6);
t_period = 2 * pi * sqrt(a*a*a/mu_mars);
tstep = 0.01*t_period;
[t, y_out] = ode113(@orbitingbody , [0:tstep:t_period]' , y0 , options );
[xm, ym, zm] = ellipsoid(0,0,0,r_mars, r_mars, r_mars,30);        % surface coordinates (for plot)
figure(); 
s = surface(xm,ym,zm,'FaceColor', [0.8510 0.5050 0.4030], 'EdgeColor','black'); 
s.FaceAlpha = 0.3;
hold on;
plot3(y_out(:,1), y_out(:,2), y_out(:,3), 'LineWidth', 2); axis equal;
grid on; xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)'); axis equal;
title('Numerically Integrated MAVEN Orbit (No Perturbations)')

%% Problem 2: Eigenvector/Eigenvalue Problem
[rotationMatrix, principalInertiaTensor] = eig(inertiaTensor);
[maxval, idx] = max(principalInertiaTensor);

%% Problem 3: Plot Body Axes and Principal Axes
figure(); hold on; axis equal; grid on;
title('MAVEN with Body and Principal Axes');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
% Plot spacecraft surfaces
for i=1:length(surfaces)
    surface = surfaces{i};
    x = surface(:,1); y = surface(:,2); z = surface(:,3);
    h = fill3(x,y,z,color);
    h.FaceAlpha = 0.25;
end
clear surface

% Plot body axes
axisScale = 3;
offset = 0.2;
quiver3(0,0,0,axisScale,0,0, 'Color', 'r', 'LineWidth', 1); text(axisScale + offset,offset, offset, 'Xb', 'Color', 'r', 'FontSize', 12);
quiver3(0,0,0,0,axisScale,0, 'Color', 'r', 'LineWidth', 1); text(offset, axisScale + offset, offset, 'Yb', 'Color', 'r', 'FontSize', 12);
quiver3(0,0,0,0,0,axisScale, 'Color', 'r', 'LineWidth', 1); text(offset, offset, axisScale + offset, 'Zb', 'Color', 'r', 'FontSize', 12);

% Plot principal axes
principalAxes = rotationMatrix * diag([1,1,1]) * axisScale;
xdir = principalAxes(1,:); ydir = principalAxes(2,:); zdir = principalAxes(3,:);
quiver3(0,0,0,xdir(1),xdir(2),xdir(3), 'Color', 'magenta', 'LineWidth', 1); text(xdir(1) + offset,xdir(2)+ offset,xdir(3)+ offset, 'Xp', 'Color', 'magenta', 'FontSize', 12);
quiver3(0,0,0,ydir(1),ydir(2),ydir(3), 'Color', 'magenta', 'LineWidth', 1); text(ydir(1)+ offset,ydir(2)+ offset,ydir(3)+ offset, 'Yp', 'Color', 'magenta', 'FontSize', 12);
quiver3(0,0,0,zdir(1),zdir(2),zdir(3), 'Color', 'magenta', 'LineWidth', 1); text(zdir(1)+ offset,zdir(2)+ offset,zdir(3)+ offset, 'Zp', 'Color', 'magenta', 'FontSize', 12);
hold off;

%% Problem 4: Euler equations in principal axes (no external torques)
Ix = principalInertiaTensor(1,1); Iy = principalInertiaTensor(2,2); Iz = principalInertiaTensor(3,3);    % principal MOI
% See eulerPrincipalAxes function.

%% Problem 5: Integrate Euler equations
% Run simulation
wx0 = 45; wy0 = -20; wz0 = -15;  % initial rotation rates (deg/s)
y0 = zeros(6,1);        % zero initial angle
y0(4) = deg2rad(wx0);   % initial rotation rate about x axis (rad/s)
y0(5) = deg2rad(wy0);   % initial rotation rate about y axis (rad/s)
y0(6) = deg2rad(wz0);   % initial rotation rate about z axis (rad/s)

options = odeset ( 'RelTol' , 1e-3 , 'AbsTol' , 1e-6);
t_period = 20*(360/wx0);
tstep = 0.05;
[t, y_out] = ode113(@eulerPrincipalAxes, [0:tstep:t_period]' , y0, options);
thetax = y_out(:,1); thetay = y_out(:,2); thetaz = y_out(:,3);      % angle (rad)
wx = y_out(:,4); wy = y_out(:,5); wz = y_out(:,6);                  % rotation rate (rad/s)

%% Problem 6: Plot rotational kinetic energy and momentum ellipsoids
% Momemtum ellipsoid
angularVelocityVector = [wx,wy,wz]';
L = principalInertiaTensor * angularVelocityVector;
Lmag = vecnorm(L); Lmag = Lmag(1);          % angular momentum magnitude (constant)
a = Lmag/Ix; b = Lmag/Iy; c = Lmag/Iz;      % Axes of momentum ellipsoid
figure(); 
[xme, yme, zme] = ellipsoid(0,0,0,a,b,c,45);
momentumEllipsoid = surf(xme,yme,zme);
momentumEllipsoid.FaceColor = 'r';
momentumEllipsoid.FaceAlpha = 0.15;

pltTitle = sprintf('Momentum & Energy Ellipsoid\n w0 = [%.1f, %.1f, %.1f] deg/s', wx0,wy0,wz0);
axis equal; title(pltTitle); hold on;
xlabel('X [rad/s]'); ylabel('Y [rad/s]'); zlabel('Z [rad/s]');

% Energy ellipsoid
T = 0.5 * dot(angularVelocityVector, L);  % rotational kinetic energy - should be constant
T = T(1);
a = sqrt(2*T/Ix); b = sqrt(2*T/Iy); c = sqrt(2*T/Iz);
[xee, yee, zee] = ellipsoid(0,0,0,a,b,c,45);
energyEllipsoid = surf(xee, yee, zee);
energyEllipsoid.FaceColor = 'b';
energyEllipsoid.FaceAlpha = 0.15;
legend('Momentum Ellipsoid', 'Energy Ellipsoid');

%% Problem 7 - Plot Polhode
% Plot polhode - evolution of angular velocity vector over time
if (norm(cross(L(:,1), [wx0 wy0 wz0]')) == 0)   % if w || L, polhode is a point
    scatter3(wx, wy, wz, 'Filled', 'go')
else
    plot3(wx,wy,wz, 'Color', 'b', 'LineWidth',3);
end
grid on;
legend('Momentum Ellipsoid', 'Energy Ellipsoid', 'Polhode');

%% Problem 8 - Plot polhode in 3x 2D planes identified by principal axes
theta = -pi:1e-5:pi;
r_yz = @(theta) sqrt((Lmag*Lmag - 2*T*Ix)./((Iy-Ix)*Iy.*cos(theta).^2 + (Iz-Ix)*Iz.*sin(theta).^2));        % as seen on Y-Z plane
r_xz = @(theta) sqrt(abs((Lmag*Lmag - 2*T*Iy)./((Ix-Iy)*Ix.*cos(theta).^2 + (Iz-Iy)*Iz.*sin(theta).^2)));   % as seen on X-Z plane
r_xy = @(theta) sqrt((Lmag*Lmag - 2*T*Iz)./((Ix-Iz)*Ix.*cos(theta).^2 + (Iy-Iz)*Iy.*sin(theta).^2));        % as seen on X-Y plane
figure();

if (norm(cross(L(:,1), [wx0 wy0 wz0]')) == 0)   % if w || L, polhode is a point
    subplot(1,3,1);
    plot(r_yz(theta).*cos(theta), r_yz(theta).*sin(theta), 'LineWidth', 1); 
    grid on; axis equal; hold on;
    scatter(wy, wz, 'filled', 'o');
    xlabel('Y [rad/s]'); ylabel('Z [rad/s]'); title('Y-Z');
    legend('Theoretical', 'Actual');
    
    subplot(1,3,2);
    plot(r_xz(theta).*cos(theta), r_xz(theta).*sin(theta), 'LineWidth', 1);
    grid on; xlim([min(wx)-1 max(wx)+1]); ylim([min(wz)-1 max(wz)+1]); hold on;
    scatter(wx, wz, 'filled', 'o');
    xlabel('X [rad/s]'); ylabel('Z [rad/s]'); title('X-Z');
    legend('Theoretical', 'Actual');
    
    subplot(1,3,3);
    plot(r_xy(theta).*cos(theta), r_xy(theta).*sin(theta), 'LineWidth', 1);  
    grid on; axis equal; hold on;
    scatter(wx, wy, 'filled', 'o');
    xlabel('X [rad/s]'); ylabel('Y [rad/s]'); title('X-Y');
    legend('Theoretical', 'Actual');

else
    subplot(1,3,1);
    plot(r_yz(theta).*cos(theta), r_yz(theta).*sin(theta), 'LineWidth', 1); 
    grid on; axis equal; hold on;
    plot(wy, wz, 'LineStyle', ':', 'LineWidth', 1);
    xlabel('Y [rad/s]'); ylabel('Z [rad/s]'); title('Y-Z');
    legend('Theoretical', 'Actual');
    
    subplot(1,3,2);
    plot(r_xz(theta).*cos(theta), r_xz(theta).*sin(theta), 'LineWidth', 1);
    grid on; xlim([min(wx) max(wx)]); ylim([min(wz) max(wz)]); hold on;
    plot(wx, wz, 'LineStyle', ':', 'LineWidth', 1);
    xlabel('X [rad/s]'); ylabel('Z [rad/s]'); title('X-Z');
    legend('Theoretical', 'Actual');
    
    subplot(1,3,3);
    plot(r_xy(theta).*cos(theta), r_xy(theta).*sin(theta), 'LineWidth', 1);  
    grid on; axis equal; hold on;
    plot(wx, wy, 'LineStyle', ':', 'LineWidth', 1);
    xlabel('X [rad/s]'); ylabel('Y [rad/s]'); title('X-Y');
    legend('Theoretical', 'Actual');
end

pltTitle = sprintf('Polhode Projection\n w0 = [%.1f, %.1f, %.1f] deg/s', wx0, wy0, wz0);
sgtitle(pltTitle);

%% Useful Functions
% ODE integration function for Euler equations (no applied forces)
function [statedot] = eulerPrincipalAxes(t , state)
    % moments of inertia about principal axes
    Ixx = 3862.4648282509110686078201979399;
    Iyy = 4187.3328357004620556836016476154;
    Izz = 4471.171637778315016475971788168;

    % State is [thetax thetay thetaz wx wy wz] ' in [rad] and [rad/s]
    theta = state (1:3);
    wx = state(4); wy = state(5); wz = state(6);
    statedot = zeros (6,1);
    statedot(1) = wx; statedot(2) = wy; statedot(3) = wz;     % thetadot = angular velocity
    statedot(4) = -(Izz - Iyy)*(wy*wz)/Ixx; % From Euler equations with no applied moment
    statedot(5) = -(Ixx - Izz)*(wz*wx)/Iyy;
    statedot(6) = -(Iyy - Ixx)*(wx*wy)/Izz;

    % Note that angle "accumulates" over time in this representation
    % Will later switch to quaternions
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

% ODE integration function for orbit propagation (2-body problem)
function [statedot] = orbitingbody (t , state)
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