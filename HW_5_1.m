%% AA 279C HW 5-1
% Natasha Evans

%% Load files
load simulationSettings
load orbitGeometry
load inertiaTensors
load massAndGeometryProperties

%% Test stability given MOI and a selected axis/axes
plotGravityGradientStability(Ixx, Iyy, Izz)

% Define initial conditions
initialRotation = deg2rad(25);
[r0, v0] = oe2mci(a,e,i,raan,argp,trueAnomaly0);        % initial position w.r.t. inertial frame
mci2rtn_initial = inertialToRTN(r0',v0');
[q10, q20, q30, q40] = dcm2Quaternion(mci2rtn_initial); % initial attitude w.r.t. inertial frame
q0 = [q10, q20, q30, q40];

% Test stability
for ax = 3
    % Run unperturbed simulation
    w0 = zeros(1,3);
    w0(ax) = initialRotation;
    y0 = [r0', v0', q0, w0];               % initial state vector
    [t, y_out] = ode113(@integrateOrbitAttitudeGG, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz, mu_mars);
    pos = y_out(:,1:3); vel = y_out(:,4:6); q = y_out(:,7:10); w = y_out(:,11:13);
    theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
    for k = 1:size(t)
        dcm = quaternion2DCM(q(k,1),q(k,2),q(k,3),q(k,4));
        [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
    end
    
    % Plot unperturbed results
    axes = ['X', 'Y', 'Z'];
    pltTitle = sprintf('Equilibrium With Rotation About %s Axis', axes(ax));
    figure(); sgtitle(pltTitle);
    subplot(2,2,1); hold on;
    plot(t, w(:,1)); plot(t,w(:,2)); plot(t,w(:,3)); grid on;
    legend('X', 'Y', 'Z');
    title('Angular Velocity (Unperturbed)');
    xlabel('Time (s)'); ylabel('Rad/s'); hold off;
    
    subplot(2,2,2); hold on;
    plot(t, theta); plot(t, phi); plot(t,psi); grid on;
    legend('Roll', 'Yaw', 'Pitch');
    title('Euler Angles (Unperturbed)');
    xlabel('Time (s)'); ylabel('Rad/s'); hold off;
    
    % Perturb equilibrium
    w0 = w0 + 0.01*initialRotation;
    y0 = [r0', v0', q0, w0];               % initial state vector
    [t, y_out] = ode113(@integrateOrbitAttitudeGG, [0:tstep:t_period]' , y0, options, Ixx, Iyy, Izz, mu_mars);
    pos = y_out(:,1:3); vel = y_out(:,4:6); q = y_out(:,7:10); w = y_out(:,11:13);
    theta = zeros(size(t)); phi = zeros(size(t)); psi = zeros(size(t));
    for k = 1:size(t)
        dcm = quaternion2DCM(q(k,1),q(k,2),q(k,3),q(k,4));
        [theta(k), phi(k), psi(k)] = dcm2Euler312(dcm);
    end
    
    % Plot perturbed results
    subplot(2,2,3); hold on;
    plot(t, w(:,1)); plot(t,w(:,2)); plot(t,w(:,3)); grid on;
    legend('X', 'Y', 'Z');
    title('Angular Velocity (Perturbed)');
    xlabel('Time (s)'); ylabel('Rad/s'); hold off;
    
    subplot(2,2,4); hold on;
    plot(t, theta); plot(t, phi); plot(t,psi); grid on;
    legend('Roll', 'Yaw', 'Pitch');
    title('Euler Angles (Perturbed)');
    xlabel('Time (s)'); ylabel('Rad/s'); hold off;
end


%% Useful Functions
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


% Convert DCM to quaternion
function [q1,q2,q3,q4] = dcm2Quaternion(A)
    q4 = (1/2) * (1 + A(1,1) + A(2,2) + A(3,3))^(1/2);
    q1 = (A(2,3) - A(3,2)) / (4*q4);
    q2 = (A(3,1) - A(1,3)) / (4*q4);
    q3 = (A(1,2) - A(2,1)) / (4*q4);
end


% Integrate orbit and attitude with gravity gradient
function [statedot] = integrateOrbitAttitudeGG(t, state, Ixx, Iyy, Izz, mu_mars)

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
    torque = gravityGradientTorque(r(1),r(2),r(3),v(1),v(2),v(3),q1,q2,q3,q4,mu_mars,Ixx,Iyy,Izz);
    statedot(11) = (torque(1)-(Izz - Iyy)*(wy*wz))/Ixx;
    statedot(12) = (torque(2)-(Ixx - Izz)*(wz*wx))/Iyy;
    statedot(13) = (torque(3)-(Iyy - Ixx)*(wx*wy))/Izz;

end


% Gravity gradient torque in princpal axes (lec7, slide 9) 
function torque = gravityGradientTorque(rx,ry,rz,vx,vy,vz,q1,q2,q3,q4,mu,Ix,Iy,Iz)
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

% Convert DCM to Euler (312)
function [theta, phi, psi] = dcm2Euler312(A)
    theta = asin(A(3,2));
    phi = atan2(A(1,2), A(2,2));
    psi = atan2(A(3,1),A(3,3));
end

% Convert Quaternion to DCM - Q4 is scalar component
function [dcm] = quaternion2DCM(Q1,Q2,Q3,Q4) 
    dcm = [Q4^2+Q1^2-Q2^2-Q3^2,   2*(Q1*Q2+Q3*Q4),        2*(Q1*Q3-Q2*Q4);
         2*(Q1*Q2-Q3*Q4),       Q4^2-Q1^2+Q2^2-Q3^2,    2*(Q2*Q3+Q1*Q4);
         2*(Q1*Q3+Q2*Q4),       2*(Q2*Q3-Q1*Q4),        Q4^2-Q1^2-Q2^2+Q3^2];
end

% Kinematic EOM for quaternion
function [qdot] = quaternionEOM(q, wx, wy, wz)
    omega = [0 wz -wy wx; -wz 0 wx wy; wy -wx 0 wz; -wx -wy -wz 0];
    q = q ./ norm(q);            % Normalize
    qdot = (1/2) * omega * q';   % Calculate derivative
end

% Calculate rotation matrix from inertial to RTN
function [dcm] = inertialToRTN(rvec,vvec)
    
    Rhat = rvec ./ norm(rvec);
    Nhat = cross(rvec,vvec) ./ norm(cross(rvec,vvec));
    That = cross(Nhat, Rhat);
    
    % Initial conditions
    dcm = [Rhat; That; Nhat];
end


% Plot stability re: gravity gradient
function [kn_xc, kt_sc, kr_sc] = plotGravityGradientStability(Ixx, Iyy, Izz) 
    % Compute coefficients
    kn_sc = (Iyy-Ixx)/Izz;
    kt_sc = (Izz-Ixx)/Iyy;
    kr_sc = (Izz-Iyy)/Ixx;
    
    % Plot regions of stable and unstable motion
    % Criteria = kr values
    kt = -1:0.01:1; kr = -1:0.01:1;
    figure(); hold on; grid on;
    xlim([-1, 1]); ylim([-1,1]);
    xlabel('K_t'); ylabel('K_r')
    criteria1 = kt;
    criteria2 = (3^(1/2).*(1 - kt).^(1/2) - 2).^2./kt;
    
    % Unstable pitch
    xPatch = [kt, fliplr(kt)];
    unstablePitch = [kr, max(ylim)*ones(size(kr))];
    patch(xPatch, unstablePitch, 'y', 'EdgeColor', 'black');
    
    % Unstable yaw/roll
    xvals = [0:0.01:1];
    xPatch = [xvals, fliplr(xvals)];
    yPatch = [-1*ones(size(xvals)), zeros(size(xvals))];
    patch(xPatch, yPatch, 'b', 'EdgeColor', 'black');
    
    xvals = [-1:0.01:0];
    xPatch = [xvals, fliplr(xvals)];
    yPatch = [zeros(size(xvals)), ones(size(xvals))];
    patch(xPatch, yPatch, 'b', 'EdgeColor', 'black');
    
    criteria2Region = criteria2(kt < 0);
    xPatch = [fliplr(kt(kt < 0)), kt(kt < 0)];
    yPatch = [-1*ones(size(criteria2Region)), criteria2Region];
    patch(xPatch, yPatch, 'b', 'EdgeColor', 'black');
    
    % Unstable yaw, pitch, and roll
    xvals = [-1:0.01:0];
    xPatch = [xvals, fliplr(xvals)];
    yPatch = [zeros(size(xvals)), ones(size(xvals))];
    patch(xPatch, yPatch, 'g', 'EdgeColor', 'black');
    
    mask = criteria2Region > kt(kt < 0);
    xPatch = [kt(mask), fliplr(kt(mask))];
    yPatch = [kr(mask), fliplr(criteria2Region(mask))];
    patch(xPatch, yPatch, 'g', 'EdgeColor', 'black');
    
    % Plot satellite location
    scatter(kt_sc, kr_sc, 'r*');
    
    % Finalize plot
    legend('Unstable Pitch', 'Unstable Yaw & Roll', '', '', 'Unstable Yaw, Pitch & Roll', '', 'Spacecraft');
    title('Gravity Gradient Stability'); hold off;

end
