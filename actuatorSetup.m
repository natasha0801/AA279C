%% Geometry setup for thrusters
load inertiaTensors rotationBodyToPrincipal
load massAndGeometryProperties.mat

% Center of mass in body frame
positions_body = [hgaXYZ; boxXYZ; spInnerLeftXYZ; spInnerRightXYZ; spOuterLeftXYZ; spOuterRightXYZ; appArmXYZ; appPlatformXYZ];
masses = [hgaMass, boxMass, solarPanelMass, solarPanelMass, solarPanelMass, solarPanelMass, appArmMass, appPlatformMass]';
cm_body = sum(masses.*positions_body) / sum(masses);
cm_principal = rotationBodyToPrincipal * cm_body';

% Position of thrusters in principal frame rel. to cm
r_t1 = [boxHeight-cm_principal(1), boxWidth/2-cm_principal(2), boxWidth/2-cm_principal(3)];
r_t2 = [boxHeight-cm_principal(1), -boxWidth/2-cm_principal(2), boxWidth/2-cm_principal(3)];
r_t3 = [boxHeight-cm_principal(1), -boxWidth/2-cm_principal(2), -boxWidth/2-cm_principal(3)];
r_t4 = [boxHeight-cm_principal(1), boxWidth/2-cm_principal(2), -boxWidth/2-cm_principal(3)];

% Direction vector of thrusters in principal axis frame
e_t1 = r_t1 ./ norm(r_t1);
e_t2 = r_t2 ./ norm(r_t2);
e_t3 = r_t3 ./ norm(r_t3);
e_t4 = r_t4 ./ norm(r_t4);
e_t1(2) = 0; e_t2(2) = 0; e_t3(2) = 0; e_t4(2) = 0;

% Visualization of thruster location
r_t1b = rotationBodyToPrincipal' * (r_t1' + cm_principal) ;
r_t2b = rotationBodyToPrincipal' * (r_t2' + cm_principal) ;
r_t3b = rotationBodyToPrincipal' * (r_t3' + cm_principal) ;
r_t4b = rotationBodyToPrincipal' * (r_t4' + cm_principal) ;
e_t1b = rotationBodyToPrincipal' * e_t1' ;
e_t2b = rotationBodyToPrincipal' * e_t2' ;
e_t3b = rotationBodyToPrincipal' * e_t3' ;
e_t4b = rotationBodyToPrincipal' * e_t4' ;

uiopen('C:\Users\natas\Documents\Classes\AA 279C\MATLAB\scBodyFramePlot.fig',1)
hold on
scatter3(r_t1b(1), r_t1b(2), r_t1b(3), 'm', 'filled');
scatter3(r_t2b(1), r_t2b(2), r_t2b(3), 'm', 'filled');
scatter3(r_t3b(1), r_t3b(2), r_t3b(3), 'm', 'filled');
scatter3(r_t4b(1), r_t4b(2), r_t4b(3), 'm', 'filled');
quiver3(r_t1b(1), r_t1b(2), r_t1b(3), e_t1b(1), e_t1b(2), e_t1b(3), 'm', 'LineWidth', 1, 'MaxHeadSize', 4);
quiver3(r_t2b(1), r_t2b(2), r_t2b(3), e_t2b(1), e_t2b(2), e_t2b(3), 'm', 'LineWidth', 1, 'MaxHeadSize', 4);
quiver3(r_t3b(1), r_t3b(2), r_t3b(3), e_t3b(1), e_t3b(2), e_t3b(3), 'm', 'LineWidth', 1, 'MaxHeadSize', 4);
quiver3(r_t4b(1), r_t4b(2), r_t4b(3), e_t4b(1), e_t4b(2), e_t4b(3), 'm', 'LineWidth', 1, 'MaxHeadSize', 4);

% Thruster position rel. to center of mass in body frame
r_t1 = r_t1 - cm_principal';
r_t2 = r_t2 - cm_principal';
r_t3 = r_t3 - cm_principal';
r_t4 = r_t4 - cm_principal';

% Save thruster geometry in principal axis frame
r_thrusters = [r_t1; r_t2; r_t3; r_t4];
e_thrusters = [e_t1; e_t2; e_t3; e_t4];
A_thrusters = [cross(r_t1, e_t1)', cross(r_t2, e_t2)', cross(r_t3, e_t3)', cross(r_t4, e_t4)'];
save thrusterGeometry r_thrusters e_thrusters A_thrusters

%% Save 3-axis RW properties
% https://satsearch.co/products/oce-technology-rw1000-reaction-wheel
A_RW = [1 0 0; 0 1 0; 0 0 1];            % mounting matrix
r_RW = 0.337/2;           % [m]
m_RW = 10;             % [kg]
I_RW = 0.5*m_RW*r_RW^2; % [kg * m^2]
maxRWTorque = 1;    % [N*m]
maxRWSpeed = 1200*(2*pi)/60;    % [rad/s]
save rwProperties A_RW r_RW m_RW I_RW maxRWTorque maxRWSpeed

%% Test Actuators
dt = 0.1;
t = 0:dt:300;
Mc = 0.001 * [10, -1, 5]' .* sin(t/20) .* t;
rwSpeeds = zeros(3,length(t));
rwTorques = zeros(3,length(t));
thrusterTorques = zeros(3,length(t));
thrust = zeros(4,length(t));
for j = 1:length(t)-1
    [rwTorques(:,j+1), rwSpeeds(:,j+1)] = ReactionWheels(Mc(:,j), rwSpeeds(:,j), dt);
    [thrusterTorques(:,j+1), thrust(:,j+1)] = Thrusters(Mc(:,j));
end

figure(); sgtitle('RW Performance');
for ax = 1:3
    subplot(4,1,ax); hold on; grid on;
    plot(t, rwTorques(ax,:)); plot(t, Mc(ax,:));
    legend('Actual', 'Command');
    xlabel('Time (s)'); ylabel('Torque (Nm)');
    stitle = sprintf("Axis %i", ax);
    title(stitle);
end
subplot(4,1,4); title('RW Speeds'); hold on;
plot(t, rwSpeeds(1,:)); plot(t,rwSpeeds(2,:)); plot(t, rwSpeeds(3,:));
xlabel('Time (s)'); ylabel('Speed (rad/s)');
legend('RX', 'RY', 'RZ');
grid on;

figure(); sgtitle('Thruster Performance');
for ax = 1:3
    subplot(4,1,ax); hold on; grid on;
    plot(t, thrusterTorques(ax,:)); plot(t, Mc(ax,:));
    legend('Actual', 'Command');
    xlabel('Time (s)'); ylabel('Torque (Nm)');
    stitle = sprintf("Axis %i", ax);
    title(stitle);
end
subplot(4,1,4); title('Thrust'); hold on;
plot(t, thrust(1,:)); plot(t,thrust(2,:)); plot(t, thrust(3,:)); plot(t, thrust(4,:));
xlabel('Time (s)'); ylabel('Thrust (N)');
legend('T1', 'T2', 'T3', 'T4')
grid on;

%% Actuator gains from Euler equations
syms Ix Iy Iz Iw w1 w2 w3 wd1 wd2 wd3 wx wy wz wdx wdy wdz n cx cy cz
I_sat = [Ix; Iy; Iz];
w_sat = [wx; wy; wz];
wd_sat = [wdx; wdy; wdz];
w_rw = [w1; w2; w3];
wd_rw = [wd1; wd2; wd3];
load rwProperties A_RW
rhs = 3*n*n*[(Iz-Iy)*cy*cz; (Ix-Iz)*cz*cx; (Iy-Ix)*cx*cy];
eul = rhs == I_sat.*wd_sat + cross(w_sat, I_sat.*w_sat) + A_RW*Iw*wd_rw + cross(w_sat, A_RW*Iw*w_rw);
disp("Euler Equations (Non-Linearized): ")
disp(eul)

syms ax adx addx ay ady addy az adz addz
eul = subs(eul, wdx, addx - n*ady);
eul = subs(eul, wdy, addy - n*adx);
eul = subs(eul, wdz, addz);
eul = subs(eul, wx, adx - n*ay);
eul = subs(eul, wy, ady - n*ax);
eul = subs(eul, wx, n);
eul = subs(eul, cy*cz,0);
eul = subs(eul, cz*cx,ay);
eul = subs(eul, cx*cy, -az);

disp("Euler Equations (Linearized): ")
disp(eul);

syms Mcx Mcy Mcz
eul = subs(eul, Iw*wd1, Mcx);
eul = subs(eul, Iw*wd2, Mcy);
eul = subs(eul, Iw*wd3, Mcz);

disp("Euler Equations (Substituting Actuator Moment): ")
disp(eul)

syms kpx kpy kpz kdx kdy kdz
eul = subs(eul, Mcx, kpx*ax + kdx*adx);
eul = subs(eul, Mcy, kpy*ay + kdy*ady);
eul = subs(eul, Mcz, kpz*az + kdz*adz);

disp("Linear control law: ")
pretty(simplify(eul))
%% Actuator Functions
function [torque, thrust] = Thrusters(Mc)
    load thrusterGeometry
    maxThrust = 1; minThrust = 0;       % N
    thrust = pinv(A_thrusters) * Mc;

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
