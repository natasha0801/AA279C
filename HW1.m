%% AA279C HW 1
% Natasha Evans
% SUID: 06743290

%% Spacecraft dimensions and masses
% Central structure dimensions
boxHeight = 2;     % [m]
boxWidth = 2.29;   % [m]
boxLength = 2.29;  % [m]
boxVolume = boxHeight * boxWidth * boxLength; % [m^3]

% Solar panel dimensions
solarPanelWidth = 2.08;         % [m]
solarPanelLength = 1.68;        % [m]
solarPanelThickness = 0.0254;   % [m]
solarPanelVolume = solarPanelWidth * solarPanelLength * solarPanelThickness; % [m^3]
solarPanelAngle = 20;   % [deg]
RtoOuterPanel = [1 0 0; 0 cosd(solarPanelAngle) -sind(solarPanelAngle); 0 sind(solarPanelAngle) cosd(solarPanelAngle)];

% HGA dimensions
hgaDiameter = 2.1;          % [m]
hgaRadius = hgaDiameter/2;  % [m]
hgaHeight = 0.871;          % [m]
hgaVolume = pi * (hgaDiameter/2)^2 * hgaHeight / 3;  % [m^3]

% Payload Platform dimensions
appArmLength = 0.95;        % [m]
appArmDiameter = 0.05;      % [m]
appPlatformSide = 0.72;     % [m]
appArmVolume = (appArmLength*appArmDiameter*appArmDiameter*pi/4); % [m^3]
appPlatformVolume = appPlatformSide^3;  % [m^3]

% Masses
primaryStructureMass = 125; % [kg]
propMass = 1640;            % [kg]
launchMass = 2454;          % [kg]
boxMass = primaryStructureMass + propMass;  % [kg]
boxDensity = boxMass / boxVolume;           % [kg/m^3]
remainingMass = launchMass - boxMass;       % [kg]
remainingVolume = hgaVolume + 4*solarPanelVolume + appArmVolume + appPlatformVolume; % [m^3]
remainingDensity = remainingMass / remainingVolume;   % [kg/m^3]
hgaMass = remainingDensity * hgaVolume;               % [kg]
solarPanelMass = remainingDensity * solarPanelVolume; % [kg]
appArmMass = remainingDensity * appArmVolume;         % [kg]
appPlatformMass = remainingDensity*appPlatformVolume; % [kg]

%% Inertia tensor calculations
% Compute MOI of each part about its own axes [kg * m^2]
I_app_arm_xx = (1/12) * appArmMass * (3*(appArmDiameter/2)^2 + appArmLength^2);
I_app_arm_yy = 0.5 * appArmMass * (appArmDiameter/2)^2;
I_app_arm_zz = I_app_arm_xx;
I_app_platform_xx = appPlatformMass * appPlatformSide^2 / 6;
I_app_platform_yy = I_app_platform_xx;
I_app_platform_zz = I_app_platform_xx;
I_sp_inner_xx = (1/12) * solarPanelMass * (solarPanelLength^2 + solarPanelThickness^2);
I_sp_inner_yy = (1/12) * solarPanelMass * (solarPanelWidth^2 + solarPanelThickness^2);
I_sp_inner_zz = (1/12) * solarPanelMass * (solarPanelLength^2 + solarPanelWidth^2);
I_hga_xx = (3/20) * hgaMass * (hgaRadius^2 + 4*hgaHeight^2);
I_hga_xx = I_hga_xx - hgaMass*(0.75*hgaHeight)^2;   % about CM
I_hga_yy = I_hga_xx;
I_hga_zz = (3/10) * hgaMass * hgaRadius^2;
I_center_xx = (1/12) * boxMass * (boxWidth^2 + boxHeight^2);
I_center_yy = (1/12) * boxMass * (boxLength^2 + boxHeight^2);
I_center_zz = (1/12) * boxMass * (boxWidth^2 + boxLength^2);
I_sp_outer = RtoOuterPanel * diag([I_sp_inner_xx; I_sp_inner_yy; I_sp_inner_zz]) * RtoOuterPanel';
I_sp_outer_xx = I_sp_outer(1,1);
I_sp_outer_yy = I_sp_outer(2,2);
I_sp_outer_zz = I_sp_outer(3,3);

% Coordinates of CM of each component relative to body frame [m]
hgaXYZ = [0,0,0.25*hgaHeight];
boxXYZ = [0,0,-boxHeight/2];
spInnerLeftXYZ = [0,-(boxWidth/2 + solarPanelLength/2), -solarPanelThickness/2];
spInnerRightXYZ = [0, (boxWidth/2 + solarPanelLength/2), -solarPanelThickness/2];
spOuterLeftXYZ = [0,-(boxWidth/2 + solarPanelLength+solarPanelLength*cosd(20)/2),solarPanelLength*sind(20)/2];
spOuterRightXYZ = [0,(boxWidth/2 + solarPanelLength+solarPanelLength*cosd(20)/2),solarPanelLength*sind(20)/2];
appArmXYZ = [(boxWidth/2 + appArmLength/2), 0, -appArmDiameter/2];
appPlatformXYZ = [(boxWidth/2 + appArmLength + appPlatformSide/2), 0, -appArmDiameter/2];

% total MOI from parallel axis theorem[kg * m^2]
Ix = 2*(I_sp_inner_xx + solarPanelMass*(spInnerLeftXYZ(2)^2 + spInnerLeftXYZ(3)^2)) ...
    + 2*(I_sp_outer_xx + solarPanelMass*(spOuterLeftXYZ(2)^2 + spOuterLeftXYZ(3)^2))...
    + I_app_arm_xx + appArmMass*(appArmXYZ(2)^2 + appArmXYZ(3)^2)...
    + I_app_platform_xx + appPlatformMass*(appPlatformXYZ(2)^2 + appPlatformXYZ(3)^2)...
    + I_center_xx + boxMass*(boxXYZ(2)^2 + boxXYZ(3))^2 ...
    + I_hga_xx + hgaMass*(hgaXYZ(2)^2 + hgaXYZ(3)^2);
Iy = 2*(I_sp_inner_yy + solarPanelMass*(spInnerLeftXYZ(1)^2 + spInnerLeftXYZ(3)^2)) ...
    + 2*(I_sp_outer_yy + solarPanelMass*(spOuterLeftXYZ(1)^2 + spOuterLeftXYZ(3)^2))...
    + I_app_arm_yy + appArmMass*(appArmXYZ(1)^2 + appArmXYZ(3)^2)...
    + I_app_platform_yy + appPlatformMass*(appPlatformXYZ(1)^2 + appPlatformXYZ(3)^2)...
    + I_center_yy + boxMass*(boxXYZ(1)^2 + boxXYZ(3))^2 ...
    + I_hga_yy + hgaMass*(hgaXYZ(1)^2 + hgaXYZ(3)^2);
Iz = 2*(I_sp_inner_zz + solarPanelMass*(spInnerLeftXYZ(1)^2 + spInnerLeftXYZ(2)^2)) ...
    + 2*(I_sp_outer_zz + solarPanelMass*(spOuterLeftXYZ(1)^2 + spOuterLeftXYZ(2)^2))...
    + I_app_arm_zz + appArmMass*(appArmXYZ(2)^2 + appArmXYZ(1)^2)...
    + I_app_platform_zz + appPlatformMass*(appPlatformXYZ(2)^2 + appPlatformXYZ(1)^2)...
    + I_center_zz + boxMass*(boxXYZ(1)^2 + boxXYZ(2))^2 ...
    + I_hga_zz + hgaMass*(hgaXYZ(1)^2 + hgaXYZ(2)^2);

% total products of inertia - all components are symmetrical except a.p.p.
Ixy = appArmMass*appArmXYZ(1)*appArmXYZ(2) + appPlatformMass*appPlatformXYZ(1)*appPlatformXYZ(2);
Ixz = appArmMass*appArmXYZ(1)*appArmXYZ(3) + appPlatformMass*appPlatformXYZ(1)*appPlatformXYZ(3);
Iyz = appArmMass*appArmXYZ(2)*appArmXYZ(3) + appPlatformMass*appPlatformXYZ(2)*appPlatformXYZ(3);

% Complete inertia tensor!
inertiaTensor = [Ix, Ixy, Ixz; Ixy, Iy, Iyz; Ixz, Iyz, Iz];

%% Discretize into surfaces (defined using points CCW)
% Central structure
boxUpper = [boxWidth/2,-boxWidth/2,0; boxWidth/2,boxWidth/2,0; -boxWidth/2,boxWidth/2,0; -boxWidth/2,-boxWidth/2,0];
boxFront = [boxWidth/2, boxWidth/2, 0; boxWidth/2, -boxWidth/2, 0; boxWidth/2, -boxWidth/2,-boxHeight; boxWidth/2,boxWidth/2, -boxHeight];
boxRight = [-boxWidth/2, boxWidth/2, 0; boxWidth/2, boxWidth/2, 0; boxWidth/2, boxWidth/2, -boxHeight; -boxWidth/2, boxWidth/2, -boxHeight];
boxLower = boxUpper; boxLower(:,1) = boxLower(:,1)*-1; boxLower = boxLower - [0,0,boxHeight;0,0,boxHeight;0,0,boxHeight;0,0,boxHeight];
boxRear = [-boxWidth/2, -boxWidth/2, 0; -boxWidth/2, boxWidth/2, 0; -boxWidth/2, boxWidth/2,-boxHeight; -boxWidth/2,-boxWidth/2, -boxHeight];
boxLeft = [boxWidth/2, -boxWidth/2, 0; -boxWidth/2, -boxWidth/2, 0; -boxWidth/2, -boxWidth/2, -boxHeight; boxWidth/2, -boxWidth/2, -boxHeight];

% Right inner solar panel
spriUpper = [solarPanelWidth/2, boxWidth/2,0; solarPanelWidth/2, boxWidth/2 + solarPanelLength,0; -solarPanelWidth/2, boxWidth/2+solarPanelLength, 0; -solarPanelWidth/2, boxWidth/2,0];
spriLower = [-solarPanelWidth/2, boxWidth/2,-solarPanelThickness; -solarPanelWidth/2, boxWidth/2 + solarPanelLength,-solarPanelThickness; solarPanelWidth/2, boxWidth/2+solarPanelLength, -solarPanelThickness; solarPanelWidth/2, boxWidth/2,-solarPanelThickness];
spriFront = [solarPanelWidth/2,boxWidth/2,0; solarPanelWidth/2,boxWidth/2,-solarPanelThickness; solarPanelWidth/2, boxWidth/2+solarPanelLength, -solarPanelThickness; solarPanelWidth/2, boxWidth/2 + solarPanelLength, 0];
spriRear = [-solarPanelWidth/2,boxWidth/2,-solarPanelThickness; -solarPanelWidth/2,boxWidth/2,0; -solarPanelWidth/2, boxWidth/2+solarPanelLength, 0; -solarPanelWidth/2, boxWidth/2 + solarPanelLength, -solarPanelThickness];
spriInner = [solarPanelWidth/2,boxWidth/2,0; -solarPanelWidth/2, boxWidth/2,0; -solarPanelWidth/2, boxWidth/2, -solarPanelThickness; solarPanelWidth/2, boxWidth/2, -solarPanelThickness];
spriOuter = [solarPanelWidth/2,boxWidth/2 + solarPanelLength, 0; -solarPanelWidth/2, boxWidth/2 + solarPanelLength, 0; -solarPanelWidth/2, boxWidth/2 + solarPanelLength, -solarPanelThickness; solarPanelWidth/2, boxWidth/2 + solarPanelLength, -solarPanelThickness];

% Left inner solar panel
spliUpper = [-solarPanelWidth/2, -boxWidth/2,0; -solarPanelWidth/2, -boxWidth/2 - solarPanelLength,0; solarPanelWidth/2, -boxWidth/2-solarPanelLength, 0; solarPanelWidth/2, -boxWidth/2,0];
spliLower = [solarPanelWidth/2, -boxWidth/2,0; solarPanelWidth/2, -boxWidth/2 - solarPanelLength,0; -solarPanelWidth/2, -boxWidth/2-solarPanelLength, 0; -solarPanelWidth/2, -boxWidth/2,0];
spliFront = [solarPanelWidth/2,-boxWidth/2,0; solarPanelWidth/2,-boxWidth/2,-solarPanelThickness; solarPanelWidth/2, -boxWidth/2-solarPanelLength, -solarPanelThickness; solarPanelWidth/2, -boxWidth/2-solarPanelLength, 0];
spliRear = [-solarPanelWidth/2,-boxWidth/2,-solarPanelThickness; -solarPanelWidth/2,-boxWidth/2,0; -solarPanelWidth/2, -boxWidth/2-solarPanelLength, 0; -solarPanelWidth/2, -boxWidth/2 - solarPanelLength, -solarPanelThickness];
spliInner = [solarPanelWidth/2,-boxWidth/2, 0; solarPanelWidth/2, -boxWidth/2, -solarPanelThickness; -solarPanelWidth/2, -boxWidth/2, -solarPanelThickness; -solarPanelWidth/2, -boxWidth/2, 0];
spliOuter = [-solarPanelWidth/2,-boxWidth/2 - solarPanelLength, 0; -solarPanelWidth/2, -boxWidth/2 - solarPanelLength, -solarPanelThickness; solarPanelWidth/2, -boxWidth/2 - solarPanelLength, -solarPanelThickness; solarPanelWidth/2, -boxWidth/2 - solarPanelLength, 0];

% Right outer solar panel
sproUpper = [-solarPanelWidth/2, boxWidth/2+solarPanelLength,0; 
    solarPanelWidth/2, boxWidth/2 + solarPanelLength,0; 
    solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelLength*cosd(20), solarPanelLength*sind(20); 
    -solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelLength*cosd(20), solarPanelLength*sind(20)];
sproLower = [solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20), -solarPanelThickness*cosd(20);
            -solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20), -solarPanelThickness*cosd(20);
            -solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20), solarPanelLength*sind(20) - solarPanelThickness*cosd(20);
            solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20), solarPanelLength*sind(20) - solarPanelThickness*cosd(20)];
sproFront = [solarPanelWidth/2, boxWidth/2 + solarPanelLength, 0;
            solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20), -solarPanelThickness*cosd(20);
            solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20), solarPanelLength*sind(20) - solarPanelThickness*cosd(20);
             solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelLength*cosd(20), solarPanelLength*sind(20)];
sproRear = [-solarPanelWidth/2, boxWidth/2 + solarPanelLength, 0;
            -solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelLength*cosd(20), solarPanelLength*sind(20);
            -solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20), solarPanelLength*sind(20) - solarPanelThickness*cosd(20);
            -solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20), -solarPanelThickness*cosd(20)];
sproInner = [solarPanelWidth/2, boxWidth/2+solarPanelLength, 0;
            -solarPanelWidth/2, boxWidth/2+solarPanelLength, 0;
            -solarPanelWidth/2, boxWidth/2+solarPanelLength+solarPanelThickness*sind(20), -solarPanelThickness*cosd(20);
            solarPanelWidth/2, boxWidth/2+solarPanelLength+solarPanelThickness*sind(20), -solarPanelThickness*cosd(20)];
sproOuter = [-solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelLength*cosd(20), solarPanelLength*sind(20); 
    solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelLength*cosd(20), solarPanelLength*sind(20);
    solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20), solarPanelLength*sind(20) - solarPanelThickness*cosd(20);
    -solarPanelWidth/2, boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20), solarPanelLength*sind(20) - solarPanelThickness*cosd(20)];

% Left outer solar panel
sploUpper = [solarPanelWidth/2, -(boxWidth/2+solarPanelLength),0; 
    -solarPanelWidth/2, -(boxWidth/2 + solarPanelLength),0; 
    -solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelLength*cosd(20)), solarPanelLength*sind(20); 
    solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelLength*cosd(20)), solarPanelLength*sind(20)];
sploLower = [-solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20)), -solarPanelThickness*cosd(20);
            solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20)), -solarPanelThickness*cosd(20);
            solarPanelWidth/2,-(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20)), solarPanelLength*sind(20) - solarPanelThickness*cosd(20);
            -solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20)), solarPanelLength*sind(20) - solarPanelThickness*cosd(20)];
sploFront = [solarPanelWidth/2, -(boxWidth/2 + solarPanelLength), 0;
            solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20)), -solarPanelThickness*cosd(20);
            solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20)), solarPanelLength*sind(20) - solarPanelThickness*cosd(20);
             solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelLength*cosd(20)), solarPanelLength*sind(20)];
sploRear = [-solarPanelWidth/2, -(boxWidth/2 + solarPanelLength), 0;
            -solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelLength*cosd(20)), solarPanelLength*sind(20);
            -solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20)), solarPanelLength*sind(20) - solarPanelThickness*cosd(20);
            -solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20)), -solarPanelThickness*cosd(20)];
sploInner = [solarPanelWidth/2, -(boxWidth/2+solarPanelLength), 0;
            -solarPanelWidth/2, -(boxWidth/2+solarPanelLength), 0;
            -solarPanelWidth/2, -(boxWidth/2+solarPanelLength+solarPanelThickness*sind(20)), -solarPanelThickness*cosd(20);
            solarPanelWidth/2, -(boxWidth/2+solarPanelLength+solarPanelThickness*sind(20)), -solarPanelThickness*cosd(20)];
sploOuter = [solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelLength*cosd(20)), solarPanelLength*sind(20); 
    -solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelLength*cosd(20)), solarPanelLength*sind(20);
    -solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20)), solarPanelLength*sind(20) - solarPanelThickness*cosd(20);
    solarPanelWidth/2, -(boxWidth/2+solarPanelLength + solarPanelThickness*sind(20) + solarPanelLength*cosd(20)), solarPanelLength*sind(20) - solarPanelThickness*cosd(20)];

% High-gain antenna
[hgaBaseX, hgaBaseY] = mapCircleToHexagon(0,0,hgaDiameter);
hga1 = [hgaBaseX(1), hgaBaseY(1), 0; hgaBaseX(2), hgaBaseY(2), 0; 0,0, hgaHeight];
hga2 = [hgaBaseX(2), hgaBaseY(2), 0; hgaBaseX(3), hgaBaseY(3), 0; 0,0, hgaHeight];
hga3 = [hgaBaseX(3), hgaBaseY(3), 0; hgaBaseX(4), hgaBaseY(4), 0; 0,0, hgaHeight];
hga4 = [hgaBaseX(4), hgaBaseY(4), 0; hgaBaseX(5), hgaBaseY(5), 0; 0,0, hgaHeight];
hga5 = [hgaBaseX(5), hgaBaseY(5), 0; hgaBaseX(6), hgaBaseY(6), 0; 0,0, hgaHeight];
hga6 = [hgaBaseX(6), hgaBaseY(6), 0; hgaBaseX(7), hgaBaseY(7), 0; 0,0, hgaHeight];
hga7 = [flip(hgaBaseX), flip(hgaBaseY), zeros(7,1)];

% Payload Platform Arm
[appArmY, appArmZ] = mapCircleToHexagon(appArmXYZ(2), appArmXYZ(3), appArmDiameter);
appArm1 = flip([boxWidth/2+appArmLength, appArmY(2), appArmZ(2); boxWidth/2, appArmY(2), appArmZ(2); 
           boxWidth/2, appArmY(1), appArmZ(1);              boxWidth/2+appArmLength, appArmY(1), appArmZ(1)]);
appArm2 = flip([boxWidth/2+appArmLength, appArmY(3), appArmZ(3); boxWidth/2, appArmY(3), appArmZ(3); 
           boxWidth/2, appArmY(2), appArmZ(2);              boxWidth/2+appArmLength, appArmY(2), appArmZ(2)]);
appArm3 = flip([boxWidth/2+appArmLength, appArmY(4), appArmZ(4); boxWidth/2, appArmY(4), appArmZ(4); 
           boxWidth/2, appArmY(3), appArmZ(3);              boxWidth/2+appArmLength, appArmY(3), appArmZ(3)]);
appArm4 = flip([boxWidth/2+appArmLength, appArmY(5), appArmZ(5); boxWidth/2, appArmY(5), appArmZ(5); 
           boxWidth/2, appArmY(4), appArmZ(4);              boxWidth/2+appArmLength, appArmY(4), appArmZ(4)]);
appArm5 = flip([boxWidth/2+appArmLength, appArmY(6), appArmZ(6); boxWidth/2, appArmY(6), appArmZ(6); 
           boxWidth/2, appArmY(5), appArmZ(5);              boxWidth/2+appArmLength, appArmY(5), appArmZ(5)]);
appArm6 = flip([boxWidth/2+appArmLength, appArmY(7), appArmZ(7); boxWidth/2, appArmY(7), appArmZ(7); 
           boxWidth/2, appArmY(6), appArmZ(6);              boxWidth/2+appArmLength, appArmY(6), appArmZ(6)]);
appArm7 = [(boxWidth/2)*ones(7,1), appArmY, appArmZ];
appArm8 = [(boxWidth/2+appArmLength)*ones(7,1),flip(appArmY),flip(appArmZ)];

% Payload Platform
xOffset = boxWidth/2 + appArmLength; zOffset = appPlatformXYZ(3);
appUpper = [xOffset, -appPlatformSide/2, zOffset+appPlatformSide/2;
    xOffset, appPlatformSide/2, zOffset+appPlatformSide/2;
    xOffset, appPlatformSide/2, zOffset-appPlatformSide/2;
    xOffset, -appPlatformSide/2, zOffset-appPlatformSide/2];
appFront = [xOffset,-appPlatformSide/2,zOffset+appPlatformSide/2;
            xOffset+appPlatformSide, -appPlatformSide/2, zOffset+appPlatformSide/2;
            xOffset+appPlatformSide, appPlatformSide/2, zOffset+appPlatformSide/2;
            xOffset, appPlatformSide/2, zOffset+appPlatformSide/2];
appRight = [xOffset+appPlatformSide, appPlatformSide/2, zOffset+appPlatformSide/2;
            xOffset+appPlatformSide, appPlatformSide/2, zOffset-appPlatformSide/2;
            xOffset, appPlatformSide/2, zOffset-appPlatformSide/2;
            xOffset, appPlatformSide/2, zOffset+appPlatformSide/2];
appLower = [xOffset+appPlatformSide, appPlatformSide/2, zOffset+appPlatformSide/2;
            xOffset+appPlatformSide, -appPlatformSide/2, zOffset+appPlatformSide/2;
            xOffset+appPlatformSide, -appPlatformSide/2, zOffset-appPlatformSide/2;
            xOffset+appPlatformSide, appPlatformSide/2, zOffset-appPlatformSide/2];
appRear = [xOffset,appPlatformSide/2,zOffset-appPlatformSide/2;
            xOffset+appPlatformSide, appPlatformSide/2, zOffset-appPlatformSide/2;
            xOffset+appPlatformSide, -appPlatformSide/2, zOffset-appPlatformSide/2;
            xOffset, -appPlatformSide/2, zOffset-appPlatformSide/2];
appLeft = [xOffset+appPlatformSide, -appPlatformSide/2, zOffset-appPlatformSide/2;
            xOffset+appPlatformSide, -appPlatformSide/2, zOffset+appPlatformSide/2;
            xOffset, -appPlatformSide/2, zOffset+appPlatformSide/2;
            xOffset, -appPlatformSide/2, zOffset-appPlatformSide/2];

% All s/c surfaces
surfaces = {boxUpper,boxFront,boxRight,boxLower,boxRear,boxLeft,...
    spriUpper,spriLower,spriFront,spriRear,spriInner,spriOuter,...
    spliUpper,spliLower,spliFront,spliRear,spliInner,spliOuter...
    sproUpper,sproLower,sproFront,sproRear,sproInner,sproOuter...
    sploUpper,sploLower,sploFront,sploRear,sploInner,sploOuter...
    hga1,hga2,hga3,hga4,hga5,hga6,hga7...
    appArm1,appArm2,appArm3,appArm4,appArm5,appArm6,appArm7,appArm8...
    appUpper,appLower,appFront,appRear,appRight,appLeft};

%% Plot surfaces and calculate geometric properties

% Surface plot
color = [0.25 0.5 0.85];
figure(); hold on; xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); axis equal; title('MAVEN Simplified Geometry');
grid on;

% Results for each individual surface
for i=1:length(surfaces)
    surface = surfaces{i};
    x = surface(:,1); y = surface(:,2); z = surface(:,3);
    h = fill3(x,y,z,color);
    h.FaceAlpha = 0.5;

    [normal, area, centroid] = surfaceProperties(surface);
    disp("Surface " + i)
    fprintf("Unit Normal: X = %.4f, Y = %.4f, Z = %.4f\n", normal(1), normal(2), normal(3));
    fprintf("Area: %.4f\n", area);
    fprintf("Centroid: X = %.4f, Y = %.4f, Z = %.4f\n", centroid(1), centroid(2), centroid(3));
    disp("-------------------")
end

% Plot body axes
quiver3(0,0,0,2,0,0, 'Color', 'r', 'LineWidth', 1); text(2.1,0.8,0, 'X', 'Color', 'r', 'FontSize', 12);
quiver3(0,0,0,0,2,0, 'Color', 'r', 'LineWidth', 1); text(0,2.1,0, 'Y', 'Color', 'r', 'FontSize', 12);
quiver3(0,0,0,0,0,2, 'Color', 'r', 'LineWidth', 1); text(0,0,2.1, 'Z', 'Color', 'r', 'FontSize', 12);
hold off;

% Barycenter calculation
volumes = [appArmVolume, appPlatformVolume, boxVolume, hgaVolume, solarPanelVolume, solarPanelVolume, solarPanelVolume, solarPanelVolume]';
coordinates = [appArmXYZ; appPlatformXYZ; boxXYZ; hgaXYZ; spOuterLeftXYZ; spInnerLeftXYZ; spInnerRightXYZ; spOuterRightXYZ];
barycenterX = sum(coordinates(:,1).*volumes)/sum(volumes);
barycenterY = sum(coordinates(:,2).*volumes)/sum(volumes);
barycenterZ = sum(coordinates(:,3).*volumes)/sum(volumes);
fprintf("Barycenter: X = %.4f , Y = %.4f, Z = %.4f m\n", barycenterX, barycenterY, barycenterZ);

%% Useful functions
% Problem 6 - surface calculations
function [unitVector, area, centroid] = surfaceProperties(surface)
    x = surface(:,1); y = surface(:,2); z = surface(:,3);

    % Compute unit normal
    vec1 = surface(2,:)-surface(1,:); vec2 = surface(3,:) - surface(2,:);
    normal = (cross(vec1,vec2));
    unitVector = normal / norm(normal);

    % Compute area
    area = 0;
    if (max(size(surface)) == 4)  % Rectangle
        base = norm(vec1); height = norm(vec2);
        area = base * height;
    elseif (max(size(surface)) == 3) % Triangle
        base = norm(vec1); diagonal = norm(vec2); height = sqrt(diagonal^2 - (base/2)^2);
        area = base * height;
    elseif (max(size(surface)) == 6) % Hexagon
        sideLength = norm(vec1);
        area = 3 * sqrt(3) * sideLength * sideLength / 2;
    end

    % Compute centroid
    centroid = [mean(x), mean(y), mean(z)];
end

% For problem 6 - discretize round surfaces into hexagons
function [x,y] = mapCircleToHexagon(xCenter, yCenter, diameter)
    increment = 360/6;
    radius = diameter/2;
    x = zeros(7,1); y = zeros(7,1);
    for i = 1:7
        x(i) = xCenter + radius*cosd(increment*(i-1));
        y(i) = yCenter + radius*sind(increment*(i-1));
    end
end