%% Gravity Gradient Torque

close all; clear;

load("InertiaData.mat")

%% TESTING!!!!!!!!

% Bigger moment of inertia just to sanity check
% I_p = [1000, 0, 0;0, 3000, 0; 0, 0, 4000];

%%% Initial Conditions

%%%% INITIALIZE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
% a   = rE + 3000;     % km
a   = rE + 525;
e   = 0.00162;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

% Convert to ECI position and velocity for orbit propagation 
oe       = [a; e; i; O; w; v];
rv_state = OE2ECI(oe, muE);

%%%% SET SIMULATION TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
numPeriods  = 1;
tspan       = 0 : 60 : T * numPeriods; % simulate once an minute?
% tspan       = 0 : 10 : 100;

%%%% INITIALIZE ATTITUDE AND ANGULAR RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial angular velocity:
w_0         = [0, 0, 0]';

% Initial attitude:
e_vec = [1;1;1]; 
e     = e_vec/ norm(e_vec);
p     = 0; % for PS4-Q1
q_0   = [e(1)*sin(p/2);
         e(2)*sin(p/2);
         e(3)*sin(p/2);
         cos(p/2)];

%%% Run Simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_0   = [q_0; w_0; rv_state];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t_q, state_out] = ode113(@(t,state) PropagateOrbit_and_Attitude_GravityGrad(state, I_p, muE), tspan, state_0, options);

%%% Plot Gravity Gradient Components and Magnitude:
gravityGrad = zeros(length(state_out),4);

for i = 1:length(state_out)

    q = state_out(i,1:4)';
    r = state_out(i,8:10)';

    gg = CalculateGravityGradient(q,r,muE,I_p);
    gg_mag = norm(gg);

    gravityGrad(i,:) = [gg, gg_mag];

end

figure()
subplot(2,2,1)
plot(t_q/T, gravityGrad(:,1),'LineWidth',2)
fontsize(16, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - X component')

subplot(2,2,2)
plot(t_q/T, gravityGrad(:,2),'LineWidth',2)
fontsize(16, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Y component')

subplot(2,2,3)
plot(t_q/T, gravityGrad(:,3),'LineWidth',2)
fontsize(16, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Z component')

subplot(2,2,4)
plot(t_q/T, gravityGrad(:,4),'LineWidth',2)
fontsize(16, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Magnitude')
% saveas(gcf,"Figures_and_Plots/PS4/GravityGradient_partB_largeI_largeR.png")

%% Part D - Initially aligned with RTN frame and initial angular velocity matches mean motion

%%%% INITIALIZE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
a   = rE + 525;
e   = 0;
i   = deg2rad(90);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

% Convert to ECI position and velocity for orbit propagation 
oe       = [a; e; i; O; w; v];
rv_state = OE2ECI(oe, muE);

%%%% SET SIMULATION TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
n           = 2*pi/T;
numPeriods  = 1;
tspan       = 0 : 60 : T * numPeriods; % simulate once an minute?
% tspan       = 0 : 10 : 100;

%%%% INITIALIZE ATTITUDE AND ANGULAR RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial angular velocity:
w_0         = [0, 0, n]';

% Initial attitude:
% Calculate RTN position:
rtn       = rv2rtn(rv_state');
R_eci_rtn = [rtn(1:3); rtn(4:6); rtn(7:9)];
q_0       = dcm2quat(R_eci_rtn);
q_0       = q_0([2 3 4 1])';

%%% Run Simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_0   = [q_0; w_0; rv_state];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t_q, state_out] = ode113(@(t,state) PropagateOrbit_and_Attitude_GravityGrad(state, I_p, muE), tspan, state_0, options);

%%% Plot Gravity Gradient Components and Magnitude:
gravityGrad = zeros(length(state_out),4);

for i = 1:length(state_out)

    q = state_out(i,1:4)';
    r = state_out(i,8:10)';

    gg = CalculateGravityGradient(q,r,muE,I_p);
    gg_mag = norm(gg);

    gravityGrad(i,:) = [gg, gg_mag];

end

figure()
subplot(2,2,1)
plot(t_q/T, gravityGrad(:,1),'LineWidth',2)
fontsize(14, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - X component')

subplot(2,2,2)
plot(t_q/T, gravityGrad(:,2),'LineWidth',2)
fontsize(14, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Y component')

subplot(2,2,3)
plot(t_q/T, gravityGrad(:,3),'LineWidth',2)
fontsize(14, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Z component')

subplot(2,2,4)
plot(t_q/T, gravityGrad(:,4),'LineWidth',2)
fontsize(14, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Magnitude')
% saveas(gcf,"Figures_and_Plots/PS4/GravityGradient_partD_zeroGG.png")

%%% Euler angles:

% Convert to Euler angles: % Going to use 3-2-1
EulerAngs = zeros(length(state_out), 3);

for i = 1:length(state_out)
    quat = state_out(i,1:4);
    [yaw, pitch, roll] = quat2angle(quat([4 1 2 3]), 'ZYX');
    EulerAngs(i,:) = [yaw, pitch, roll];
end

figure()
hold on;
fontsize(14,'points')
plot(t_q/T, rad2deg(EulerAngs(:,1)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs(:,2)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs(:,3)), 'LineWidth', 2)
grid on;
legend('Yaw', 'Pitch', 'Roll')
title("Spacecraft Euler Angles (3-2-1)")
xlabel('Orbits')
ylabel('Angle [deg]')
% saveas(gcf,"Figures_and_Plots/PS4/GravityGradient_partD_zeroGG_eulerAngs.png")


% Sanity check: RTN Euler Angles should match the spacecraft Euler
% angles...

rtn_series = rv2rtn(state_out(:,8:end));

EulerAngs_rtn = zeros(length(state_out), 3);

for i = 1:length(state_out)
    dcm = [rtn_series(i,1:3); rtn_series(i,4:6); rtn_series(i,7:9)];
    [yaw, pitch, roll] = dcm2angle(dcm,'ZYX');
    EulerAngs_rtn(i,:) = [yaw, pitch, roll];
end

figure()
hold on;
fontsize(14,'points')
plot(t_q/T, rad2deg(EulerAngs_rtn(:,1)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs_rtn(:,2)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs_rtn(:,3)), 'LineWidth', 2)
grid on;
legend('Yaw', 'Pitch', 'Roll')
title("RTN Euler Angles (3-2-1)")
xlabel('Orbits')
ylabel('Angle [deg]')

% Another sanity check: RTN to PA Euler angles should be zero....

EulerAngs_rtn_pa = zeros(length(state_out), 3);

for i = 1:length(state_out)
    R_eci_rtn = [rtn_series(i,1:3); rtn_series(i,4:6); rtn_series(i,7:9)];
    quat = state_out(i,1:4);
    R_eci_pa  = quat2dcm(quat([4 1 2 3]));
    R_rtn_pa  = R_eci_rtn' * R_eci_pa;
    [yaw, pitch, roll] = dcm2angle(R_rtn_pa,'ZYX');
    EulerAngs_rtn_pa(i,:) = [yaw, pitch, roll];
end

figure()
hold on;
fontsize(14,'points')
plot(t_q/T, rad2deg(EulerAngs_rtn_pa(:,1)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs_rtn_pa(:,2)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs_rtn_pa(:,3)), 'LineWidth', 2)
grid on;
legend('Yaw', 'Pitch', 'Roll')
title("RTN to PA Euler Angles (3-2-1)")
xlabel('Orbits')
ylabel('Angle [deg]')
ylim([-90 90])
saveas(gcf,"Figures_and_Plots/PS4/GravityGradient_partD_zeroGG_eulerAngsRTN2PA.png")

%% Part E - Initially align BODY AXIS with RTN frame and initial angular velocity matches mean motion

%%%% INITIALIZE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
a   = rE + 525;
e   = 0;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

% Convert to ECI position and velocity for orbit propagation 
oe       = [a; e; i; O; w; v];
rv_state = OE2ECI(oe, muE);

%%%% SET SIMULATION TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
n           = 2*pi/T;
numPeriods  = 3;
tspan       = 0 : 60 : T * numPeriods; % simulate once an minute?
% tspan       = 0 : 10 : 100;

%%%% INITIALIZE ATTITUDE AND ANGULAR RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial angular velocity:
w_0         = [0, 0, n]';

% Initial attitude:
% Calculate RTN position:
rtn       = rv2rtn(rv_state');
R_eci_rtn = [rtn(1:3); rtn(4:6); rtn(7:9)];
R_0       = R_b_p * R_eci_rtn; 
q_0       = dcm2quat(R_0);
q_0       = q_0([2 3 4 1])';

%%% Run Simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_0   = [q_0; w_0; rv_state];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t_q, state_out] = ode113(@(t,state) PropagateOrbit_and_Attitude_GravityGrad(state, I_p, muE), tspan, state_0, options);

%%% Plot Gravity Gradient Components and Magnitude:
gravityGrad = zeros(length(state_out),4);

for i = 1:length(state_out)

    q = state_out(i,1:4)';
    r = state_out(i,8:10)';

    gg = CalculateGravityGradient(q,r,muE,I_p);
    gg_mag = norm(gg);

    gravityGrad(i,:) = [gg, gg_mag];

end

figure()
subplot(2,2,1)
plot(t_q/T, gravityGrad(:,1),'LineWidth',2)
fontsize(14, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - X component')

subplot(2,2,2)
plot(t_q/T, gravityGrad(:,2),'LineWidth',2)
fontsize(14, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Y component')

subplot(2,2,3)
plot(t_q/T, gravityGrad(:,3),'LineWidth',2)
fontsize(14, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Z component')

subplot(2,2,4)
plot(t_q/T, gravityGrad(:,4),'LineWidth',2)
fontsize(14, 'points')
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Magnitude')
% saveas(gcf,"Figures_and_Plots/PS4/GravityGradient_partE_body.png")

%%% Euler angles:

% Convert to Euler angles: % Going to use 3-2-1
EulerAngs = zeros(length(state_out), 3);

for i = 1:length(state_out)
    quat = state_out(i,1:4);
    [yaw, pitch, roll] = quat2angle(quat([4 1 2 3]), 'ZYX');
    EulerAngs(i,:) = [yaw, pitch, roll];
end

figure()
hold on;
fontsize(14,'points')
plot(t_q/T, rad2deg(EulerAngs(:,1)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs(:,2)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs(:,3)), 'LineWidth', 2)
grid on;
legend('Yaw', 'Pitch', 'Roll')
title("Propagated Euler Angles (3-2-1)")
xlabel('Orbits')
ylabel('Angle [deg]')
saveas(gcf,"Figures_and_Plots/PS4/GravityGradient_partE_eulerAngs.png")


rtn_series = rv2rtn(state_out(:,8:end));
EulerAngs_rtn_pa = zeros(length(state_out), 3);

for i = 1:length(state_out)
    R_eci_rtn = [rtn_series(i,1:3); rtn_series(i,4:6); rtn_series(i,7:9)];
    quat = state_out(i,1:4);
    R_eci_pa  = quat2dcm(quat([4 1 2 3]));
    R_rtn_pa  = R_eci_rtn' * R_eci_pa;
    [yaw, pitch, roll] = dcm2angle(R_rtn_pa,'ZYX');
    EulerAngs_rtn_pa(i,:) = [yaw, pitch, roll];
end

figure()
hold on;
fontsize(14,'points')
plot(t_q/T, rad2deg(EulerAngs_rtn_pa(:,1)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs_rtn_pa(:,2)), 'LineWidth', 2)
plot(t_q/T, rad2deg(EulerAngs_rtn_pa(:,3)), 'LineWidth', 2)
grid on;
legend('Yaw', 'Pitch', 'Roll')
title("RTN to PA Euler Angles (3-2-1)")
xlabel('Orbits')
ylabel('Angle [deg]')
ylim([-90 90])
saveas(gcf,"Figures_and_Plots/PS4/GravityGradient_partE_RTN2PA.png")