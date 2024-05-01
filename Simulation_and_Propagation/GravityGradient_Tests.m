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
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - X component')

subplot(2,2,2)
plot(t_q/T, gravityGrad(:,2),'LineWidth',2)
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Y component')

subplot(2,2,3)
plot(t_q/T, gravityGrad(:,3),'LineWidth',2)
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Z component')

subplot(2,2,4)
plot(t_q/T, gravityGrad(:,4),'LineWidth',2)
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Magnitude')
% saveas(gcf,"Figures_and_Plots/PS4/GravityGradient_partB_largeI_largeR.png")

%% Part C - Initially aligned with RTN frame and initial angular velocity matches mean motion

%%%% INITIALIZE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
% a   = rE + 3000;     % km
a   = rE + 525;
% e   = 0.00162; 
e = 0;
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
% Calculate RTN position:
rtn       = rv2rtn(rv_state');
R_eci_rtn = [rtn(1:3); rtn(4:6); rtn(7:9)];




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
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - X component')

subplot(2,2,2)
plot(t_q/T, gravityGrad(:,2),'LineWidth',2)
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Y component')

subplot(2,2,3)
plot(t_q/T, gravityGrad(:,3),'LineWidth',2)
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Z component')

subplot(2,2,4)
plot(t_q/T, gravityGrad(:,4),'LineWidth',2)
grid on;
xlabel('Orbits')
ylabel('Gravity Gradient [N*m]')
title('Gravity Gradient - Magnitude')