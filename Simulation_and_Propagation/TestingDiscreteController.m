% Discrete Control torque w/ no actuator model
close all; clear;

load("InertiaData.mat")
load("cygnss.mat")

% Run attitude propagation without any perturbations to show ideal attitude
% propagation and then run with perturbations to show attitude error

%%%% INITIALIZE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
% a   = rE + 3000;     % km
a   = rE + 525;
e   = 0;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

% Calculate ECI and RTN position and velocity for orbit propagation 
oe          = [a; e; i; O; w; v];
rv_state    = OE2ECI(oe, muE);
RTNout      = rv2rtn(rv_state');
A_eci_rtn   = [RTNout(1:3)', RTNout(4:6)', RTNout(7:9)' ]';


%%%% INITIALIZE ATTITUDE AND ANGULAR RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the initial attitude such that principal z is aligned with -R, x is
% aligned with N and y aligned with T
% add some initial error to the rates (TODO add some error to attitude too)
R_eci_nadirRTN = [0 0 1; 0 1 0; -1 0 0] * A_eci_rtn;
q_0 = dcm2quat(R_eci_nadirRTN);
q_0 = q_0([2 3 4 1])';
w_0 = [sqrt(muE/(a^3)) + deg2rad(0.5), deg2rad(0.2), deg2rad(0.5)]';

% Arbitrary initial attitude, for inertial pointing
% q_0 = [0, 0, 0, 1]';
% w_0 = deg2rad([0.5, 0.1, 1]');
% w = w_0; % for the loop

R_i_p = quat2dcm(q_0([4 1 2 3])');

%%%% SET SIMULATION TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_period = 0.1; % sampling frequency is 10 Hz
T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
% numPeriods  = 0.5;
numPeriods = .1;
t_span       = 0 : sample_period : T * numPeriods; % simulate once an minute?
% t_span      = 0:0.5:200;

%%%% INITIALIZE EPOCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];

%%% INITIALIZE REACTION WHEEL QUANTITIES %%%%%%%%
A = [1/sqrt(2), 1/sqrt(3),  1/sqrt(3);...
     0        , 1/sqrt(3), -1/sqrt(3);...
    -1/sqrt(2), 1/sqrt(3),  1/sqrt(3)];
% A     = eye(3,3);
% A     = 1/sqrt(3) * [-1, 1, 1, -1; -1, -1, 1, 1; 1, 1, 1, 1];
% Astar = (A'* A)^(-1) * A';
Astar = pinv(A);
Lw_0  = [0; 0; 0];

% timeStep = .5;

%%%% RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state_0 = [q_0; w_0; rv_state];
M_vec = [0;0;0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
gravityGrad     = 1;
magTorque       = 1;
aeroDrag        = 1;
SRP             = 1;
control         = 1;
discreteControl = 1;
% [t_out, state_out] = ode113(@(t,state) PropagateOrbit_Attitude_wPert_wControl(state, I_p, muE, cygnss, A, Astar, timeStep, initialEpoch, gravityGrad, magTorque, aeroDrag, SRP, control), t_span, state_0, options);

% Get initial control signal
% [az, ay, ax] = quat2angle( state_0([4 1 2 3])', 'ZYX' );
% Mc_rw = controlTorque_inertial_discrete(I_p, [ax; ay; az], state_0(5:7)');
% % TODO add magnetoquers once this is working
% Mact = ComputeActuatorTorque(Lw_0, Mc_rw, w_0, A, Astar);



%% Plot results
Mcontrol  = zeros(length(t_span), 3);
Mactuator = zeros(length(t_span), 3);
Mact_rot  = zeros(length(t_span), 3); 
all_angs  = zeros(length(t_span), 3); 
all_rates = zeros(length(t_span), 3); 
all_angs_rtn = zeros(length(t_span), 3); 

last_time = 0;


% for ii = 1:length(t_out)
for i = 2:length(t_span) % skip time step 0 

    t_span_step = [last_time, t_span(i)];

    % Calculate the actual reaction wheel torque that will be used for control
    reacWheelTorque    = Control2ReacWheelTorque(Lw_0, I_p, R_i_p, w_0, A, Astar, rv_state);
    reacWheelTorque_pa = A * reacWheelTorque;
    
    % Now integrate to find the new angular momentum of the reaction
    % wheels if we were to apply that control torque given our current
    % state and reset Lw_0 so at the next time step we have our
    % updated reaction wheel angular momentum. 
    [t_prop, Lw_out] = ode113(@(t,Lw) Control2ReacWheelTorque(Lw,...
    I_p, R_i_p, w_0, A, Astar, rv_state), t_span_step, Lw_0, options);
    
    Lw_0 = Lw_out(end,:)';

    % TODO add reaction wheel control here
    magnetorquerTorque_pa = [0;0;0];

    appliedTorque_pa = reacWheelTorque_pa + magnetorquerTorque_pa;


    [t_all, state_all] = ode113(@(t,state) IntegratedODE(state_0, I_p, muE, cygnss, t_span(i), ...
        initialEpoch, appliedTorque_pa, gravityGrad, magTorque, aeroDrag, SRP), t_span_step, state_0, options);
    
    ti  = t_all(end);
    qi  = state_all(end, 1:4)';
    wi  = state_all(end, 5:7)';
    rv_state = state_all(end,8:13)';
    w_0 = wi; % for the loop
    state_0 = state_all(end,:)';

    R_i_p = quat2dcm(qi([4 1 2 3])' );

    [az, ay, ax] = quat2angle( qi([4 1 2 3])', 'ZYX' );
    all_angs(i, :) = [ax, ay, az];

    RTNout      = rv2rtn(rv_state');
    R_eci_rtn   = [RTNout(1:3)', RTNout(4:6)', RTNout(7:9)' ]';
    R_rtn_p = R_i_p * R_eci_rtn';
    [az_r, ay_r, ax_r] = dcm2angle( R_rtn_p, 'ZYX' );
    all_angs_rtn(i,:) = [az_r, ay_r, ax_r];

    all_rates(i, :) = wi;

    % M1 = sin(0.1*ti);
    % M2 = sin(0.2*ti);
    % M3 = cos(0.2*ti);
    % Mc = [M1; M2; M3];
    
    % control_angs = [R_i_p(2,3), -R_i_p(1,3), R_i_p(1,2)]'; % assume a 3-2-1 rotation
    % control_angs = [ax; ay; az]; % should get same result as above
    % Mc_rw        = controlTorque_inertial(I_p, control_angs, wi);
    % 
    % Mact = ComputeActuatorTorque(Lw_0, Mc_rw, wi, A, Astar)';
    % Mrot = A * Mact';


    Mactuator(i,:) = reacWheelTorque_pa';
    % Mact_rot(i,:)  = Mrot';
    % Mcontrol(i,:)  = Mc_rw';

    last_time = t_span(i); % so we don't have to run the integrator for long

end

t_hr = t_span/3600;

% figure()
% plot(t_hr, Mcontrol(:,1), 'LineWidth', 2)
% grid on;
% hold on;
% plot(t_hr, Mcontrol(:,2), 'LineWidth', 2);
% plot(t_hr, Mcontrol(:,3), 'LineWidth', 2);
% xlabel('Time [hours]')
% ylabel('Control Torque')
% title('Input Control Torque')
% legend('Mx', 'My', 'Mz')
% fontsize(14,'points')
% 
% figure()
% plot(t_hr, Mact_rot(:,1), 'LineWidth', 2);
% hold on;
% grid on;
% plot(t_hr, Mact_rot(:,2), 'LineWidth', 2);
% plot(t_hr, Mact_rot(:,3), 'LineWidth', 2);
% legend('M_{x,a}', 'M_{y,a}', 'M_{z,a}')
% fontsize(14, 'points')
% xlabel('Time [hours]')
% ylabel('Actuator Torque')
% title('Actuator Torque')

figure()
plot(t_span, Mactuator(:,1), 'LineWidth', 2);
hold on;
grid on;
plot(t_span, Mactuator(:,2), 'LineWidth', 2);
plot(t_span, Mactuator(:,3), 'LineWidth', 2);
legend('M_{1,a}', 'M_{2,a}', 'M_{3,a}')
fontsize(14, 'points')
xlabel('Time, seconds')
ylabel('Actuator Torque')
title('Actuator Torque')
% 
% figure()
% subplot(3,1,1)
% plot(Mcontrol(:,1), Mactuator(:,1), 'LineWidth', 2)
% grid on;
% xlabel('M_{control, x}')
% ylabel('M_{Reaction Wheel 1}')
% fontsize(14,'points')
% 
% subplot(3,1,2)
% plot(Mcontrol(:,2), Mactuator(:,2), 'LineWidth', 2)
% grid on;
% xlabel('M_{control, y}')
% ylabel('M_{Reaction Wheel 2}')
% fontsize(14,'points')
% 
% subplot(3,1,3)
% plot(Mcontrol(:,3), Mactuator(:,3), 'LineWidth', 2)
% grid on;
% xlabel('M_{control, z}')
% ylabel('M_{Reaction Wheel 3}')
% fontsize(14,'points')

% figure()
% title('Rotated')
% subplot(3,1,1)
% plot(Mcontrol(:,1), Mact_rot(:,1), 'LineWidth', 2)
% grid on;
% xlabel('M_{control, x}')
% ylabel('M_{actuator, x}')
% fontsize(14,'points')
% 
% subplot(3,1,2)
% plot(Mcontrol(:,2), Mact_rot(:,2), 'LineWidth', 2)
% grid on;
% xlabel('M_{control, y}')
% ylabel('M_{actuator, y}')
% fontsize(14,'points')
% 
% subplot(3,1,3)
% plot(Mcontrol(:,3), Mact_rot(:,3), 'LineWidth', 2)
% grid on;
% xlabel('M_{control, z}')
% ylabel('M_{actuator, z}')
% fontsize(14,'points')

figure 
hold on
grid on
plot(t_span, rad2deg(all_angs(:,1)), 'LineWidth', 2 )
plot(t_span, rad2deg(all_angs(:,2)), 'LineWidth', 2 )
plot(t_span, rad2deg(all_angs(:,3)), 'LineWidth', 2 )
legend("X rotation", "Y rotation", "Z rotation")
ylabel("Euler angle, degrees")
xlabel("Time, seconds")
title("Attitude, Euler angles")


figure 
hold on
grid on
plot(t_span, rad2deg(all_angs_rtn(:,1)), 'LineWidth', 2 )
plot(t_span, rad2deg(all_angs_rtn(:,2)), 'LineWidth', 2 )
plot(t_span, rad2deg(all_angs_rtn(:,3)), 'LineWidth', 2 )
legend("X rotation", "Y rotation", "Z rotation")
ylabel("Euler angle, degrees")
xlabel("Time, seconds")
title("Attitude in RTN frame, Euler angles")

figure 
hold on
grid on
plot(t_span, rad2deg(all_rates(:,1)), 'LineWidth', 2 )
plot(t_span, rad2deg(all_rates(:,2)), 'LineWidth', 2 )
plot(t_span, rad2deg(all_rates(:,3)), 'LineWidth', 2 )
legend("X rate", "Y rate", "Z rate")
ylabel("Body rates, degrees/s")
xlabel("Time, seconds")
title("Rates")