% Attitude Control Error with No Controller
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
state       = OE2ECI(oe, muE);
RTNout      = rv2rtn(state');
A_eci_rtn   = [RTNout(1:3)', RTNout(4:6)', RTNout(7:9)' ]';


%%%% INITIALIZE ATTITUDE AND ANGULAR RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the initial attitude such that principal z is aligned with -R, x is
% aligned with N and y aligned with T
R_eci_nadirRTN = [0 0 1; 0 1 0; -1 0 0] * A_eci_rtn;
q_0 = dcm2quat(R_eci_nadirRTN);
q_0 = q_0([2 3 4 1])';
w_0 = [ sqrt(muE/(a^3)), 0, 0]';

%%%% SET SIMULATION TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
numPeriods  = 2;
t_span       = 0 : 60 : T * numPeriods; % simulate once an minute?

%%%% INITIALIZE EPOCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];

%%%% RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qw_0 = [q_0; w_0];
M_vec = [0;0;0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_nopert, state_nopert] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

%%%% Calculate Euler Angles with zero perturbations %%%%%%%%%%%%%%%%%%%%%%
eulerAngs_noPert = zeros(length(t_nopert), 3);
for i = 1:length(t_nopert)

    q_i_pert   = state_nopert(i,1:4);
    dcm_i = quat2dcm(q_i_pert([4 1 2 3]));
    [yaw, pitch, roll] = quat2angle(q_i_pert([4 1 2 3]),'ZYX');

    eulerAngs_noPert(i,:) = [yaw, pitch, roll];


end

t_hr = t_nopert / 3600;

figure()
subplot(3,1,1)
plot(t_hr, rad2deg(eulerAngs_noPert(:,1)), 'LineWidth', 2)
grid on
xlabel('Time [hr]')
ylabel('Yaw [degrees]')
fontsize(14,'points')
ylim([0 50])
title('Yaw')

subplot(3,1,2)
plot(t_hr, rad2deg(eulerAngs_noPert(:,2)), 'LineWidth', 2)
grid on;
xlabel('Time [hr]')
ylabel('Pitch [degrees]')
fontsize(14,'points')
ylim([-70 -40])
title('Pitch')

subplot(3,1,3)
plot(t_hr, rad2deg(eulerAngs_noPert(:,3)), 'LineWidth', 2)
grid on;
xlabel('Time [hr]')
ylabel('Roll [degrees]')
title('Roll')
fontsize(14,'points')

%% Now include perturbations:

state_0     = [q_0; w_0; state];
gravityGrad = 1;
magTorque   = 1;
aeroDrag    = 1;
SRP         = 1;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t_pert, state_pert] = ...
    ode113(@(t,state) PropagateOrbit_and_Attitude_wPerturbations(state, I_p, muE, cygnss, t, initialEpoch, gravityGrad, magTorque, aeroDrag, SRP),...
    t_span, state_0, options);

eulerAngs_pert  = zeros(length(t_nopert), 3);
eulerAngs_delta = zeros(length(t_nopert), 3);

for j = 1 : length(t_pert)

    q_i_pert   = state_pert(j,1:4);
    dcm_i      = quat2dcm(q_i_pert([4 1 2 3]));
    [yaw, pitch, roll] = quat2angle(q_i_pert([4 1 2 3]),'ZYX');

    eulerAngs_pert(j,:) = [yaw, pitch, roll];

    
    q_i_nopert = state_nopert(j,1:4);
    q_delta    = quatmultiply(quatinv(q_i_nopert([4 1 2 3])), q_i_pert([4 1 2 3]));
    [yaw_d, pitch_d, roll_d] = quat2angle(q_delta, 'ZYX');
    eulerAngs_delta(j,:)     = [yaw_d, pitch_d, roll_d];

end

figure()
subplot(3,1,1)
plot(t_hr, rad2deg(eulerAngs_pert(:,1)), 'LineWidth', 2)
grid on
xlabel('Time [hr]')
ylabel('Yaw [degrees]')
fontsize(14,'points')
% ylim([0 50])
title('Yaw')

subplot(3,1,2)
plot(t_hr, rad2deg(eulerAngs_pert(:,2)), 'LineWidth', 2)
grid on;
xlabel('Time [hr]')
ylabel('Pitch [degrees]')
fontsize(14,'points')
% ylim([-70 -40])
title('Pitch')

subplot(3,1,3)
plot(t_hr, rad2deg(eulerAngs_pert(:,3)), 'LineWidth', 2)
grid on;
xlabel('Time [hr]')
ylabel('Roll [degrees]')
title('Roll')
fontsize(14,'points')

figure()
subplot(3,1,1)
plot(t_hr, rad2deg(eulerAngs_delta(:,1)), 'LineWidth', 2)
grid on
xlabel('Time [hr]')
ylabel('\Delta [degrees]')
fontsize(14,'points')
% ylim([0 50])
title('Delta Yaw')

subplot(3,1,2)
plot(t_hr, rad2deg(eulerAngs_delta(:,2)), 'LineWidth', 2)
grid on;
xlabel('Time [hr]')
ylabel('\Delta [degrees]')
fontsize(14,'points')
% ylim([-70 -40])
title('Delta Pitch')

subplot(3,1,3)
plot(t_hr, rad2deg(eulerAngs_delta(:,3)), 'LineWidth', 2)
grid on;
xlabel('Time [hr]')
ylabel('\Delta [degrees]')
title('Delta Roll')
fontsize(14,'points')

