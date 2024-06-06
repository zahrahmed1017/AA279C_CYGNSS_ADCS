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
w_0 = [sqrt(muE/(a^3)), 0, 0]';

%%%% SET SIMULATION TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
numPeriods  = 0.5;
% t_span       = 0 : 10 : T * numPeriods; % simulate once an minute?
t_span      = 0:0.5:200;

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


%%%% RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state_0 = [q_0; w_0; Lw_0];
M_vec = [0;0;0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, state_out] = ode113(@(t,state) PropagateAttitude_Quat_wControl(state, t, I_p, A, Astar), t_span, state_0, options);

%% Plot results
Mcontrol  = zeros(length(t_out), 3);
Mactuator = zeros(length(t_out), 3);
Mact_rot  = zeros(length(t_out), 3); 

for ii = 1:length(t_out)
    
    ti  = t_out(ii);
    wi  = state_out(ii, 5:7)';
    Lwi = state_out(ii, 8:10)';
    
    M1 = sin(0.1*ti);
    M2 = sin(0.2*ti);
    M3 = cos(0.2*ti);
    Mc = [M1; M2; M3];

    Mact = ComputeActuatorTorque(Lwi, Mc, wi, A, Astar)';
    Mrot = A*Mact';

    Mactuator(ii,:) = Mact';
    Mact_rot(ii,:)  = Mrot';
    Mcontrol(ii,:)  = Mc';

end

t_min = t_out/60;

figure()
plot(t_min, Mcontrol(:,1), 'LineWidth', 2)
grid on;
hold on;
plot(t_min, Mcontrol(:,2), 'LineWidth', 2);
plot(t_min, Mcontrol(:,3), 'LineWidth', 2);
xlabel('Time [min]')
ylabel('Control Torque')
title('Input Control Torque')
legend('Mx', 'My', 'Mz')
fontsize(14,'points')

figure()
plot(t_min, Mact_rot(:,1), 'LineWidth', 2);
hold on;
grid on;
plot(t_min, Mact_rot(:,2), 'LineWidth', 2);
plot(t_min, Mact_rot(:,3), 'LineWidth', 2);
legend('M_{x,a}', 'M_{y,a}', 'M_{z,a}')
fontsize(14, 'points')
xlabel('Time [min]')
ylabel('Actuator Torque')
title('Actuator Torque')

figure()
plot(t_min, Mactuator(:,1), 'LineWidth', 2);
hold on;
grid on;
plot(t_min, Mactuator(:,2), 'LineWidth', 2);
plot(t_min, Mactuator(:,3), 'LineWidth', 2);
legend('M_{1,a}', 'M_{2,a}', 'M_{3,a}')
fontsize(14, 'points')
xlabel('Time [min]')
ylabel('Actuator Torque')
title('Actuator Torque')

figure()
subplot(3,1,1)
plot(Mcontrol(:,1), Mactuator(:,1), 'LineWidth', 2)
grid on;
xlabel('M_{control, x}')
ylabel('M_{Reaction Wheel 1}')
fontsize(14,'points')

subplot(3,1,2)
plot(Mcontrol(:,2), Mactuator(:,2), 'LineWidth', 2)
grid on;
xlabel('M_{control, y}')
ylabel('M_{Reaction Wheel 2}')
fontsize(14,'points')

subplot(3,1,3)
plot(Mcontrol(:,3), Mactuator(:,3), 'LineWidth', 2)
grid on;
xlabel('M_{control, z}')
ylabel('M_{Reaction Wheel 3}')
fontsize(14,'points')

figure()
title('Rotated')
subplot(3,1,1)
plot(Mcontrol(:,1), Mact_rot(:,1), 'LineWidth', 2)
grid on;
xlabel('M_{control, x}')
ylabel('M_{actuator, x}')
fontsize(14,'points')

subplot(3,1,2)
plot(Mcontrol(:,2), Mact_rot(:,2), 'LineWidth', 2)
grid on;
xlabel('M_{control, y}')
ylabel('M_{actuator, y}')
fontsize(14,'points')

subplot(3,1,3)
plot(Mcontrol(:,3), Mact_rot(:,3), 'LineWidth', 2)
grid on;
xlabel('M_{control, z}')
ylabel('M_{actuator, z}')
fontsize(14,'points')
