% MEKF Initialization
clear;
load("InertiaData.mat")

q_0 = [0; 0; 0; 1]; % aligned with inertial, scalar last
w_0 = [0; deg2rad(5); 0];

state_0 = [q_0; w_0];
cov_0   = [1e-3, 0,    0,    0,    0,    0;...
           0,    1e-3, 0,    0,    0,    0;...
           0,    0,    1e-3, 0,    0,    0;...
           0,    0,    0,    1e-6^2,0,   0;...
           0,    0,    0,    0,    1e-6^2,0;...
           0,    0,    0,    0,    0,     1e-6^2];
processNoise_0 = cov_0./100;
sigma_ss  = deg2rad(0.3);
sigma_mag = 0.1; 
measCov   = [sigma_ss^2, 0,          0,          0,           0,           0;...
             0,          sigma_ss^2, 0,          0,           0,           0;...
             0,          0,          sigma_ss^2, 0,           0,           0;...
             0,          0,          0,          sigma_mag^2, 0,           0;...
             0,          0,          0,          0,           sigma_mag^2, 0;...
             0,          0,          0,          0,           0,           sigma_mag^2];
dt = 5; %seconds
tspan = 0:dt:1000;

% filter = MEKF(state_0, cov_0, processNoise_0, dt, I_p);

%% Run Filter 

q_mekf   = zeros(length(tspan), 4);
cov_mekf = zeros(length(tspan), 3);

if size(state_0,2) ~= 1
    state_0 = state_0';
end

obj.state         = state_0;
obj.errorState    = [0; 0; 0; state_0(5:7)];
obj.stateCov      = cov_0;
obj.processNoise  = processNoise_0; 
obj.measureNoise  = measCov;
obj.dt            = dt;
obj.I             = I_p;

for i = 1:length(tspan)

    %%%%%%%%%%% STEP 1 filter.Propagate() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 1.1 QuatKin_Linear_Discrete
    w = obj.errorState(4:6);
    w1 = w(1);
    w2 = w(2);
    w3 = w(3);
    skewSym      = [0,  w3, -w2, w1;...
                   -w3, 0,   w1, w2;...
                    w2, -w1, 0,  w3;...
                   -w1, -w2, -w3, 0];
    obj.quatKin = eye(4,4) + (0.5).* skewSym .* obj.dt;

    %%%% 1.2 PropagateQuat
    Omega = obj.quatKin;
    qt_1  = obj.state(1:4);
    qt    = Omega * qt_1;
    obj.state(1:4) = qt;  

    obj.errorState(1:3) = [0; 0; 0];

   %%%%%%%%%%%%% STEP 2 filter.TimeUpdate() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%% 2.1. Update STM

   %%%% Update A11
    if size(obj.errorState,2) ~= 1
        obj.errorState = obj.errorState';
    end

    wx = crossMatrix(obj.errorState(4:6));
    obj.A11 = eye(3,3) - (1/2) * obj.dt .* wx ; 

   %%%% Update A12
    obj.A12 = (1/4) * obj.dt .* eye(3,3);

   %%%% Update A22
    Ix = obj.I(1,1);
    Iy = obj.I(2,2);
    Iz = obj.I(3,3);
    wx = obj.errorState(4);
    wy = obj.errorState(5);
    wz = obj.errorState(6);

    eulEqLin = [ 0,                 (Iy - Iz)/Ix * wz, (Iy - Iz)/Ix * wy;...
                (Iz - Ix)/Iy * wz,   0,                (Iz - Ix)/Iy * wx;...
                (Ix - Iy)/Iz * wy,   (Ix - Iy)/Iz * wx, 0];
    obj.A22  = eye(3,3) + obj.dt .* eulEqLin; 

    %%%% Update STM
    obj.STM = [obj.A11, obj.A12; ...
               zeros(3,3), obj.A22];

    %%%% Update Angular Velocity
    obj.errorState(4:6) = obj.A22 * obj.errorState(4:6);
            
    %%%% Update State Covariance
    obj.stateCov = obj.STM * obj.stateCov * obj.STM' + obj.processNoise;

    cov_mekf(i,1) = sqrt(obj.stateCov(1,1));
    cov_mekf(i,2) = sqrt(obj.stateCov(2,2));
    cov_mekf(i,3) = sqrt(obj.stateCov(3,3));

   %%%%%%%%%%%%% STEP 3: filter.MeasurementUpdate() %%%%%%%%%%%%%%%%%%%%%%%

   %%% Get true measurements - 3n x 1 vector where n is the number of
   %%% measurements
   y = getTrueMeasurements();

   %%% Get modeled measurements - 3n x 1 vector where n is the number of
   %%% measurements
   z = getModeledMeasurements();

   %%% Calculate Sensitivity Matrix - 3n x 6 
   H = zeros(length(z), length(obj.errorState));
   for i = 1:length(z)/3

   H = 

   %%% Calculate Kalman Gain 

   %%%%%%%%%%%%% STEP 4: filter.Reset() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Push error state to absolute state
    % MRP -> Quat
    % angular velocity in angular state transfers to absolute state
    q_old          = obj.state(1:4); % scalar last here
    errState       = [obj.errorState(1:3); 1];
    q_new          = quatMul(errState, q_old);
    obj.state(1:4) = q_new / norm(q_new);
    obj.state(5:7) = obj.errorState(4:6);
    q_mekf(i, :)   = obj.state(1:4)';  

end



%% Compare to numerical propagation (torque-free)

qw_0 = [q_0; w_0];
M_vec = [0; 0; 0];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), tspan, qw_0, options);

figure 
hold on
plot(t_q, qw_prop(:,1), 'LineWidth', 2)
plot(t_q, qw_prop(:,2), 'LineWidth', 2)
plot(t_q, qw_prop(:,3), 'LineWidth', 2)
plot(t_q, qw_prop(:,4), 'LineWidth', 2)
grid on;
plot(t_q, sqrt(qw_prop(:,1).^2 + qw_prop(:,2).^2 + qw_prop(:,3).^2 + qw_prop(:,4).^2), 'LineWidth', 2)
legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
title("Numerically Propagated quaternions")
fontsize(15,'points')
xlabel('Time [s]')
% saveas(gcf, "Figures_and_Plots/quat_integration_results_q.png")

figure 
hold on
plot(t_q, q_mekf(:,1),  'LineWidth', 2)
plot(t_q, q_mekf(:,2),  'LineWidth', 2)
plot(t_q, q_mekf(:,3),  'LineWidth', 2)
plot(t_q, q_mekf(:,4),  'LineWidth', 2)
grid on;
plot(t_q, sqrt(q_mekf(:,1).^2 + q_mekf(:,2).^2 + q_mekf(:,3).^2 + q_mekf(:,4).^2), 'LineWidth', 2)
legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
ylim([-1.5 1.5])
fontsize(15,'points')
title("MEKF quaternions")
xlabel('Time [s]')
% saveas(gcf, "Figures_and_Plots/quat_integration_results_q.png")


figure()
title("Quaternion Errors")

subplot(4,1,1)
plot(t_q, q_mekf(:,1) - qw_prop(:,1),  'LineWidth', 2)
ylabel('\Delta q_1')
fontsize(15,'points')
grid on;
subplot(4,1,2)
plot(t_q, q_mekf(:,2) - qw_prop(:,2),  'LineWidth', 2)
ylabel('\Delta q_2')
fontsize(15,'points')
grid on;
subplot(4,1,3)
plot(t_q, q_mekf(:,3) - qw_prop(:,3),  'LineWidth', 2)
ylabel('\Delta q_3')
fontsize(15,'points')
grid on;
subplot(4,1,4)
plot(t_q, q_mekf(:,4) - qw_prop(:,4),  'LineWidth', 2)
ylabel('\Delta q_4')
fontsize(15,'points')
grid on;
xlabel('Time [s]')
% plot(t_q, sqrt(q_mekf(:,1).^2 + q_mekf(:,2).^2 + q_mekf(:,3).^2 + q_mekf(:,4).^2), 'LineWidth', 2)
% legend("q_1", "q2", "q_3", "q_4", "q_{mag}")
% title("MEKF quaternions")

%% Plot covariance errors
n = 50;
figure()
subplot(3,1,1)
plot(tspan(1:n), cov_mekf(1:n,1), 'LineWidth', 2)
hold on;
grid on;
plot(tspan(1:n), -cov_mekf(1:n,1), 'LineWidth', 2)
xlabel('Time [s]')
ylabel('\sigma_{P1}')
fontsize(15,'points')

subplot(3,1,2)
plot(tspan(1:n), cov_mekf(1:n,2), 'LineWidth', 2)
hold on;
grid on;
plot(tspan(1:n), -cov_mekf(1:n,2), 'LineWidth', 2)
xlabel('Time [s]')
ylabel('\sigma_{P2}')
fontsize(15,'points')

subplot(3,1,3)
plot(tspan(1:n), cov_mekf(1:n,3), 'LineWidth', 2)
hold on;
grid on;
plot(tspan(1:n), -cov_mekf(1:n,3), 'LineWidth', 2)
xlabel('Time [s]')
ylabel('\sigma_{P3}')
fontsize(15,'points')

