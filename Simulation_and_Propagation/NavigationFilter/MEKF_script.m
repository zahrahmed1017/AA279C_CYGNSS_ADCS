% MEKF Initialization
clear;
load("InertiaData.mat")

q_0 = [0; 0; 0; 1]; % aligned with inertial, scalar last
w_0 = [0; deg2rad(5); 0];

state_0 = [q_0; w_0];
cov_0   = [1e-6, 0,    0,    0,    0,    0;...
           0,    1e-6, 0,    0,    0,    0;...
           0,    0,    1e-6, 0,    0,    0;...
           0,    0,    0,    1e-6^2,0,   0;...
           0,    0,    0,    0,    1e-6^2,0;...
           0,    0,    0,    0,    0,     1e-6^2];
processNoise_0 = cov_0./100;
dt = 5; %seconds
tspan = 0:dt:1000;

% filter = MEKF(state_0, cov_0, processNoise_0, dt, I_p);

%% Run Measurement Update

q_mekf = zeros(length(tspan), 4);

if size(state_0,2) ~= 1
    state_0 = state_0';
end

obj.state         = state_0;
obj.errorState    = [0; 0; 0; state_0(5:7)];
obj.stateCov      = cov_0;
obj.processNoise  = processNoise_0; 
obj.dt            = dt;
obj.I             = I_p;

for i = 1:length(tspan)

   % A. filter.Propagate();
    % A1. QuatKin_Linear_Discrete
    w = obj.errorState(4:6);
    w1 = w(1);
    w2 = w(2);
    w3 = w(3);
    skewSym      = [0,  w3, -w2, w1;...
                   -w3, 0,   w1, w2;...
                    w2, -w1, 0,  w3;...
                   -w1, -w2, -w3, 0];
    obj.quatKin = eye(4,4) + (0.5).* skewSym .* obj.dt;

    % A2. PropagateQuat
    Omega = obj.quatKin;
    qt_1  = obj.state(1:4);
    qt    = Omega * qt_1;
    obj.state(1:4) = qt;  

    obj.errorState(1:3) = [0; 0; 0];

   % B. filter.TimeUpdate();

   % B1. Update STM

   % B11 Update A11
    if size(obj.errorState,2) ~= 1
        obj.errorState = obj.errorState';
    end

    wx = crossMatrix(obj.errorState(4:6));
    obj.A11 = eye(3,3) - (1/2) * obj.dt .* wx ; 

   % B12 Update A12
    obj.A12 = (1/4) * obj.dt .* eye(3,3);

   % B13 Update A22
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

    % Update STM
    obj.STM = [obj.A11, obj.A12; ...
               zeros(3,3), obj.A22];

    % Update Angular Velocity
    obj.errorState(4:6) = obj.A22 * obj.errorState(4:6);
            
   % Update State Covariance
    obj.stateCov = obj.STM * obj.stateCov * obj.STM' + obj.processNoise; 

   % filter.Reset();
    % Push error state to absolute state
    % MRP -> Quat
    % angular velocity in angular state transfers to absolute state
    q_old   = obj.state(1:4); % scalar last here
    dcm_old = quat2dcm(q_old([4 1 2 3])'); % scalar first for matlab functions
    dcm_err = rod2dcm(obj.errorState(1:3)');
    dcm_new = dcm_err * dcm_old;
    q_new   = dcm2quat(dcm_new); % scalar first for matlab functions
    % obj.state(1:4) = q_new([2 3 4 1])'; % Putting scalar last again
    obj.state(1:4) = q_old ./ norm(q_old);
    obj.state(5:7) = obj.errorState(4:6);
            
    % Output

    q_mekf(i, :) = obj.state(1:4)';

end



%% Compare to numerical propagation (torque-free)

qw_0 = [q_0; w_0];
M_vec = [0; 0; 0];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), tspan, qw_0, options);

% % Match quaternion signs:
% for j = 1:length(tspan)
% 
%     if qw_prop(j,2) > 0
%         if q_mekf(j,1) < 0
%             q_mekf(j,:) = q_mekf * -1;
%         end
%     end
% end

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
