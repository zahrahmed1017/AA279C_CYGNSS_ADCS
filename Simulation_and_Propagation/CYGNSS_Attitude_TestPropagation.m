% PS3 - Q6: 
% - Go back to your original satellite inertia tensor. 
% - Numerically integrate Euler AND Kinematic equations from arbitrary initial conditions
% (warning: stay far from singularity of adopted parameterization). 
% Multiple revolutions. The output is the evolution of the attitude 
% parameters over time. These attitude parameters describe orientation of 
% principal axes relative to inertial axes.

close all; clear;

load("InertiaData.mat")

%% Initial conditions
% w_0 = [0.001, 0, deg2rad(5)]'; % rad/s
% w_0 = [ 0, deg2rad(1.5), deg2rad(5)]'; % rad/s (pitch, roll, yaw)
w_0 = [ 0, 0, deg2rad(5)]'; % rad/s (pitch, roll, yaw)
M_vec = [0, 0, 0]';

% initial minute rotation p = 0.001 rad about z axis (to avoid singularity)
e = [1;1;1] / sqrt(3);
% p = 0.5;
p = 0; % for PS4-Q1
q_0 = [e(1)*sin(p/2);
       e(2)*sin(p/2);
       e(3)*sin(p/2);
       cos(p/2)];
dcm_0 = quaternion2dcm(q_0);
% t_span = [0, 10*60];
t_span = 0:3:10*60;
% t_span = 0:10:2500;

q_0_check = dcm2quaternion(dcm_0);


% TODO write a function for this
% ang_0 = quat2eul([q_0(4); q_0(1:3)]', 'ZXZ'); % matlab quaternion convention puts scalar first
ang_0 = dcm2eulerAng(quaternion2dcm(q_0));
% for sanity check: does it convert back the same?
q_test = dcm2quaternion(eulerAng2dcm(ang_0));

%% Propagate the quaternions

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

qw_0 = [q_0; w_0];

[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

figure 
hold on
plot(t_q, qw_prop(:,1), 'LineWidth',2)
plot(t_q, qw_prop(:,2), 'LineWidth', 2)
plot(t_q, qw_prop(:,3), 'LineWidth', 2)
plot(t_q, qw_prop(:,4), 'LineWidth', 2)
grid on;
plot(t_q, sqrt(qw_prop(:,1).^2 + qw_prop(:,2).^2 + qw_prop(:,3).^2 + qw_prop(:,4).^2), 'LineWidth', 2)
legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
title("Propagated quaternions")
saveas(gcf, "Figures_and_Plots/quat_integration_results_q.png")


% plot rates
figure()
hold on;
plot(t_q, rad2deg(qw_prop(:,5)), 'LineWidth',2)
plot(t_q, rad2deg(qw_prop(:,6)), 'LineWidth', 2)
plot(t_q, rad2deg(qw_prop(:,7)), 'LineWidth', 2)
legend('\omega_x', '\omega_y', '\omega_z')
title('Angular rates in PA coordinates')
grid on; 
xlabel("Time, s")
ylabel("Rate, ^\circ /s")
ylim([-0.5, 5.5])
saveas(gcf, "Figures_and_Plots/PS4/Q1a_ang_vel.png")

% check that the two sequences match, convert quat to euler angles
eulerAngs = zeros(length(qw_prop), 3);
for i = 1:length(eulerAngs)
    % dcm = quaternion2dcm(qw_prop(i,1:4));
    % eul = dcm2eulerAng(dcm);
    % [phi, theta, psi] = dcm2angle(dcm,'ZXZ','Robust');
    quat = qw_prop(i,1:4);
    [phi, theta, psi] = quat2angle(quat([4 1 2 3]), 'ZXZ');
    eulerAngs(i,:) = [phi, theta, psi];
end

% plot orientation as Euler angles
figure()
hold on;
plot(t_q, eulerAngs(:,1),'LineWidth', 2);
plot(t_q, eulerAngs(:,2), 'LineWidth', 2);
plot(t_q, eulerAngs(:,3), 'LineWidth', 2);
legend('\phi', '\theta', '\psi')
title('Quaternion Propagation shown as Euler Angles')
grid on; 
% saveas(gcf, "Figures_and_Plots/quat_integration_results_EA.png")
saveas(gcf, "Figures_and_Plots/PS4/Q1a_euler_angs.png"); % for PS4

save("Data/PropAttitude_Quat_Data.mat", 't_q', 'qw_prop', 'eulerAngs')




%% Propagate the DCM:

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

dcm_w_0 = [dcm_0(:,1); dcm_0(:,2); dcm_0(:,3); w_0];

[t_dcm, dcm_w_prop] = ode113(@(t,dcm_w) PropagateAttitude_DCM(dcm_w, M_vec, I_p), t_span, dcm_w_0, options);

eulerAngs_dcm = zeros(length(dcm_w_prop), 3);
for i = 1:length(eulerAngs_dcm)
    dcm_vec = dcm_w_prop(i,1:9)';
    dcm = [dcm_vec(1:3), dcm_vec(4:6), dcm_vec(7:9)]; 
    [phi, theta, psi] = dcm2angle(dcm,'ZXZ','Robust');
    eulerAngs_dcm(i,:) = [phi, theta, psi];
end

figure()
hold on;
plot(t_dcm, eulerAngs_dcm(:,1),'LineWidth', 2);
plot(t_dcm, eulerAngs_dcm(:,2), 'LineWidth', 2);
plot(t_dcm, eulerAngs_dcm(:,3), 'LineWidth', 2);
legend('\phi', '\theta', '\psi')
title('DCM Propagation shown as Euler Angles')
grid on; 
saveas(gcf, "Figures_and_Plots/DCM_integration_results.png")



%% Plot the difference between the quaternion and DCM propagation

phi_diff   = eulerAngs_dcm(:,1) - eulerAngs(:,1);
theta_diff = eulerAngs_dcm(:,2) - eulerAngs(:,2);
psi_diff   = eulerAngs_dcm(:,3) - eulerAngs(:,3);

figure()
subplot(3,1,1)
plot(t_q, phi_diff, 'LineWidth', 2)
title('Phi Error')
grid on;
subplot(3,1,2)
plot(t_q, theta_diff, 'LineWidth', 2)
title('Theta Error')
grid on;
subplot(3,1,3)
plot(t_q, psi_diff, 'LineWidth', 2)
title('Psi Error')
grid on;
saveas(gcf, "Figures_and_Plots/compare_quat_DCM.png")

%% Propagate the Euler angles
% 
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% 
% ang_w_0 = [ang_0; w_0];
% 
% [t_ang, ang_w_prop] = ode113(@(t,ang_w) PropagateAttitude_EulerAng(ang_w, M_vec, I_p), t_span, ang_w_0, options);
% 
% figure()
% hold on;
% plot(t_ang, ang_w_prop(:,1),'LineWidth', 2);
% plot(t_ang, ang_w_prop(:,2), 'LineWidth', 2);
% plot(t_ang, ang_w_prop(:,3), 'LineWidth', 2);
% grid on; 


%% Checking Utility functions for converting attitude representations:

% % First checking the quaternion2dcm function:
% dcm_1 = quaternion2dcm(q_0);
% q1 = q_0(1);
% q2 = q_0(2);
% q3 = q_0(3);
% q4 = q_0(4);
% dcm_2 = [q4^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4);...
%          2*(q1*q2 - q3*q4),          q4^2 - q1^2 + q2^2 - q3^2,  2*(q2*q3 + q1*q4);...
%          2*(q1*q3 + q2*q4),          2*(q2*q3 - q1*q4),          q4^2 + q1^2 - q2^2 - q3^2];
% 
% 
% 
% % Second check the dcm2eulerAng function:
% angs1 = dcm2eulerAng(dcm_1); % Our function
% [phi, theta, psi] = dcm2angle(dcm_1, "ZXZ"); % MATLAB function
% 
% % Third check the eulerAng2dcm:
% dcm_check  = eulerAng2dcm(angs1); % Our function
% dcm_check2 = angle2dcm(angs1(1),angs1(2),angs1(3),'ZXZ'); % MATLAB function
% 
% % Fourth check dcm2quaternion:
% q_check1 = dcm2quaternion(dcm_check);
% q_check2 = dcm2quat(dcm_check);