%% PS4 - Q2 - Stability tests

close all; clear;

load("InertiaData.mat")
load("OrbitPropData_equil.mat")

%% Universal stuff

M_vec = [0, 0, 0]';

e = [1;1;1] / sqrt(3);
p = 0; % for PS4-Q1
q_0 = [e(1)*sin(p/2);
       e(2)*sin(p/2);
       e(3)*sin(p/2);
       cos(p/2)];

w_bar = deg2rad(5);
d_w = 0.001;

t_span = 0:1:15*60;

options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

%% Rotation about x axis

w_0 = [ w_bar, 0, 0]'; % rad/s (pitch, roll, yaw)

w_0_p = w_0 + d_w;

qw_0 = [q_0; w_0_p];
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

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

saveas(gcf, "Figures_and_Plots/PS4/Q2a_ang_vel_xpert.png")

% get euler angles
eulerAngs = zeros(length(qw_prop), 3);
for i = 1:length(eulerAngs)
    quat = qw_prop(i,1:4);
    [phi, theta, psi] = quat2angle(quat([4 1 2 3]), 'ZYX');
    eulerAngs(i,:) = [phi, theta, psi];
end

% plot orientation as Euler angles
figure()
hold on;
plot(t_q, rad2deg(eulerAngs(:,1)),'LineWidth', 2);
plot(t_q, rad2deg(eulerAngs(:,2)), 'LineWidth', 2);
plot(t_q, rad2deg(eulerAngs(:,3)), 'LineWidth', 2);
legend('\phi', '\theta', '\psi')
title('Quaternion Propagation shown as ZYX Euler Angles')
grid on; 
xlabel("Time, s")
ylabel("Angle, ^\circ")
saveas(gcf, "Figures_and_Plots/PS4/Q2a_euler_angs_xpert.png"); % for PS4

%% Rotation about y axis

w_0 = [  0, w_bar, 0]'; % rad/s (pitch, roll, yaw)

w_0_p = w_0 + d_w;

qw_0 = [q_0; w_0_p];
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

% rotate rates to inertial frame
w_i =zeros(length(qw_prop), 3);


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
saveas(gcf, "Figures_and_Plots/PS4/Q2a_ang_vel_ypert.png")

% get euler angles
eulerAngs = zeros(length(qw_prop), 3);
for i = 1:length(eulerAngs)
    quat = qw_prop(i,1:4);
    [phi, theta, psi] = quat2angle(quat([4 1 2 3]), 'ZYX');
    eulerAngs(i,:) = [phi, theta, psi];
end

% plot orientation as Euler angles
figure()
hold on;
plot(t_q, rad2deg(eulerAngs(:,1)),'LineWidth', 2);
plot(t_q, rad2deg(eulerAngs(:,2)), 'LineWidth', 2);
plot(t_q, rad2deg(eulerAngs(:,3)), 'LineWidth', 2);
legend('\phi', '\theta', '\psi')
title('Quaternion Propagation shown as ZYX Euler Angles')
grid on; 
xlabel("Time, s")
ylabel("Angle, ^\circ")
saveas(gcf, "Figures_and_Plots/PS4/Q2a_euler_angs_ypert.png"); % for PS4

%% Rotation about z axis

w_0 = [  0, 0, w_bar]'; % rad/s (pitch, roll, yaw)

w_0_p = w_0 + d_w;

qw_0 = [q_0; w_0_p];
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

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
saveas(gcf, "Figures_and_Plots/PS4/Q2a_ang_vel_zpert.png")

% get euler angles
eulerAngs = zeros(length(qw_prop), 3);
for i = 1:length(eulerAngs)
    quat = qw_prop(i,1:4);
    [phi, theta, psi] = quat2angle(quat([4 1 2 3]), 'ZYX');
    eulerAngs(i,:) = [phi, theta, psi];
end

% plot orientation as Euler angles
figure()
hold on;
plot(t_q, rad2deg(eulerAngs(:,1)),'LineWidth', 2);
plot(t_q, rad2deg(eulerAngs(:,2)), 'LineWidth', 2);
plot(t_q, rad2deg(eulerAngs(:,3)), 'LineWidth', 2);
legend('\phi', '\theta', '\psi')
title('Quaternion Propagation shown as ZYX Euler Angles')
grid on; 
xlabel("Time, s")
ylabel("Angle, ^\circ")
saveas(gcf, "Figures_and_Plots/PS4/Q2a_euler_angs_zpert.png"); % for PS4

