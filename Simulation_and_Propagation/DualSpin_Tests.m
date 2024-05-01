%% Dual Spin Satellite:

% PS4 - Q3:

%{
Adding a momentum wheel or rotor (dual-spin satellite)
 a. Re-program Euler equations to include a generic momentum wheel or rotor
    with rotation axis aligned with one of the principal axes of inertia. 
 b. Numerically integrate Euler AND Kinematic equations from equilibrium 
    initial condition. Verify that integration is correct as from previous 
    tests (conservation laws, rotations, etc.).
 c. Verify equilibrium and its stability similar to previous pset.
 d. Use the stability condition to make attitude motion stable for rotation
    about intermediate moment of inertia by changing moment of inertia and/or 
    angular velocity of the momentum wheel or rotor.
 e. Try to make rotation about another arbitrary axis (potentially relevant
    to your project) stable through a generic momentum wheel or rotor
%}


close all; clear; 

load("InertiaData.mat")

%% Define Initial Conditions

% Reaction Wheel
m_rotor = 3; %kg
R_rotor = 0.2; % radius [m]
I_rotor = 0.5 * m_rotor * (R_rotor)^2;
r_rotor = [0,0,1]; % rotor axis of rotation

% Angular velocity
w_0     = [0; deg2rad(2); deg2rad(5); deg2rad(50)]; % [wx, wy, wz, wr]

% Initial Attitiude
e       = [1;1;1] / sqrt(3);
p       = 0; % for PS4-Q1, spacecraft initially aligned with inertial frame;
q_0     = [e(1)*sin(p/2);
           e(2)*sin(p/2);
           e(3)*sin(p/2);
           cos(p/2)];
dcm_0   = quaternion2dcm(q_0);
t_span  = 0:0.5:10*60;

% External Torques
M_vec = [0, 0, 0];
Mr    = 0;

%% Propagate the quaternions

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
qw_0    = [q_0; w_0];

[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_DualSpin(qw, M_vec, I_p, Mr, I_rotor, r_rotor), t_span, qw_0, options);

figure 
hold on
fontsize(14,'points')
plot(t_q, qw_prop(:,1), 'LineWidth',2)
plot(t_q, qw_prop(:,2), 'LineWidth', 2)
plot(t_q, qw_prop(:,3), 'LineWidth', 2)
plot(t_q, qw_prop(:,4), 'LineWidth', 2)
grid on;
plot(t_q, sqrt(qw_prop(:,1).^2 + qw_prop(:,2).^2 + qw_prop(:,3).^2 + qw_prop(:,4).^2), 'LineWidth', 2)
legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
title("Propagated quaternions")
xlabel('Time [s]')
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partB_quaternion.png")

% Convert to Euler angles: % Going to use 3-2-1
EulerAngs = zeros(length(qw_prop), 3);

for i = 1:length(qw_prop)
    quat = qw_prop(i,1:4);
    [yaw, pitch, roll] = quat2angle(quat([4 1 2 3]), 'ZYX');
    EulerAngs(i,:) = [yaw, pitch, roll];
end

figure()
hold on;
fontsize(14,'points')
plot(t_q, EulerAngs(:,1), 'LineWidth', 2)
plot(t_q, EulerAngs(:,2), 'LineWidth', 2)
plot(t_q, EulerAngs(:,3), 'LineWidth', 2)
grid on;
legend('Yaw', 'Pitch', 'Roll')
title("Propagated Euler Angles (3-2-1)")
xlabel('Time [s]')
ylabel('Angle [rad]')
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partB_eulerAngs.png")

%% Validate correctness of propagation by checking if angular momentum is conserved in inertial frame

w_sat    = qw_prop(:,5:7);
w_rot    = qw_prop(:,8);
L_sat_pa = (I_p * w_sat')';
L_rot_pa = (I_rotor * w_rot * r_rotor);

L_tot_pa = L_sat_pa + L_rot_pa;

[n, ~] = size(qw_prop);
L_tot_i = zeros(n, 3);
for i = 1:n

    q_i_pa = qw_prop(i,1:4)'; % attitude quaternion for this time step
    A_i_pa = quaternion2dcm(q_i_pa); % rotation from inertial frame to PA

    L_tot_i(i, :) = (A_i_pa' * L_tot_pa(i, :)')'; % transpose A because we want to go from PA to inertial
 
end

% Angular momentum over time
figure 
hold on
grid on
fontsize(14,'points')
plot(t_q, rad2deg(L_tot_i(:,1)), 'LineWidth',2)
plot(t_q, rad2deg(L_tot_i(:,2)), 'LineWidth',2)
plot(t_q, rad2deg(L_tot_i(:,3)), 'LineWidth',2)
title("Components of angular momentum in inertial coordinates")
xlabel("Time, s")
ylabel("Magnitude, kg m^2/s")
legend("L_1", "L_2", "L_3")
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partB_constantL.png")

%% Verify Equilibrium

% For equilibrium, w of the satellite should be aligned with w of rotor:
r_rotor = [0,0,1]; % rotor axis of rotation
w_sat   = [0; 0; deg2rad(5)];
w_r     = deg2rad(50); 

% Angular velocity
w_0     = [w_sat; w_r]; % [wx, wy, wz, wr]

% Using the same other initial conditions as first part of this script
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
qw_0    = [q_0; w_0];

[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_DualSpin(qw, M_vec, I_p, Mr, I_rotor, r_rotor), t_span, qw_0, options);

% Confirm the satellite is in equilibrium by looking at the angular
% velocity in inertial frame (components should be constant)
[n, ~] = size(qw_prop);
w_pa  = qw_prop(:,5:7);
w_i_i = zeros(n, 3);
for i = 1:n
    q_i_pa = qw_prop(i,1:4)'; % attitude quaternion for this time step
    A_i_pa = quaternion2dcm(q_i_pa); % rotation from inertial frame to PA
    w_i_i(i, :) = (A_i_pa' * rad2deg(w_pa(i, :))')'; 
end

% Angular velocity over time
figure 
subplot(3,1,1)
fontsize(14, 'points')
grid on;
plot(t_q, w_i_i(:,1), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
title('w_x')

subplot(3,1,2)
fontsize(14, 'points')
grid on;
plot(t_q, w_i_i(:,2), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
title('w_y')

subplot(3,1,3)
grid on;
plot(t_q, w_i_i(:,3), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
fontsize(14, 'points')
ylim([0 10])
title('w_z')

saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partC_equilibrium.png")

% Convert to Euler angles: % Going to use 3-2-1
EulerAngs = zeros(length(qw_prop), 3);

for i = 1:length(qw_prop)
    quat = qw_prop(i,1:4);
    [yaw, pitch, roll] = quat2angle(quat([4 1 2 3]), 'ZYX');
    EulerAngs(i,:) = [yaw, pitch, roll];
end

figure()
hold on;
fontsize(14,'points')
plot(t_q, EulerAngs(:,1), 'LineWidth', 2)
plot(t_q, EulerAngs(:,2), 'LineWidth', 2)
plot(t_q, EulerAngs(:,3), 'LineWidth', 2)
grid on;
legend('Yaw', 'Pitch', 'Roll')
title("Propagated Euler Angles (3-2-1)")
xlabel('Time [s]')
ylabel('Angle [rad]')
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partC_equilibrium_eulerAngs.png")

% figure 
% hold on
% plot(t_q, qw_prop(:,1), 'LineWidth',2)
% plot(t_q, qw_prop(:,2), 'LineWidth', 2)
% plot(t_q, qw_prop(:,3), 'LineWidth', 2)
% plot(t_q, qw_prop(:,4), 'LineWidth', 2)
% grid on;
% plot(t_q, sqrt(qw_prop(:,1).^2 + qw_prop(:,2).^2 + qw_prop(:,3).^2 + qw_prop(:,4).^2), 'LineWidth', 2)
% legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
% title("Propagated quaternions")

%% Stability

% For equilibrium, w of the satellite should be aligned with w of rotor:
r_rotor = [0,0,1]; % rotor axis of rotation
delta_w = deg2rad(0.1);
w_sat   = [delta_w; delta_w; deg2rad(5) + delta_w];
w_r     = deg2rad(50);

% Angular velocity
w_0     = [w_sat; w_r]; % [wx, wy, wz, wr]

% Using the same other initial conditions as first part of this script
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
qw_0    = [q_0; w_0];

[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_DualSpin(qw, M_vec, I_p, Mr, I_rotor, r_rotor), t_span, qw_0, options);

% Confirm the satellite is in equilibrium by looking at the angular
% velocity in inertial frame (components should be constant)

[n, ~] = size(qw_prop);
w_pa  = qw_prop(:,5:7);
w_i_i = zeros(n, 3);
for i = 1:n
    q_i_pa = qw_prop(i,1:4)'; % attitude quaternion for this time step
    A_i_pa = quaternion2dcm(q_i_pa); % rotation from inertial frame to PA
    w_i_i(i, :) = (A_i_pa' * rad2deg(w_pa(i, :))')'; 
end

% Angular velocity over time
figure 
subplot(3,1,1)
fontsize(14, 'points')
grid on;
plot(t_q, w_i_i(:,1), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
title('w_x')

subplot(3,1,2)
fontsize(14, 'points')
grid on;
plot(t_q, w_i_i(:,2), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
title('w_y')

subplot(3,1,3)
grid on;
plot(t_q, w_i_i(:,3), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
fontsize(14, 'points')
ylim([0 10])
title('w_z')
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partC_stability.png")

% figure 
% hold on
% plot(t_q, qw_prop(:,1), 'LineWidth',2)
% plot(t_q, qw_prop(:,2), 'LineWidth', 2)
% plot(t_q, qw_prop(:,3), 'LineWidth', 2)
% plot(t_q, qw_prop(:,4), 'LineWidth', 2)
% grid on;
% plot(t_q, sqrt(qw_prop(:,1).^2 + qw_prop(:,2).^2 + qw_prop(:,3).^2 + qw_prop(:,4).^2), 'LineWidth', 2)
% legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
% title("Propagated quaternions")

% Convert to Euler angles: % Going to use 3-2-1
EulerAngs = zeros(length(qw_prop), 3);

for i = 1:length(qw_prop)
    quat = qw_prop(i,1:4);
    [yaw, pitch, roll] = quat2angle(quat([4 1 2 3]), 'ZYX');
    EulerAngs(i,:) = [yaw, pitch, roll];
end

figure()
hold on;
fontsize(14,'points')
plot(t_q, EulerAngs(:,1), 'LineWidth', 2)
plot(t_q, EulerAngs(:,2), 'LineWidth', 2)
plot(t_q, EulerAngs(:,3), 'LineWidth', 2)
grid on;
legend('Yaw', 'Pitch', 'Roll')
title("Propagated Euler Angles (3-2-1)")
xlabel('Time [s]')
ylabel('Angle [rad]')
saveas(gcf,"Figures_and_Plots/PS4/DualSpin_partC_stability_eulerAngs.png")

%% Stability condition to make it rotational motion stable about intermediate moment of inertia

% For equilibrium, w of the satellite should be aligned with w of rotor:
r_rotor = [0,1,0]; % rotor axis of rotation
delta_w = deg2rad(15);
w_sat   = [0 + delta_w; deg2rad(5) + delta_w; 0 + delta_w];
w_r     = -2 ; % -0.2 rad/s is unstable, deg2rad(50) is stable

% Angular velocity
w_0     = [w_sat; w_r]; % [wx, wy, wz, wr]

% Check stability:
K1 = (I_p(2,2) - I_p(3,3)) * w_sat(3) / I_rotor;
K2 = (I_p(1,1) - I_p(3,3)) * w_sat(3) / I_rotor;

fprintf('Stability criteria 1: %.4g \n',K1)
fprintf('Stability criteria 1: %.4g \n',K2)
fprintf('Rotor angular velocity: %.4g \n', w_r)

% Using the same other initial conditions as first part of this script
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
qw_0    = [q_0; w_0];

[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_DualSpin(qw, M_vec, I_p, Mr, I_rotor, r_rotor), t_span, qw_0, options);

% Confirm the satellite is in equilibrium by looking at the angular
% velocity in inertial frame (components should be constant)

[n, ~] = size(qw_prop);
w_pa  = qw_prop(:,5:7);
w_i_i = zeros(n, 3);
for i = 1:n

    q_i_pa = qw_prop(i,1:4)'; % attitude quaternion for this time step
    A_i_pa = quaternion2dcm(q_i_pa); % rotation from inertial frame to PA
    w_i_i(i, :) = (A_i_pa' * rad2deg(w_pa(i, :))')'; 
end

% Angular velocity over time
figure 
subplot(3,1,1)
fontsize(14, 'points')
grid on;
plot(t_q, w_i_i(:,1), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
title('w_x')

subplot(3,1,2)
fontsize(14, 'points')
grid on;
plot(t_q, w_i_i(:,2), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
title('w_y')

subplot(3,1,3)
grid on;
plot(t_q, w_i_i(:,3), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
fontsize(14, 'points')
title('w_z')
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partD_intermediate_unstable.png")

% Convert to Euler angles: % Going to use 3-2-1
EulerAngs = zeros(length(qw_prop), 3);

for i = 1:length(qw_prop)
    quat = qw_prop(i,1:4);
    [pitch, yaw, roll] = quat2angle(quat([4 1 2 3]), 'YZX');
    EulerAngs(i,:) = [yaw, pitch, roll];
end

figure 
subplot(3,1,1)
fontsize(14, 'points')
grid on;
plot(t_q, rad2deg(EulerAngs(:,1)), 'LineWidth',2)
xlabel("Time, s")
ylabel("Angle, deg")
title('Yaw')

subplot(3,1,2)
fontsize(14, 'points')
grid on;
plot(t_q, rad2deg(EulerAngs(:,2)), 'LineWidth',2)
xlabel("Time, s")
ylabel("Angle, deg")
title('Pitch')

subplot(3,1,3)
grid on;
plot(t_q, rad2deg(EulerAngs(:,3)), 'LineWidth',2)
xlabel("Time, s")
ylabel("Angle, deg")
fontsize(14, 'points')
title('Roll')
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partD_intermediate_eulerAngs_unstable.png")

%% Arbitrary Axis:

% For equilibrium, w of the satellite should be aligned with w of rotor:
r_rotor = [1,0,0]; % rotor axis of rotation
r_rotor_hat = r_rotor / norm(r_rotor);
w_sat   = [deg2rad(10); 0; deg2rad(10)];
w_r     = 10; %deg2rad(50); % -0.02 rad/s is unstable, deg2rad(50) is stable

% Angular velocity
e       = [-1;0;-1] / sqrt(3);
p       = deg2rad(90); % for PS4-Q1, spacecraft initially aligned with inertial frame;
q_0     = [e(1)*sin(p/2);
           e(2)*sin(p/2);
           e(3)*sin(p/2);
           cos(p/2)];
w_0     = [w_sat; w_r]; % [wx, wy, wz, wr]

% Using the same other initial conditions as first part of this script
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
qw_0    = [q_0; w_0];

[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_DualSpin(qw, M_vec, I_p, Mr, I_rotor, r_rotor_hat), t_span, qw_0, options);

% Check stability:
K1 = (I_p(2,2) - I_p(3,3)) * w_sat(3) / I_rotor;
K2 = (I_p(1,1) - I_p(3,3)) * w_sat(3) / I_rotor;

fprintf('Stability criteria 1: %.4g \n', K1)
fprintf('Stability criteria 2: %.4g \n', K2)
fprintf('Rotor angular velocity: %.4g \n', w_r)

% Confirm the satellite is in equilibrium by looking at the angular
% velocity in inertial frame (components should be constant)

[n, ~] = size(qw_prop);
w_pa  = qw_prop(:,5:7);
w_i_i = zeros(n, 3);
for i = 1:n

    q_i_pa = qw_prop(i,1:4)'; % attitude quaternion for this time step
    A_i_pa = quaternion2dcm(q_i_pa); % rotation from inertial frame to PA
    w_i_i(i, :) = (A_i_pa' * rad2deg(w_pa(i, :))')'; 
end

% Angular velocity over time
figure 
subplot(3,1,1)
fontsize(14, 'points')
grid on;
plot(t_q, w_i_i(:,1), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
title('w_x')

subplot(3,1,2)
fontsize(14, 'points')
grid on;
plot(t_q, w_i_i(:,2), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
title('w_y')

subplot(3,1,3)
grid on;
plot(t_q, w_i_i(:,3), 'LineWidth',2)
xlabel("Time, s")
ylabel("Magnitude, deg/s")
fontsize(14, 'points')
title('w_z')
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partE_arbitrary.png")

% Convert to Euler angles: % Going to use 3-2-1
EulerAngs = zeros(length(qw_prop), 3);

for i = 1:length(qw_prop)
    quat = qw_prop(i,1:4);
    [pitch, yaw, roll] = quat2angle(quat([4 1 2 3]), 'YZX');
    EulerAngs(i,:) = [yaw, pitch, roll];
end

figure()
subplot(3,1,1)
fontsize(14, 'points')
grid on;
plot(t_q, rad2deg(EulerAngs(:,1)), 'LineWidth',2)
xlabel("Time, s")
ylabel("Angle, deg")
title('Yaw')

subplot(3,1,2)
fontsize(14, 'points')
grid on;
plot(t_q, rad2deg(EulerAngs(:,2)), 'LineWidth',2)
xlabel("Time, s")
ylabel("Angle, deg")
title('Pitch')

subplot(3,1,3)
grid on;
plot(t_q, rad2deg(EulerAngs(:,3)), 'LineWidth',2)
xlabel("Time, s")
ylabel("Angle, deg")
fontsize(14, 'points')
title('Roll')
saveas(gcf, "Figures_and_Plots/PS4/DualSpin_partE_arbitrary_eulerAngs.png")