close all; clear;

load("InertiaData.mat")

%% Scenario 1: Assume that 2 components of the initial angular velocities 
% are zero, and that the principal axes are aligned with the inertial frame
% (e.g., zero Euler angles). Verify that during the simulation the 2 
% components of angular velocity remain zero, and that the attitude 
% represents a pure rotation about the rotation axis (e.g., linearly 
% increasing Euler angle). Plot velocities and angles.


w_0 = [ 0, 0, deg2rad(5)]'; % rad/s (pitch, roll, yaw)
M_vec = [0, 0, 0]';

e = [1;1;1] / sqrt(3);
q_0 = [e(1)*sin(p/2);
       e(2)*sin(p/2);
       e(3)*sin(p/2);
       cos(p/2)];

t_span = 0:3:10*60;

% propagate attitude
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
qw_0 = [q_0; w_0];
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

% rotate to inertial frame



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


