close all; clear;

load("InertiaData.mat")
load("OrbitPropData_equil.mat")

%% Scenario 1: 
% Assume that 2 components of the initial angular velocities 
% are zero, and that the principal axes are aligned with the inertial frame
% (e.g., zero Euler angles). Verify that during the simulation the 2 
% components of angular velocity remain zero, and that the attitude 
% represents a pure rotation about the rotation axis (e.g., linearly 
% increasing Euler angle). Plot velocities and angles.


w_0 = [ 0, 0, deg2rad(5)]'; % rad/s (pitch, roll, yaw)
M_vec = [0, 0, 0]';

e = [1;1;1] / sqrt(3);
p = 0; % for PS4-Q1
q_0 = [e(1)*sin(p/2);
       e(2)*sin(p/2);
       e(3)*sin(p/2);
       cos(p/2)];

t_span = 0:3:10*60;

% propagate attitude
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
qw_0 = [q_0; w_0];
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

% these attitudes already represent the rotation from inertial to body

% plot rates - original
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

% test: rotate body rates into inertial
% w_i_pa = qw_prop(:,5:7)'; % quantity in inertial frame, but expressed in PA
% [n, ~] = size(qw_prop);
% w_i_i = zeros(3, n);
% for i = 1:n
% 
%     q_i_pa = qw_prop(i,1:4)'; % attitude quaternion for this time step
%     A_i_pa = quaternion2dcm(q_i_pa); % rotation from inertial frame to PA
%     w_i_i(:, i) = A_i_pa' * w_i_pa(:, i); 
% end
% ^ not needed because it returns the same quantities

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

%% Scenario 2: 
% Repeat a. by setting the initial attitude to match the RTN frame. Set the
% initial angular velocity to be non-zero only about N. Show the evolution 
% of attitude motion in the RTN frame & give an interpretation of the results
