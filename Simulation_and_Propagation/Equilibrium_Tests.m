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

t_span = 0:1:10*60;

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
% saveas(gcf, "Figures_and_Plots/quat_integration_results_EA.png")
saveas(gcf, "Figures_and_Plots/PS4/Q1a_euler_angs.png"); % for PS4

%% Scenario 2: 
% Repeat a. by setting the initial attitude to match the RTN frame. Set the
% initial angular velocity to be non-zero only about N. Show the evolution 
% of attitude motion in the RTN frame & give an interpretation of the results

% Propagate orbit over the time span --> R T N vectors
% Define Orbital Properties
a   = 6903;         % km
e   = 0.00162;
% e = 0;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]
oe = [a; e; i; O; w; v];
% Convert to ECI position and velocity for orbit propagation 
state = OE2ECI(oe, muE);
initial     = state;
options     = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[tout,RVout] = ode113(@(t,State) PropagateOrbit(State, muE),...
                     t_span, initial, options);
RTNout = rv2rtn(RVout);

% get initial rotation from inertial to RTN (1st row from RTNout)
A_eci_rtn = [RTNout(1, 1:3)', RTNout(1, 4:6)', RTNout(1, 7:9)' ]';
% A_eci_rtn

% set this as the initial attitude 
% (rotation from inertial to PA is same as rotation from inertial to RTN)
% q_init = dcm2quaternion(A_eci_rtn); % want to represent an initial rotation from inertial
q_init = dcm2quat(A_eci_rtn);
q_init = q_init([2 3 4 1])';

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% w_0 = [ 0, 0, deg2rad(5)]'; % rad/s about the N axis
w_0 = [ 0, 0, sqrt(muE/(a^3))]';
qw_0 = [q_init; w_0];
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

% plot the attitude motion in the RTN frame
% i.e., the rotation from RTN to PA and body rates
[n,~] = size(qw_prop);
euler_angs_rtn = zeros(n, 3);
body_rates_pa = zeros(n, 3);
body_rates_rtn = zeros(n, 3);

for i=1:n

    % get A_rtn_pa
    A_rtn_eci =  [RTNout(i, 1:3)', RTNout(i, 4:6)', RTNout(i, 7:9)' ];
    % A_eci_pa = quaternion2dcm(qw_prop(i,1:4)); % something wring with
    % this function 
    quat_i = qw_prop(i,1:4);
    A_eci_pa = quat2dcm(quat_i([4, 1, 2, 3]) );
    A_rtn_pa = A_eci_pa * A_rtn_eci;

    % convert to quaternion
    q_i = dcm2quaternion(A_rtn_pa); % we actually don't need this

    % save euler angles
    [phi_i, theta_i, psi_i] = dcm2angle(A_rtn_pa, 'ZYX');
    angs_i = [phi_i, theta_i, psi_i] ;
    euler_angs_rtn(i,:) = angs_i;

    % save body rates in PA
    body_rates_pa(i,:) = qw_prop(i,5:7);

    % save body rates in RTN
    body_rates_rtn(i,:) = A_rtn_pa' * qw_prop(i,5:7)';


end

figure 
hold on
plot(t_q, rad2deg(euler_angs_rtn(:,1)),'LineWidth', 2);
plot(t_q, rad2deg(euler_angs_rtn(:,2)), 'LineWidth', 2);
plot(t_q, rad2deg(euler_angs_rtn(:,3)), 'LineWidth', 2);
legend('\phi', '\theta', '\psi')
title('Quaternion Propagation shown as ZYX Euler Angles Attitude in RTN frame')
grid on; 
xlabel("Time, s")
ylabel("Angle, ^\circ")
saveas(gcf, "Figures_and_Plots/PS4/Q1b_euler_angs.png")

figure 
hold on
plot(t_q, rad2deg(body_rates_rtn(:,1)),'LineWidth', 2);
plot(t_q, rad2deg(body_rates_rtn(:,2)), 'LineWidth', 2);
plot(t_q, rad2deg(body_rates_rtn(:,3)), 'LineWidth', 2);
legend('\omega_x', '\omega_y', '\omega_z')
title('Body rates in RTN frame')
grid on; 
ylim([-0.5, 5.5])
xlabel("Time, s")
ylabel("Rate, ^\circ /s")
saveas(gcf, "Figures_and_Plots/PS4/Q1b_ang_vel_rtn.png")


figure 
hold on
plot(t_q, rad2deg(body_rates_pa(:,1)),'LineWidth', 2);
plot(t_q, rad2deg(body_rates_pa(:,2)), 'LineWidth', 2);
plot(t_q, rad2deg(body_rates_pa(:,3)), 'LineWidth', 2);
legend('\omega_x', '\omega_y', '\omega_z')
title('Body rates in PA frame')
grid on; 
ylim([-0.5, 5.5])
xlabel("Time, s")
ylabel("Rate, ^\circ /s")
saveas(gcf, "Figures_and_Plots/PS4/Q1b_ang_vel_pa.png")

