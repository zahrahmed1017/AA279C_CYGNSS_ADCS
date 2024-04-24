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
w_0 = [ deg2rad(5), deg2rad(0.5), deg2rad(1.5)]'; % rad/s
M_vec = [0, 0, 0]';

% initial minute rotation p = 0.001 rad about z axis (to avoid singularity)
e = [1;1;1] / sqrt(3);
p = 0.5;
q_0 = [e(1)*sin(p/2);
       e(2)*sin(p/2);
       e(3)*sin(p/2);
       cos(p/2)];

% TODO write a function for this
% ang_0 = quat2eul([q_0(4); q_0(1:3)]', 'ZXZ'); % matlab quaternion convention puts scalar first

ang_0 = dcm2eulerAng(quaternion2dcm(q_0));

% for sanity check: does it convert back the same?
q_test = dcm2quaternion(eulerAng2dcm(ang_0));

t_span = [0, 10*60];

%% Propagate the quaternions

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

qw_0 = [q_0; w_0];

[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

figure 
hold on
plot(t_q, qw_prop(:,1))
plot(t_q, qw_prop(:,2))
plot(t_q, qw_prop(:,3))
plot(t_q, qw_prop(:,4))
grid on;
plot(t_q, sqrt(qw_prop(:,1).^2 + qw_prop(:,2).^2 + qw_prop(:,3).^2 + qw_prop(:,4).^2), 'LineWidth', 2)
legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
title("Propagated quaternions")

% % check that the two sequences match, convert quat to euler angles
% eulerAngs = zeros(length(qw_prop), 3);
% for i = 1:length(eulerAngs)
%     dcm = quaternion2dcm(qw_prop(i,1:4));
%     % eul = dcm2eulerAng(dcm);
%     [phi, theta, psi] = dcm2angle(dcm,'ZXZ','Robust');
%     eulerAngs(i,:) = [phi, theta, psi];
% end



%% Propagate the Euler angles

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

ang_w_0 = [ang_0; w_0];

[t_q, ang_w_prop] = ode113(@(t,ang_w) PropagateAttitude_EulerAng(ang_w, M_vec, I_p), t_span, ang_w_0, options);

%% Checking Utility functions for converting attitude representations:

% First checking the quaternion2dcm function:
dcm_1 = quaternion2dcm(q_0);
q1 = q_0(1);
q2 = q_0(2);
q3 = q_0(3);
q4 = q_0(4);
dcm_2 = [q4^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4);...
         2*(q1*q2 - q3*q4),          q4^2 - q1^2 + q2^2 - q3^2,  2*(q2*q3 + q1*q4);...
         2*(q1*q3 + q2*q4),          2*(q2*q3 - q1*q4),          q4^2 + q1^2 - q2^2 - q3^2];



% Second check the dcm2eulerAng function:
angs1 = dcm2eulerAng(dcm_1); % Our function
[phi, theta, psi] = dcm2angle(dcm_1, "ZXZ"); % MATLAB function

% Third check the eulerAng2dcm:
dcm_check  = eulerAng2dcm(angs1); % Our function
dcm_check2 = angle2dcm(angs1(1),angs1(2),angs1(3),'ZXZ'); % MATLAB function

% Fourth check dcm2quaternion:
q_check1 = dcm2quaternion(dcm_check);
q_check2 = dcm2quat(dcm_check);