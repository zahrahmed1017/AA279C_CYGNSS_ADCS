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
ang_0 = dcm2eulerAng(quat2dcm(q_0));

% for sanity check: does it convert back the same?
q_test = dcm2quat(eulerAng2dcm(ang_0));

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
plot(t_q, sqrt(qw_prop(:,1).^2 + qw_prop(:,2).^2 + qw_prop(:,3).^2 + qw_prop(:,4).^2))
legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
title("Propagated quaternions")

% check that the two sequences match

%% Propagate the Euler angles

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

ang_w_0 = [ang_0; w_0];

[t_q, ang_w_prop] = ode113(@(t,ang_w) PropagateAttitude_EulerAng(ang_w, M_vec, I_p), t_span, ang_w_0, options);

