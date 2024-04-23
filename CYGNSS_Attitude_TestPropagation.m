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
w_0 = [0.001, 0, 5]'; % rad/s
M_vec = [0, 0, 0]';

% initial minute rotation p = 0.001 rad about z axis (to avoid singularity)
e = [0;0;1];
p = 0.001;
q_0 = [e(1)*sin(p/2);
       e(2)*sin(p/2);
       e(3)*sin(p/2);
       cos(p/2)];

% TODO write a function for this
% ang_0 = quat2eul([q_0(4); q_0(1:3)]', 'ZXZ'); % matlab quaternion convention puts scalar first

t_span = [0, 10*60];

%% Propagate the quaternions

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

qw_0 = [q_0; w_0];

[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), t_span, qw_0, options);

figure 
hold on
plot(t_q, qw_prop(1))
plot(t_q, qw_prop(2))
plot(t_q, qw_prop(3))
% legend("q_1", "q_2", "q_3", "q_4")
title("Propagated quaternions")

% check that the two sequences match


