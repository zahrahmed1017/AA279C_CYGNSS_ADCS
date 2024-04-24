%% Load data

close all; clear;

load PropAttitude_Quat_Data.mat
load InertiaData.mat

%% Compute angular velocity and momentum in inertial coordinates

% Get angular momentum , expressed in PA coordinates 
w_i_pa = qw_prop(:,5:7)'; % quantity in inertial frame, but expressed in PA
% coordinates

% This should multiply each column in the 3xN time series array by the matrix, 
% returning another 3xN time series array
L_i_pa = I_p * w_i_pa; 

[n, ~] = size(qw_prop);
L_i_i = zeros(3, n);
w_i_i = zeros(3, n);
for i = 1:n

    q_i_pa = qw_prop(i,1:4)'; % attitude quaternion for this time step
    A_i_pa = quaternion2dcm(q_i_pa); % rotation from inertial frame to PA

    L_i_i(:, i) = A_i_pa' * L_i_pa(:, i); % transpose A because we want to go from PA to inertial
    w_i_i(:, i) = A_i_pa' * w_i_pa(:, i); 
end


%% Plot them

% Angular momentum over time
figure 
hold on
plot(t_q, rad2deg(L_i_i(1,:)))
plot(t_q, rad2deg(L_i_i(2,:)))
plot(t_q, rad2deg(L_i_i(3,:)))
title("Components of angular momentum in inertial coordinates")
xlabel("Time, s")
ylabel("Magnitude, deg/s")
legend("L_1", "L_2", "L_3")
ylim([-3 9])

% Herpelhode trajectory in 3D
figure
hold on
plot3(rad2deg(w_i_i(1,:)), rad2deg(w_i_i(2,:)), rad2deg(w_i_i(3,:)))
view(3)
quiver3( 0, 0, 0,  rad2deg(L_i_i(1,1)), rad2deg(L_i_i(2,2)),  rad2deg(L_i_i(3,3)))
axis equal
grid on


%% Prove that herpelhode is in plane perpendicular to angular momentum

test_vec = [];
for j=1:10
    ind1 = rand_index(n);
    ind2 = rand_index(n);
    plane_vec = w_i_i(:,ind2) - w_i_i(:,ind1); % this should exist in the plane perpendicular to L_i
    test_dot = dot( plane_vec, L_i_i(:,1) );
    test_vec = [test_vec, test_dot];
end



function ind = rand_index(n)
% n -> number of entries in time series

% want a number between 1 and n
ind = round(1 + (n-1).*rand(1,1));
end