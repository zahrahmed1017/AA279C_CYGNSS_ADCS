
close all, clc, clear

load("InertiaData.mat")

%% Initial Conditions

w_init = [ 0, deg2rad(5), deg2rad(0.0001)]'; % PA frame
M = [0,0,0];

%% Energy Ellipsoid

T = RotationalKineticEnergy(I_p, w_init);

a_e = sqrt(2*T / I_p(1,1));
b_e = sqrt(2*T / I_p(2,2));
c_e = sqrt(2*T / I_p(3,3));

%% Momentum Ellipsoid

L_vec = AngularMomentum(I_p, w_init);
L = sqrt(dot(L_vec,L_vec));

a_m = L / I_p(1,1);
b_m = L / I_p(2,2);
c_m = L / I_p(3,3);


%% Plot

figure
hold on

% energy ellipsoid
[X_e, Y_e, Z_e] = ellipsoid(0,0,0, a_e, b_e, c_e, 75);
surf(X_e, Y_e, Z_e, 'FaceColor', 'b', 'FaceAlpha', 0.5)

% momentum ellipsoid
[X_m, Y_m, Z_m] = ellipsoid(0,0,0, a_m, b_m, c_m, 75);
surf(X_m, Y_m, Z_m, 'FaceColor', 'r', 'FaceAlpha', 0.5)

axis equal
ylim([-.2,.2])
xlim([-.2,.2])
zlim([-.2,.2])