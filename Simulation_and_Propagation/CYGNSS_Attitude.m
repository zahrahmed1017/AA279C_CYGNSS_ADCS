close all, clc, clear

load("InertiaData.mat")

% Initial conditions

w_init = [ 0, deg2rad(5), deg2rad(0.0001)]; % PA frame
M = [0,0,0];

tspan = [0 120*5]; % seconds
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, w_prop] = ode113(@(t,w) PropagateAttitude(w, M, I_p), tspan, w_init, options);

% plot
figure

subplot(3, 1, 1)
plot(t, rad2deg(w_prop(:,1)))
ylabel("\omega_x, ^\circ/s")

subplot(3, 1, 2)
plot(t, rad2deg(w_prop(:,2)))
ylabel("\omega_y, ^\circ/s")

subplot(3, 1, 3)
plot(t, rad2deg(w_prop(:,3)))
ylabel("\omega_z, ^\circ/s")
xlabel("Time, s")

sgtitle("Torque-free Angular Rates")
