close all, clc, clear

% Test attitude propagation script w/ a simple shape
% Short cylinder w/ radius 1 m, height 1 m, mass 1 kg; z axis is through
% the faces
r = 1;
h = 1;
m = 1;
Ix = (1/12)*m*(3*r^2 + h^2);
Iy = (1/12)*m*(3*r^2 + h^2);
Iz = (1/2)*m*r^2;
I_p = diag([Ix, Iy, Iz]);
w_init = [deg2rad(0.5), deg2rad(0.5), deg2rad(5)];
M = [0,0,0];

tspan = [0 120*5]; % seconds
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, w_prop] = ode113(@(t,w) PropagateAttitude(w, M, I_p), tspan, w_init, options);

w_xy = sqrt(w_prop(:,1).^2 + w_prop(:,2).^2);

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

sgtitle("Cylinder Angular Rates")
saveas(gcf, "Figures_and_Plots/AttitudeTest_CylRates.png")

figure
subplot(2,1,1)
plot(t, rad2deg(w_prop(:,3)))
ylabel("\omega_z")
ylim([min(rad2deg(w_prop(:,3))) - 1, min(rad2deg(w_prop(:,3))) + 1 ])

subplot(2,1,2)
plot(t, rad2deg(w_xy))
ylabel("\omega_{xy}")
ylim([min(rad2deg(w_xy)) - 1, min(rad2deg(w_xy)) + 1 ])

sgtitle("\omega_z and \omega_{xy} Magnitudes")
saveas(gcf, "Figures_and_Plots/AttitudeTest_CylRateMags.png")


% rectangular prism w/ width(x) 9 m, depth(y) 4 m, height(z) 1 m, mass 1 kg; z axis is through
% the faces
width = 9;
depth = 4;
height = 1;
m = 1;
Ix = (1/12)*m*(height^2 + depth^2);
Iy = (1/12)*m*(width^2 + height^2);
Iz = (1/12)*m*(width^2 + depth^2);
I_p = diag([Ix, Iy, Iz]);
w_init = [deg2rad(0.0001),  deg2rad(5), deg2rad(0.0001)];
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

sgtitle("Rec. Prism Angular Rates")
saveas(gcf, "Figures_and_Plots/AttitudeTest_PrismRates.png")

