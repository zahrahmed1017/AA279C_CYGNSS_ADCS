close all, clc, clear

save_figs = true;
save_error = false; % Separate condition if we want to save error plots because want to propagate this one for longer

%% CYLINDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

minutes = 5;
tspan = 0:minutes*60; % seconds

% Numerical Solution
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, w_prop_num] = ode113(@(t,w) PropagateEuler(w, M, I_p), tspan, w_init, options);
w_xy_num = sqrt(w_prop_num(:,1).^2 + w_prop_num(:,2).^2);

% Analytical Solution
w_prop_an = Analytical_Euler_TF_AS(w_init,tspan,Ix,Iz);

% Angular Momentum Vector
Lx = w_prop_num(:,1) .* Ix;
Ly = w_prop_num(:,2) .* Iy;
Lz = w_prop_num(:,3) .* Iz;

%% Plots 

% Individual w components 
figure()
fontsize(14,"points")

subplot(3, 1, 1)
fontsize(14,"points")
plot(t, rad2deg(w_prop_num(:,1)),'LineWidth',2)
ylabel("\omega_x, ^\circ/s")
grid on;

subplot(3, 1, 2)
fontsize(14,"points")
plot(t, rad2deg(w_prop_num(:,2)),'LineWidth',2)
ylabel("\omega_y, ^\circ/s")
grid on;

subplot(3, 1, 3)
fontsize(14,"points")
plot(t, rad2deg(w_prop_num(:,3)),'LineWidth',2)
ylabel("\omega_z, ^\circ/s")
xlabel("Time, s")
grid on;

sgtitle("Cylinder Angular Rates")
if save_figs
    saveas(gcf, "Figures_and_Plots/AttitudeTest_CylRates.png")
end

% wz vs. wxy
figure
fontsize(14,"points")
subplot(2,1,1)
fontsize(14,"points")
plot(t, rad2deg(w_prop_num(:,3)),'LineWidth',2)
ylabel("\omega_z")
ylim([min(rad2deg(w_prop_num(:,3))) - 1, min(rad2deg(w_prop_num(:,3))) + 1 ])
grid on;

subplot(2,1,2)
fontsize(14,"points")
plot(t, rad2deg(w_xy_num),'LineWidth',2)
ylabel("\omega_{xy}")
ylim([min(rad2deg(w_xy_num)) - 1, min(rad2deg(w_xy_num)) + 1 ])
xlabel("Time, s")
grid on;

sgtitle("\omega_z and \omega_{xy} Magnitudes")
if save_figs
    saveas(gcf, "Figures_and_Plots/AttitudeTest_CylRateMags.png")
end

% Analytical vs Numerical 
figure()
title("Angular Velocity for Torque-Free, Axially Symmetric Object")
subplot(3,1,1)
fontsize(14,"points")
plot(t, w_prop_num(:,1),'LineWidth',2)
hold on;
grid on;
plot(t, w_prop_an(:,1),'LineWidth',2);
ylabel("\omega_x")
xlabel("Time, s")
legend("Numerical", "Analytical")

subplot(3,1,2)
fontsize(14,"points")
plot(t, w_prop_num(:,2),'LineWidth',2)
hold on;
grid on;
plot(t, w_prop_an(:,2),'LineWidth',2)
ylabel("\omega_y")
xlabel("Time, s")
legend("Numerical", "Analytical")

subplot(3,1,3)
fontsize(14,"points")
plot(t, w_prop_num(:,3),'LineWidth',2)
hold on;
grid on;
plot(t, w_prop_an(:,3),'LineWidth',2)
ylabel("\omega_z")
xlabel("Time, s")
legend("Numerical", "Analytical")

if save_figs
    saveas(gcf, "Figures_and_Plots/AttitudeTest_AnalyticVsNumeric.png")
end

% Analytical vs Numerical Error Plot
figure()
subplot(3,1,1)
plot(t, w_prop_num(:,1) - w_prop_an(:,1), 'LineWidth',2)
grid on;
ylabel('Error in \omega_x')
xlabel('Time, s')

subplot(3,1,2)
plot(t, w_prop_num(:,2) - w_prop_an(:,2), 'LineWidth',2)
grid on;
ylabel('Error in \omega_y')
xlabel('Time, s')

subplot(3,1,3)
plot(t, w_prop_num(:,3) - w_prop_an(:,3), 'LineWidth',2)
grid on;
ylabel('Error in \omega_z')
xlabel('Time, s')

if save_error
    saveas(gcf,"Figures_and_Plots/AttitudeTest_NumericErrors.png")
end

% Angular Velocity and Angular Momentum 
figure()
subplot(1,2,1)
plot3(w_prop_num(:,1), w_prop_num(:,2), w_prop_num(:,3), 'LineWidth', 2)
xlabel('\omega_x')
ylabel('\omega_y')
zlabel('\omega_z')
title('Angular Velocity Vector Tip')
fontsize(14,"points")
grid on;
hold on;

% Plane properties
x = w_prop_num(:,1);
y = w_prop_num(:,2);
z = w_prop_num(:,3);
zPlane = z(1); % z-value of the first point
xRange = [min(x), max(x)]; 
yRange = [min(y), max(y)]; 

% Coordinates for the corners of the plane
X = [xRange(1), xRange(2), xRange(2), xRange(1)];
Y = [yRange(1), yRange(1), yRange(2), yRange(2)];
Z = [zPlane, zPlane, zPlane, zPlane];

% Plot the plane
fill3(X, Y, Z, 'c', 'FaceAlpha', 0.1); % 'c' is cyan, adjust transparency with 'FaceAlpha'

subplot(1,2,2)
plot3(Lx, Ly, Lz, 'LineWidth', 2)
xlabel('L_x')
ylabel('L_y')
zlabel('L_z')
title('Angular Momentum Vector Tip')
fontsize(14,"points")
grid on;
hold on;

% Plane properties
zPlane = Lz(1); % z-value of the first point
xRange = [min(Lx), max(Lx)];
yRange = [min(Ly), max(Ly)]; 

% Coordinates for the corners of the plane
X = [xRange(1), xRange(2), xRange(2), xRange(1)];
Y = [yRange(1), yRange(1), yRange(2), yRange(2)];
Z = [zPlane, zPlane, zPlane, zPlane];

% Plot the plane
fill3(X, Y, Z, 'c', 'FaceAlpha', 0.1); % 'c' is cyan, adjust transparency with 'FaceAlpha'


% if save_figs
%     saveas(gcf, "Figures_and_Plots/AttitudeTest_AngularMomentum.png")
% end



%%%%%%%%%%%%%%%% RECTANGULAR PRISM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rectangular prism w/ width(x) 9 m, depth(y) 4 m, height(z) 1 m, mass 1 kg; z axis is through
% the faces
% width = 9;
% depth = 4;
% height = 1;
% m = 1;
% Ix = (1/12)*m*(height^2 + depth^2);
% Iy = (1/12)*m*(width^2 + height^2);
% Iz = (1/12)*m*(width^2 + depth^2);
% I_p = diag([Ix, Iy, Iz]);
% w_init = [deg2rad(0.0001),  deg2rad(5), deg2rad(0.0001)];
% M = [0,0,0];
% 
% tspan = [0 120*5]; % seconds
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% [t, w_prop] = ode113(@(t,w) PropagateEuler(w, M, I_p), tspan, w_init, options);
% 
% % plot
% figure
% 
% subplot(3, 1, 1)
% plot(t, rad2deg(w_prop(:,1)))
% ylabel("\omega_x, ^\circ/s")
% 
% subplot(3, 1, 2)
% plot(t, rad2deg(w_prop(:,2)))
% ylabel("\omega_y, ^\circ/s")
% 
% subplot(3, 1, 3)
% plot(t, rad2deg(w_prop(:,3)))
% ylabel("\omega_z, ^\circ/s")
% xlabel("Time, s")
% 
% sgtitle("Rec. Prism Angular Rates")
% saveas(gcf, "Figures_and_Plots/AttitudeTest_PrismRates.png")

