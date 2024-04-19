
close all, clc, clear

load("InertiaData.mat")

save_plots = false;
fig_num = 1; % use this when generating plots for many w cases

%% Initial Conditions

w_init = [ 0, deg2rad(5), deg2rad(0.0001)]'; % Hyperbola in XZ frame deg2rad(0.0001)
% w_init = [ 0, deg2rad(5), deg2rad(5)]'; % ellipse in XY plane
% w_init = [ deg2rad(5), deg2rad(0.5), deg2rad(1.5)]'; % ellipse in YZ plane
% w_init = [ deg2rad(5), 0, 0]'; % constant w vector -> overlap is a point


M = [0,0,0];

tspan = [0 1800*5]; % seconds

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

%% Get w trajectory

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, w_prop] = ode113(@(t,w) PropagateAttitude(w, M, I_p), tspan, w_init, options);


%% Plot energy and momentum ellipsoids individually for the mems

% energy ellipsoid
eng_only_fig = figure;
grid on;
hold on;
[X_e, Y_e, Z_e] = ellipsoid(0,0,0, a_e, b_e, c_e, 50);
surf(X_e, Y_e, Z_e, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k')
plot3(w_prop(:,1), w_prop(:,2), w_prop(:,3), 'cyan', 'LineWidth', 2)
axis equal
xlabel("\omega_x")
ylabel("\omega_y")
zlabel("\omega_z")
view(45, 30)
fontsize(14,"points")
hold off;

% momentum ellipsoid
mom_only_fig = figure;
grid on;
hold on;
[X_m, Y_m, Z_m] = ellipsoid(0,0,0, a_m, b_m, c_m, 50);
surf(X_m, Y_m, Z_m, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'k')
plot3(w_prop(:,1), w_prop(:,2), w_prop(:,3), 'cyan', 'LineWidth', 2)
axis equal
xlabel("\omega_x")
ylabel("\omega_y")
zlabel("\omega_z")
view(45, 30)
fontsize(14,"points")
hold off;

if save_plots
    saveas(eng_only_fig, sprintf("Figures_and_Plots/energyEllipsoid_3D_%i.png",fig_num))
    saveas(mom_only_fig, sprintf("Figures_and_Plots/momentumEllipsoid_3D_%i.png",fig_num))
end

%% Plot

ell_fig = figure;
hold on
grid on

% energy ellipsoid
[X_e, Y_e, Z_e] = ellipsoid(0,0,0, a_e, b_e, c_e, 75);
surf(X_e, Y_e, Z_e, 'FaceColor', 'b', 'FaceAlpha', 0.25, 'EdgeColor', 'b', 'EdgeAlpha', '0.75')

% momentum ellipsoid
[X_m, Y_m, Z_m] = ellipsoid(0,0,0, a_m, b_m, c_m, 75);
surf(X_m, Y_m, Z_m, 'FaceColor', 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'r', 'EdgeAlpha', '0.75')

% w trajectory
plot3(w_prop(:,1), w_prop(:,2), w_prop(:,3), 'cyan', 'LineWidth', 2)
plot3(w_prop(1,1), w_prop(1,2), w_prop(1,3), 'Marker', '.', 'MarkerSize', 20, 'Color', 'green', 'LineStyle', 'none')

legend("Energy Ellipsoid", "Momentum Ellipsoid", '\omega Trajectory', 'Initial \omega', 'Location', 'best')

axis equal
xlabel("\omega_x")
ylabel("\omega_y")
zlabel("\omega_z")
fontsize(14,"points")

% set axes to encompass whole ellipsoid(s)
max_wx = max([a_e, a_m]);
max_wy = max([b_e, b_m]);
max_wz = max([c_e, c_m]);
factor = 2; % padding to allow room for legend, etc

xlim([-max_wx*factor, max_wx*factor])
ylim([-max_wy*factor, max_wy*factor])
zlim([-max_wz*factor, max_wz*factor])

if save_plots
    view(45,30)
    saveas(ell_fig, sprintf("Figures_and_Plots/ellipsoids_3D_%i.png",fig_num))

    view(0,90) % XY plane
    saveas(ell_fig, sprintf("Figures_and_Plots/ellipsoids_XY_%i.png",fig_num))

    view(0,0) % XZ plane
    saveas(ell_fig, sprintf("Figures_and_Plots/ellipsoids_XZ_%i.png",fig_num))

    view(90,0) % YZ plane
    saveas(ell_fig, sprintf("Figures_and_Plots/ellipsoids_YZ_%i.png",fig_num))

end


%% Plot angular rates
rate_fig = figure;

subplot(3, 1, 1)
plot(t, rad2deg(w_prop(:,1)), 'LineWidth', 2)
ylabel("\omega_x, ^\circ/s")

subplot(3, 1, 2)
plot(t, rad2deg(w_prop(:,2)), 'LineWidth', 2)
ylabel("\omega_y, ^\circ/s")

subplot(3, 1, 3)
plot(t, rad2deg(w_prop(:,3)), 'LineWidth', 2)
ylabel("\omega_z, ^\circ/s")
xlabel("Time, s")

fontsize(14,"points")
sgtitle("Torque-free Angular Rates")

% Putting all components on one plot to easily show difference in magnitude
rate_single_fig = figure;
grid on;
hold on;
plot(t, rad2deg(w_prop(:,1)), 'LineWidth', 2)
plot(t, rad2deg(w_prop(:,2)), 'LineWidth', 2)
plot(t, rad2deg(w_prop(:,3)), 'LineWidth', 2)
legend("\omega_x, ^\circ/s","\omega_y, ^\circ/s","\omega_z, ^\circ/s")
xlim([0 1000])
xlabel("Time, s")
ylabel("\omega, deg/s")
fontsize(14,"points")
hold off;

if save_plots
    saveas(rate_fig, sprintf("Figures_and_Plots/AngularRates_%i.png",fig_num))
    saveas(rate_single_fig,sprintf("Figures_and_Plots/AngularRates_singleplot_%i.png", fig_num))
end

