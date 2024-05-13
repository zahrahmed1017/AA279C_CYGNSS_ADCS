%% PS5 - Q1 - Gravity Gradient Stability

close all; clear;

load("InertiaData.mat")

%% Part a - Calculate k coefficients; plot stable/unstable regions

% create the stability plot - shade UNSTABLE regions
figure
hold on

alpha_patch = 0.2;

% K_R >= K_T 
p1_x = [-1, 1, -1, -1];
p1_y = [-1, 1, 1, -1];
p1 = patch(p1_x, p1_y, 'y');
p1.FaceVertexAlphaData = alpha_patch;
p1.FaceAlpha = 'flat';

% K_R * K_T <= 0
p2_x = [-1, -1, 0, 0, 1, 1, -1];
p2_y = [0, 1, 1, -1, -1, 0, 0];
p2 = patch(p2_x, p2_y, 'b');
p2.FaceVertexAlphaData = alpha_patch;
p2.FaceAlpha = 'flat';

% 1 + 3*K_T + K_R*K_T <= 4*sqrt(K_R*K_T)
p3_x = linspace(-1, -0.0505, 100);
p3_y = (-3* p3_x.^2 + 4*sqrt(3)*sqrt(p3_x.^2 - p3_x.^3) + 7*p3_x) ./ ...
    (p3_x.^2);
p3_x = [p3_x, -1, -1];
p3_y = [p3_y, -1, 0];
p3 = patch(p3_x, p3_y, 'r');
p3.FaceVertexAlphaData = alpha_patch;
p3.FaceAlpha = 'flat';

axis equal
ylim([-1, 1])
xlim([-1, 1])
xlabel("$K_T$", 'Interpreter', 'latex')
ylabel("$K_R$", 'Interpreter', 'latex')



% First K coefficients: z aligned w/ N

Ix = I_p(1,1);
Iy = I_p(2,2);
Iz = I_p(3,3);

K_N1 = (Iy - Ix) / Iz;
K_R1 = (Iz - Iy) / Ix;
K_T1 = (Iz - Ix) / Iy;

plot(K_T1, K_R1, 'kx')

% Second K coefficients: x aligned w/ T, Y aligned w/ N

K_N2 = (Iz - Ix) / Iy;
K_R2 = (Iy - Ix) / Iz;
K_T2 = (Iy - Iz) / Ix;

plot(K_T2, K_R2, 'k^')

% Third K coefficients: X aligned w/ N; Z aligned w/ T

K_N3 = (Iz - Iy) / Ix;
K_R3 = (Ix - Iz) / Iy;
K_T3 = (Ix - Iy) / Iz;
plot(K_T3, K_R3, 'ksquare')

% 4th K coefficients: adjust the inertia parameters

Ix2 = 1.2;
Iy2 = 1.3;
Iz2 = 1.5;
I_p2 = [Ix2, 0, 0; 0, Iy2, 0; 0, 0, Iz2];

K_N4 = (Iz2 - Iy2) / Ix2;
K_R4 = (Ix2 - Iz2) / Iy2;
K_T4 = (Ix2 - Iy2) / Iz2;
plot(K_T4, K_R4, 'k*')

box on

legend("$K_R \geq K_T$", "$K_R K_T \leq 0$", "$1 + 3 K_T + K_R K_T \leq 4\sqrt{K_R*K_T}$", ...
     "$Z_p \parallel N$; $X_p \parallel R$",...
     "$Y_p \parallel N$; $Z_p \parallel -R$", ...
     "$X_p \parallel N$; $Y_p \parallel R$", ...
     "$X_p \parallel N$; $Y_p \parallel R$; new inertia tensor", ...
     'Interpreter', 'latex', 'location', 'eastoutside')

saveas(gcf, "Figures_and_Plots/PS5/Q1a_StabilityPlot.png");


%% Part b - Try to reproduce stable & unstable motion - Z aligned w/ N case

% load("OrbitPropData_equil.mat")

scenario = 4;
% 1 -> 1st K coeffs
% 2 -> 2nd K coeffs
% 3 -> 3rd K coeffs
% 4 -> 4th K coeffs
% note: only implemented 1, 2, 4

rE  = 6378; % km
% a   = rE + 3000;     % km
a   = rE + 525;
e   = 0.00162;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

% Convert to ECI position and velocity for orbit propagation 
oe       = [a; e; i; O; w; v];
rv_state = OE2ECI(oe, muE);

%%%% SET SIMULATION TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
numPeriods  = 40;
% tspan       = 0 : 1 : T * numPeriods; % simulate once an minute?
% tspan       = 0 : 10 : 100;
tspan      = [0, T*numPeriods];

%%%% INITIALIZE ATTITUDE AND ANGULAR RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For altered inertia tensor:
if scenario == 4
    I_p = I_p2;
end

% Initial angular velocity:
n           = sqrt(muE / a^3);
pert = 0.01*n;

switch scenario
    case 1
        % rotate about Z
        w_0         = [pert, pert, n + pert]';
        % % Initial attitude: XYZ aligned w/ RTN
        rtn       = rv2rtn(rv_state');
        R_eci_rtn = [rtn(1:3); rtn(4:6); rtn(7:9)];
        q_0       = dcm2quat(R_eci_rtn);
        q_0       = q_0([2 3 4 1])';

    case 2
        % rotate about Y
        w_0         = [pert, n + pert, pert]';

        
        % Initial attitude: Y aligned w/ N; Z aligned w/ R
        rtn       = rv2rtn(rv_state');
        R_eci_rtn = [rtn(1:3); rtn(4:6); rtn(7:9)];
        R_rtn_pa  = [0, 1, 0; 
                     0, 0, 1;
                     1, 0, 0];
        q_0       = dcm2quat(R_rtn_pa * R_eci_rtn);
        q_0       = q_0([2 3 4 1])';


    case {3,4}
        % rotate about X
        w_0         = [n + pert, pert,  pert]';

        % Initial attitude: X aligned w/ N; Z aligned w/ T
        rtn       = rv2rtn(rv_state');
        R_eci_rtn = [rtn(1:3); rtn(4:6); rtn(7:9)];
        R_rtn_pa  = [0, 0, 1; 
                     1, 0, 0;
                     0, -1, 0];
        q_0       = dcm2quat(R_rtn_pa * R_eci_rtn);
        q_0       = q_0([2 3 4 1])';

end


state_0   = [q_0; w_0; rv_state];

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);

[t_q, state_out] = ode113(@(t,state) PropagateOrbit_and_Attitude_GravityGrad(state, I_p, muE), tspan, state_0, options);

% Plot rates
figure
hold on
plot(t_q, rad2deg(state_out(:,5)), 'LineWidth', 2)
plot(t_q, rad2deg(state_out(:,6)), 'LineWidth', 2)
plot(t_q, rad2deg(state_out(:,7)), 'LineWidth', 2)
xlabel( 'Time, s')
ylabel("Body rate, ^\circ/s")
title("Body rates in PA frame")
legend('\omega_x', '\omega_y', '\omega_z')


switch scenario
    case 1
        saveas(gcf, "Figures_and_Plots/PS5/Q1b_PArates_case1.png");
    case 2
        saveas(gcf, "Figures_and_Plots/PS5/Q1b_PArates_case2.png");
    case 3
        saveas(gcf, "Figures_and_Plots/PS5/Q1b_PArates_case3.png");
    case 4
        saveas(gcf, "Figures_and_Plots/PS5/Q1b_PArates_case4.png");
end

% plot angles

% get euler angles
eulerAngs = zeros(length(state_out), 3);
for i = 1:length(eulerAngs)
    quat = state_out(i,1:4);
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
title('ZYX Euler Angles')
grid on; 
xlabel("Time, s")
ylabel("Angle, ^\circ")

switch scenario
    case 1
        saveas(gcf, "Figures_and_Plots/PS5/Q1b_EulerAngles_case1.png");
    case 2
        saveas(gcf, "Figures_and_Plots/PS5/Q1b_EulerAngles_case2.png");
    case 3
        saveas(gcf, "Figures_and_Plots/PS5/Q1b_EulerAngles_case3.png");
    case 4
        saveas(gcf, "Figures_and_Plots/PS5/Q1b_EulerAngles_case4.png");
end








