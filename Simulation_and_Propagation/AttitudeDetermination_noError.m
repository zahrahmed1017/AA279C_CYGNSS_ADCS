%% Testing attitude determination implementations
close all; clear; clc;

load("InertiaData.mat")
load("cygnss.mat")

%% Attitude estimation with perfect sensors

%%%% INITIALIZE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
h   = 525;
a   = rE + h;
e   = 0.0;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

% Convert to ECI position and velocity for orbit propagation 
oe       = [a; e; i; O; w; v];
rv_state = OE2ECI(oe, muE);
% Initial attitude:
e_vec = [1;1;1]; 
e     = e_vec/ norm(e_vec);
p     = 0; % for PS4-Q1
q_0   = [e(1)*sin(p/2);
         e(2)*sin(p/2);
         e(3)*sin(p/2);
         cos(p/2)];
dcm_0 = quat2dcm(q_0([4 1 2 3])');
% Initial angular velocity:
w_0         = [deg2rad(1), deg2rad(1), deg2rad(2)]';


%%%% Initial Epoch for Testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];
timeStep     = 0; %seconds
fractionDay  = timeStep/86400;
calday       = datetime(initialEpoch(1),initialEpoch(2), initialEpoch(3),0,0,timeStep);
gmst         = CAL2GMST(initialEpoch(1),initialEpoch(2), initialEpoch(3), fractionDay);

T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
numPeriods  = 0.25;
tspan       = 0 : 1 : T * numPeriods; % simulate once an minute?

% propagate orbit (no perturbations for now)
state_0   = [q_0; w_0; rv_state];
gravityGrad = 0;
magTorque   = 0;
aeroDrag    = 0;
SRP         = 0;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t_out, state_out] = ...
    ode113(@(t,state) PropagateOrbit_and_Attitude_wPerturbations(state, I_p, muE, cygnss, t, initialEpoch, gravityGrad, magTorque, aeroDrag, SRP),...
    tspan, state_0, options);

% create some fictitious star unit vectors
star1_i = [1/sqrt(2), 1/sqrt(2), 0]';
star2_i = [0, 1, 0]';

DCM_est_det1 = []; % this is the set of DCMs from the deterministic method w/ 3 unit vecs
DCM_est_det2 = []; % this is the set of DCMs from the deterministic method w/ 2 unit vecs
DCM_est_q = []; % this is the set of DCMs from the q method
q_est_all = [q_0]; % the estimated quaternions from the rate sensor

angs_est1 = [];
angs_est2 = [];
angs_est3 = [];
[yaw_est4, pitch_est4, roll_est4] = quat2angle(q_0([4 1 2 3])', 'ZYX');
angs_est4 = [yaw_est4, pitch_est4, roll_est4];
angs_gt = [];

q_last = q_0; 

for i=1:length(t_out)
    q_now = state_out(i,1:4);
    [yaw_gt, pitch_gt, roll_gt] = quat2angle(q_now([4 1 2 3]), 'ZYX');
    angs_gt = [angs_gt; yaw_gt, pitch_gt, roll_gt];

    % Get fake sensor measurements -- magnetic field vector and sun vector

    [~,B_norm, B_vec_i] = CalculateMagneticTorque(state_out(i,8:10),calday, gmst);
    B_vec_i = B_vec_i/norm(B_vec_i);

    [sun_vec_i, ~] = CalculateSunPositionECI(calday);
    sun_vec_i = sun_vec_i/norm(sun_vec_i);

    % rotate into PA frame at this time step
    R_i_p = quat2dcm(state_out(i,[4 1 2 3]));
    B_vec_p = R_i_p * B_vec_i;
    sun_vec_p = R_i_p * sun_vec_i;
    star1_p =  R_i_p * star1_i;
    star2_p =  R_i_p * star2_i;

    % Deterministic method w/ 3 unit vectors
    % V1 = [B_vec_i, sun_vec_i, star1_i];
    % M1 = [B_vec_p, sun_vec_p, star1_p];
    % R_est1 = M1 * inv(V1);
    R_est1 = deterministicAtt3(B_vec_p, sun_vec_p, star1_p, B_vec_i, sun_vec_i, star1_i);
    DCM_est_det1 = cat(3, DCM_est_det1, R_est1);
    
    % quat_est1 = dcm2quaternion(R_est1);
    % [yaw_est1, pitch_est1, roll_est1] = quat2angle(quat_est1([4 1 2 3])', "ZYX");
    [yaw_est1, pitch_est1, roll_est1] = dcm2angle(R_est1, 'ZYX');
    angs_est1 = [angs_est1; yaw_est1, pitch_est1, roll_est1];
        

    % Deterministic method w/ 2 unit vectors and "dummy"
    % p_p = B_vec_p;
    % q_p = cross(B_vec_p, sun_vec_p) / norm(cross(B_vec_p, sun_vec_p));
    % r_p = cross(p_p, q_p);
    % p_i = B_vec_i;
    % q_i = cross(B_vec_i, sun_vec_i) / norm(cross(B_vec_i, sun_vec_i));
    % r_i = cross(p_i, q_i);
    % V2 = [p_i, q_i, r_i];
    % M2 = [p_p, q_p, r_p];
    % R_est2 = M2 * inv(V2);
    R_est2 = deterministicAtt2(B_vec_p, sun_vec_p, B_vec_i, sun_vec_i);
    DCM_est_det2 = cat(3, DCM_est_det2, R_est2);
    [yaw_est2, pitch_est2, roll_est2] = dcm2angle(R_est2, 'ZYX');
    angs_est2 = [angs_est2; yaw_est2, pitch_est2, roll_est2];
        

    % q-method
    weights = [1 1 2 2];
    weights = weights/norm(weights);
    % W = [sqrt(weights); sqrt(weights); sqrt(weights)] .* [B_vec_p, sun_vec_p, star1_p, star2_p];
    % U = [sqrt(weights); sqrt(weights); sqrt(weights)] .* [B_vec_i, sun_vec_i, star1_i, star2_i];
    % B = W*U';
    % S = B + B';
    % Z = [B(2,3)-B(3,2), B(3,1)-B(1,3), B(1,2)-B(2,1)]';
    % sigma = trace(B);
    % K = [S-eye(3)*sigma, Z; 
    %      Z',       sigma ];
    % [v,d] = eig(K);
    % [~,ind] = max([d(1,1), d(2,2), d(3,3), d(4,4)]);
    % q_est = v(:,ind);
    % R_est3 = quat2dcm(q_est([4 1 2 3])');

    R_est3 = qMethod(weights, [B_vec_p, sun_vec_p, star1_p, star2_p], [B_vec_i, sun_vec_i, star1_i, star2_i]);
    [yaw_est3, pitch_est3, roll_est3] = dcm2angle(R_est3, 'ZYX');
    angs_est3 = [angs_est3; yaw_est3, pitch_est3, roll_est3];
        

    % rate measurement and attitude propagation
    w_meas = state_out(i,5:7); % no error added yet
    q_dot_est = QuaternionKinematics(q_last, w_meas);
    if i < length(t_out)
        dt = t_out(i+1) - t_out(i);
        q_next = q_last + q_dot_est * dt; % integrate as you might onboard
        q_next = q_next / norm(q_next); % renormalize
        q_est_all = [q_est_all q_next];
        q_last = q_next;
        [yaw_est4, pitch_est4, roll_est4] = quat2angle(q_last([4 1 2 3])', 'ZYX');
        angs_est4 = [angs_est4; yaw_est4, pitch_est4, roll_est4];
    end

    

end

%% Plot everything

% deterministic w/ 3
figure
hold on
plot(t_out, rad2deg(angs_est1(:,1)), 'b-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_est1(:,2)), 'r-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_est1(:,3)), 'g-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,1)), 'color', '#9696ff', 'LineStyle', '--', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,2)), 'color', '#fc9696', 'LineStyle', '--', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,3)), 'color', '#96ff96', 'LineStyle', '--', 'LineWidth', 2)
legend("Ground truth yaw", ...
       "Ground truth pitch", ...
       "Ground truth roll", ...
       "Estimated yaw", ...
       "Estimated pitch", ...
       "Estimated roll")
grid on
title("Attitudes estimated w/ deterministic algorithm, 3 references")
xlabel("Time, s")
ylabel("Angle, deg")
saveas(gcf, "Figures_and_Plots/PS6/AttDet_det3.png")

% deterministic w/ 2
figure
hold on
plot(t_out, rad2deg(angs_est2(:,1)), 'b-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_est2(:,2)), 'r-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_est2(:,3)), 'g-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,1)), 'color', '#9696ff', 'LineStyle', '--', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,2)), 'color', '#fc9696', 'LineStyle', '--', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,3)), 'color', '#96ff96', 'LineStyle', '--', 'LineWidth', 2)
legend("Ground truth yaw", ...
       "Ground truth pitch", ...
       "Ground truth roll", ...
       "Estimated yaw", ...
       "Estimated pitch", ...
       "Estimated roll")
grid on
title("Attitudes estimated w/ deterministic algorithm, 2 references")
xlabel("Time, s")
ylabel("Angle, deg")
saveas(gcf, "Figures_and_Plots/PS6/AttDet_det2.png")

% stochastic w/ 4
figure
hold on
plot(t_out, rad2deg(angs_est3(:,1)), 'b-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_est3(:,2)), 'r-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_est3(:,3)), 'g-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,1)), 'color', '#9696ff', 'LineStyle', '--', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,2)), 'color', '#fc9696', 'LineStyle', '--', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,3)), 'color', '#96ff96', 'LineStyle', '--', 'LineWidth', 2)
legend("Ground truth yaw", ...
       "Ground truth pitch", ...
       "Ground truth roll", ...
       "Estimated yaw", ...
       "Estimated pitch", ...
       "Estimated roll")
grid on
title("Attitudes estimated w/ stochastic algorithm, 4 references")
xlabel("Time, s")
ylabel("Angle, deg")
saveas(gcf, "Figures_and_Plots/PS6/AttDet_stoc.png")

% based on rates
figure
hold on
plot(t_out, rad2deg(angs_est4(:,1)), 'b-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_est4(:,2)), 'r-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_est4(:,3)), 'g-', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,1)), 'color', '#9696ff', 'LineStyle', '--', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,2)), 'color', '#fc9696', 'LineStyle', '--', 'LineWidth', 2)
plot(t_out, rad2deg(angs_gt(:,3)), 'color', '#96ff96', 'LineStyle', '--', 'LineWidth', 2)
legend("Ground truth yaw", ...
       "Ground truth pitch", ...
       "Ground truth roll", ...
       "Estimated yaw", ...
       "Estimated pitch", ...
       "Estimated roll")
grid on
title("Attitudes estimated w/ rate measurements")
xlabel("Time, s")
ylabel("Angle, deg")
saveas(gcf, "Figures_and_Plots/PS6/AttDet_rate.png")

