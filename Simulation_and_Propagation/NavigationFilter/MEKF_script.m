close all;
clear;
load("InertiaData.mat")
load("cygnss.mat")
load("measurementCov.mat")

%% Ground Truth Orbit/Attitude Propagation for Sensor Measurements

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
w_0         = [deg2rad(1), deg2rad(0), deg2rad(1)]';


%%%% Initial Epoch for Testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];
timeStep     = 0; %seconds
fractionDay  = timeStep/86400;
calday       = datetime(initialEpoch(1),initialEpoch(2), initialEpoch(3),0,0,timeStep);
gmst         = CAL2GMST(initialEpoch(1),initialEpoch(2), initialEpoch(3), fractionDay);

T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
numPeriods  = 0.1;
dt          = 2; %seconds
tspan       = 0 : dt : T * numPeriods; % 

% propagate orbit (no perturbations for now)
state_0     = [q_0; w_0; rv_state];
gravityGrad = 1;
magTorque   = 1;
aeroDrag    = 1;
SRP         = 1;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t_out, state_out] = ...
    ode113(@(t,state) PropagateOrbit_and_Attitude_wPerturbations(state, I_p, muE, cygnss, t, initialEpoch, gravityGrad, magTorque, aeroDrag, SRP),...
    tspan, state_0, options);


%% MEKF Initialization
% q_0 = [0; 0; 0; 1]; % aligned with inertial, scalar last
% w_0 = [0; deg2rad(5); 0];

% Initial State and Covariance
state_0        = [q_0; w_0];

cov_0          = [1e-3, 0,    0,    0,    0,    0;...
                   0,    1e-3, 0,    0,    0,    0;...
                   0,    0,    1e-3, 0,    0,    0;...
                   0,    0,    0,    1e-12,0,   0;...
                   0,    0,    0,    0,    1e-12,0;...
                   0,    0,    0,    0,    0,     1e-12];
cov_0 = cov_0 * 50;


% cov_0 = eye(6,6) * 0.1;

%% Adding noise to initial state estimate (rate and ref. quaternion)

state_0(5:7)        = state_0(5:7) + cov_0(4:6, 4:6) * randn(3,1);

init_err_mrp        = cov_0(1:3, 1:3) * randn(3,1); 
init_ref_quat       = quatMul(state_0(1:4), mrp2quat(init_err_mrp));
state_0(1:4)        = init_ref_quat;


% Initial Process Noise 
% processNoise_0 = cov_0./100;
% processNoise_0 = cov_0/100;
processNoise_0 = [1e-3, 0,    0,    0,    0,    0;...
                   0,    1e-3, 0,    0,    0,    0;...
                   0,    0,    1e-3, 0,    0,    0;...
                   0,    0,    0,    1e-12,0,   0;...
                   0,    0,    0,    0,    1e-12,0;...
                   0,    0,    0,    0,    0,     1e-12] ;



% Initialize Sun Sensor
ss_ang_noise   = deg2rad(0.1); % 1 sigma
% ss_ang_noise   = deg2rad(1); % 1 sigma
ss_ang_quant   = deg2rad(0.125);
ss             = SunSensor(ss_ang_noise, ss_ang_quant, calday);
% For debug: use a basic sensor model
ss_noise_simple = deg2rad(1);

% Initialize Magnetometer
mag_v_bias     = [0, 0, 0]'; % init as zero because we assume good calibration
mag_v_noise    = 0.025; % volts
mag_F_matrix   = 1e5 * eye(3); % V/T
mag            = MagSensor(mag_v_bias, mag_v_noise, mag_F_matrix, calday, gmst);
mag_noise_simple = deg2rad(1);

% Initialize measurement covariance

% try: tweaking the measurement covariance (esp. for magnetometer)
% not really helping so far
meas_cov = eye(6);
meas_cov(1:3, 1:3) = meas_cov(1:3, 1:3) * ss_noise_simple;
meas_cov(4:6, 4:6) = meas_cov(4:6, 4:6) * mag_noise_simple;

measCov = meas_cov;

% measCov        = [ss_ang_noise^2, 0,          0,          0,               0,               0;...
%                   0,          ss_ang_noise^2, 0,          0,               0,               0;...
%                   0,          0,          ss_ang_noise^2, 0,               0,               0;...
%                   0,          0,          0,              mag_v_noise^2,   0,               0;...
%                   0,          0,          0,              0,               mag_v_noise^2,   0;...
%                   0,          0,          0,              0,               0,               mag_v_noise^2];
% dt = 5; %seconds
% tspan = 0:dt:1000;
% filter = MEKF(state_0, cov_0, processNoise_0, dt, I_p);


%% Run Filter 

q_mekf      = zeros(length(tspan), 4);
w_mekf      = zeros(length(tspan), 3);
mrp_mekf    = zeros(length(tspan), 3);
meas_mekf   = zeros(length(tspan), 6);
model_mekf  = zeros(length(tspan), 6);
prefit_err  = zeros(length(tspan), 6);
postfit_err = zeros(length(tspan), 6);
cov_mekf    = zeros(length(tspan), 3);
kgain_mekf  = zeros(length(tspan), 1);

if size(state_0,2) ~= 1
    state_0 = state_0';
end

obj.state         = state_0;
obj.errorState    = [0; 0; 0; state_0(5:7)];
obj.stateCov      = cov_0;
obj.processNoise  = processNoise_0; 
obj.measCov       = measCov;
obj.dt            = dt;
obj.I             = I_p;

for i = 1:length(tspan)-1

    %%%%%%%%%%% STEP 1 filter.Propagate() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 1.1 QuatKin_Linear_Discrete
    w = obj.errorState(4:6);
    w1 = w(1);
    w2 = w(2);
    w3 = w(3);
    skewSym     = [0,  w3, -w2, w1;...
                   -w3, 0,   w1, w2;...
                    w2, -w1, 0,  w3;...
                   -w1, -w2, -w3, 0];
    obj.quatKin = eye(4,4) + (0.5).* skewSym .* obj.dt;

    %%%% 1.2 PropagateQuat
    Omega = obj.quatKin;
    qt_1  = obj.state(1:4);
    qt    = Omega * qt_1;
    obj.state(1:4) = qt/norm(qt);  

    obj.errorState(1:3) = [0; 0; 0];

   %%%%%%%%%%%%% STEP 2 filter.TimeUpdate() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%% 2.1. Update STM

   %%%% Update A11
    if size(obj.errorState,2) ~= 1
        obj.errorState = obj.errorState';
    end

    wx = crossMatrix(obj.errorState(4:6));
    obj.A11 = eye(3,3) - (1/2) * obj.dt .* wx ; 

   %%%% Update A12
    obj.A12 = (1/4) * obj.dt .* eye(3,3);

   %%%% Update A22
    Ix = obj.I(1,1);
    Iy = obj.I(2,2);
    Iz = obj.I(3,3);
    wx = obj.errorState(4);
    wy = obj.errorState(5);
    wz = obj.errorState(6);

    eulEqLin = [ 0,                 (Iy - Iz)/Ix * wz, (Iy - Iz)/Ix * wy;...
                (Iz - Ix)/Iy * wz,   0,                (Iz - Ix)/Iy * wx;...
                (Ix - Iy)/Iz * wy,   (Ix - Iy)/Iz * wx, 0];
    obj.A22  = eye(3,3) + (0.5* obj.dt .* eulEqLin);

    %%%% Update STM
    obj.STM = [obj.A11, obj.A12; ...
               zeros(3,3), obj.A22];

    %%%% Update Angular Velocity
    obj.errorState(4:6) = obj.A22 * obj.errorState(4:6);
    w_mekf(i,:)         = obj.errorState(4:6)';
            
    %%%% Update State Covariance
    obj.stateCov = obj.STM * obj.stateCov * obj.STM' + obj.processNoise;

    cov_mekf(i,1) = sqrt(obj.stateCov(1,1));
    cov_mekf(i,2) = sqrt(obj.stateCov(2,2));
    cov_mekf(i,3) = sqrt(obj.stateCov(3,3));

   % %%%%%%%%%%%%% STEP 3: filter.MeasurementUpdate() %%%%%%%%%%%%%%%%%%%%%%%

   %%% Get true measurements - 3n x 1 vector where n is the number of
   %%% measurements
   R_i_p    = quat2dcm(state_out(i+1,[4 1 2 3])); % Ground truth attitude

   % % % commenting these out for now, doing simple model instead
   % sun_meas = ss.get_measurement(R_i_p); 
   % mag_meas = mag.get_measurement(R_i_p, state_out(i+1,8:10));
   % y        = [sun_meas; mag_meas];
   % meas_mekf(i,:) = y';

   %%% Get modeled measurements - 3n x 1 vector where n is the number of
   %%% measurements
   quat_est = obj.state(1:4)';
   R_est    = quat2dcm(quat_est([4 1 2 3]));
   % R_est = quat2dcm(state_out(i+1, [4 1 2 3]));

   [sun_vec_i, ~]       = CalculateSunPositionECI(calday);
   sun_vec_i            = sun_vec_i/norm(sun_vec_i);
   [~, B_norm, B_vec_i] = CalculateMagneticTorque(state_out(i+1,8:10),calday, gmst);
   B_vec_i              = B_vec_i/norm(B_vec_i);

   % for SIMPLE measurement simulation
   sun_meas = (-crossMatrix(randn(3,1)*ss_noise_simple) + eye(3)) * (R_i_p *sun_vec_i) ;
   sun_meas = sun_meas / norm(sun_meas);
   mag_meas = (-crossMatrix(randn(3,1)*mag_noise_simple) + eye(3)) * (R_i_p * B_vec_i) ;
   mag_meas = mag_meas / norm(mag_meas);
   y        = [sun_meas; mag_meas];
   meas_mekf(i,:) = y';

   z1 = R_est * sun_vec_i; % Sun sensor modeled measurement
   z2 = R_est * B_vec_i; % Magnetometer modeled measurement
   z  = [z1; z2];
   model_mekf(i,:) = z';

   %%% PRE-FIT RESIDUALS %%%%%%%%%
   prefit_err(i,:) = y' - z';

   %%% Calculate Sensitivity Matrix - 3n x 6 
   z1x = crossMatrix(z1);
   z2x = crossMatrix(z2);
   H   = zeros(length(z), length(obj.errorState));
   H(1:3, 1:3) = z1x;
   H(4:6, 1:3) = z2x;

   %%% Calculate Kalman Gain 
   K_t = obj.stateCov * H' * ((H * obj.stateCov * H') + obj.measCov)^(-1);

   %%% Update state
   obj.errorState = obj.errorState + K_t * (y - z);
   mrp_mekf(i,:)  = obj.errorState(1:3)';

   %%% Update covariance 
   obj.stateCov   = ((eye(6,6) - K_t * H) * obj.stateCov * (eye(6,6) - K_t * H)') + ...
       (K_t * obj.measCov * K_t');

    % cov_mekf(i,1) = sqrt(obj.stateCov(1,1));
    % cov_mekf(i,2) = sqrt(obj.stateCov(2,2));
    % cov_mekf(i,3) = sqrt(obj.stateCov(3,3));

   %%%%%%%%%%%% STEP 4: filter.Reset() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Push error state to absolute state
    % MRP -> Quat
    % angular velocity in angular state transfers to absolute state
    q_old             = obj.state(1:4); % scalar last here
    errState_quat     = mrp2quat(obj.errorState(1:3));
    % errState_quat     = [obj.errorState(1:3)/2; 1];
    % q_new             = quatMul(errState_quat, q_old);
    q_new             = quatMul(q_old, errState_quat);
    obj.state(1:4)    = q_new / norm(q_new);
    obj.state(5:7)    = obj.errorState(4:6);
    q_mekf(i, :)      = obj.state(1:4)';  

   %%%%%%%%%%%%% POST-FIT RESIDUAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   quat_est = obj.state(1:4)';
   R_est    = quat2dcm(quat_est([4 1 2 3]));

   [sun_vec_i, ~]      = CalculateSunPositionECI(calday);
   sun_vec_i           = sun_vec_i/norm(sun_vec_i);
   [~,B_norm, B_vec_i] = CalculateMagneticTorque(state_out(i+1,8:10),calday, gmst);
   B_vec_i             = B_vec_i/norm(B_vec_i);

   z1_post = R_est * sun_vec_i; % Sun sensor modeled measurement
   z2_post = R_est * B_vec_i;   % Magnetometer modeled measurement
   z_post  = [z1_post; z2_post];

   postfit_err(i,:) = y' - z_post';

end



%% Compare to numerical propagation (torque-free)

qw_0  = [q_0; w_0];
M_vec = [0; 0; 0];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_q, qw_prop] = ode113(@(t,qw) PropagateAttitude_Quat(qw, M_vec, I_p), tspan, qw_0, options);

figure 
hold on
plot(t_q, qw_prop(:,1), 'LineWidth', 2)
plot(t_q, qw_prop(:,2), 'LineWidth', 2)
plot(t_q, qw_prop(:,3), 'LineWidth', 2)
plot(t_q, qw_prop(:,4), 'LineWidth', 2)
grid on;
plot(t_q, sqrt(qw_prop(:,1).^2 + qw_prop(:,2).^2 + qw_prop(:,3).^2 + qw_prop(:,4).^2), 'LineWidth', 2)
legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
title("Numerically Propagated quaternions")
fontsize(15,'points')
xlabel('Time [s]')
% saveas(gcf, "Figures_and_Plots/quat_integration_results_q.png")

figure 
hold on
plot(t_q, q_mekf(:,1),  'LineWidth', 2)
plot(t_q, q_mekf(:,2),  'LineWidth', 2)
plot(t_q, q_mekf(:,3),  'LineWidth', 2)
plot(t_q, q_mekf(:,4),  'LineWidth', 2)
grid on;
plot(t_q, sqrt(q_mekf(:,1).^2 + q_mekf(:,2).^2 + q_mekf(:,3).^2 + q_mekf(:,4).^2), 'LineWidth', 2)
legend("q_1", "q_2", "q_3", "q_4", "q_{mag}")
ylim([-1 1])
fontsize(15,'points')
title("MEKF quaternions")
xlabel('Time [s]')
% saveas(gcf, "Figures_and_Plots/quat_integration_results_q.png")


figure()
title("Quaternion Errors")
subplot(4,1,1)
plot(t_q, q_mekf(:,1) - qw_prop(:,1),  'LineWidth', 2)
ylabel('\Delta q_1')
fontsize(15,'points')
grid on;
subplot(4,1,2)
plot(t_q, q_mekf(:,2) - qw_prop(:,2),  'LineWidth', 2)
ylabel('\Delta q_2')
fontsize(15,'points')
grid on;
subplot(4,1,3)
plot(t_q, q_mekf(:,3) - qw_prop(:,3),  'LineWidth', 2)
ylabel('\Delta q_3')
fontsize(15,'points')
grid on;
subplot(4,1,4)
plot(t_q, q_mekf(:,4) - qw_prop(:,4),  'LineWidth', 2)
ylabel('\Delta q_4')
fontsize(15,'points')
grid on;
xlabel('Time [s]')
% plot(t_q, sqrt(q_mekf(:,1).^2 + q_mekf(:,2).^2 + q_mekf(:,3).^2 + q_mekf(:,4).^2), 'LineWidth', 2)
% legend("q_1", "q2", "q_3", "q_4", "q_{mag}")
% title("MEKF quaternions")

figure()
subplot(3,1,1)
plot(t_q, w_mekf(:,1), 'LineWidth', 2)
hold on;
plot(t_q, w_mekf(:,2), 'LineWidth', 2)
plot(t_q, w_mekf(:,3), 'LineWidth', 2)
title('MEKF angular velocity')
legend('wx', 'wy', 'wz')
fontsize(14,'points')
subplot(3,1,2)
plot(t_q, qw_prop(:,5), 'LineWidth', 2)
hold on
plot(t_q, qw_prop(:,6), 'LineWidth', 2)
plot(t_q, qw_prop(:,7), 'LineWidth', 2)
title('Numerical prop angular velocity')
legend('wx', 'wy', 'wz')
fontsize(14,'points')
subplot(3,1,3)
plot(t_q, w_mekf(:,1) - qw_prop(:,5),  'LineWidth', 2)
hold on;
plot(t_q, w_mekf(:,2) - qw_prop(:,6),  'LineWidth', 2)
plot(t_q, w_mekf(:,3) - qw_prop(:,7),  'LineWidth', 2)
title('Error angular velocity')
legend('wx', 'wy', 'wz')
fontsize(14,'points')

%% Plot covariance errors
n = length(tspan);

figure()
subplot(3,1,1)
plot(tspan, cov_mekf(:,1), 'LineWidth', 2)
hold on;
grid on;
plot(tspan, -cov_mekf(:,1),  'LineWidth', 2)
plot(tspan, mrp_mekf(:,1), 'LineWidth', 2)
legend('+\sigma', '-\sigma', 'p1')
xlabel('Time [s]')
ylabel('\sigma_{P1}')
fontsize(15,'points')

subplot(3,1,2)
plot(tspan(1:n), cov_mekf(1:n,2), 'LineWidth', 2)
hold on;
grid on;
plot(tspan(1:n), -cov_mekf(1:n,2), '-r', 'LineWidth', 2)
plot(tspan(1:n), mrp_mekf(1:n,2), 'LineWidth', 2)
legend('+\sigma', '-\sigma', 'p2')
xlabel('Time [s]')
ylabel('\sigma_{P2}')
fontsize(15,'points')

subplot(3,1,3)
plot(tspan(1:n), cov_mekf(1:n,3), 'LineWidth', 2)
hold on;
grid on;
plot(tspan(1:n), -cov_mekf(1:n,3), 'LineWidth', 2)
plot(tspan(1:n), mrp_mekf(1:n,3), 'LineWidth', 2)
legend('+\sigma', '-\sigma', 'p3')
xlabel('Time [s]')
ylabel('\sigma_{P3}')
fontsize(15,'points')

%% Plot pre-fit and post-fit residuals:

figure()
subplot(3,2,1)
plot(tspan, abs(prefit_err(:,1)), 'LineWidth', 2)
hold on;
plot(tspan, abs(postfit_err(:,1)), 'LineWidth', 2)
title('Error Sun-sensor x-component')
legend('Prefit Residual', 'Postfit Residual')
fontsize(15,'points')
xlabel('Time [sec]')
ylabel('Residual [rad]')

subplot(3,2,3)
plot(tspan, abs(prefit_err(:,2)), 'LineWidth', 2)
hold on;
plot(tspan, abs(postfit_err(:,2)), 'LineWidth', 2)
title('Error Sun-sensor y-component')
legend('Prefit Residual', 'Postfit Residual')
fontsize(15,'points')
xlabel('Time [sec]')
ylabel('Residual [rad]')

subplot(3,2,5)
plot(tspan, abs(prefit_err(:,3)), 'LineWidth', 2)
hold on;
plot(tspan, abs(postfit_err(:,3)), 'LineWidth', 2)
title('Error Sun-sensor z-component')
legend('Prefit Residual', 'Postfit Residual')
fontsize(15,'points')
xlabel('Time [sec]')
ylabel('Residual [rad]')

subplot(3,2,2)
plot(tspan, abs(prefit_err(:,4)), 'LineWidth', 2)
hold on;
plot(tspan, abs(postfit_err(:,4)), 'LineWidth', 2)
title('Error Magnetometer x-component')
legend('Prefit Residual', 'Postfit Residual')
fontsize(15,'points')
xlabel('Time [sec]')
ylabel('Residual [rad]')

subplot(3,2,4)
plot(tspan, abs(prefit_err(:,5)), 'LineWidth', 2)
hold on;
plot(tspan, abs(postfit_err(:,5)), 'LineWidth', 2)
title('Error Magnetometer y-component')
legend('Prefit Residual', 'Postfit Residual')
fontsize(15,'points')
xlabel('Time [sec]')
ylabel('Residual [rad]')

subplot(3,2,6)
plot(tspan, abs(prefit_err(:,6)), 'LineWidth', 2)
hold on;
plot(tspan, abs(postfit_err(:,6)), 'LineWidth', 2)
title('Error Magnetometer z-component')
legend('Prefit Residual', 'Postfit Residual')
fontsize(15,'points')
xlabel('Time [sec]')
ylabel('Residual [rad]')

%% Pre-fit and post fit stats

% Pre-fit error
std_pre1 = std(prefit_err(:,1));
std_pre2 = std(prefit_err(:,2));
std_pre3 = std(prefit_err(:,3));
std_pre4 = std(prefit_err(:,4));
std_pre5 = std(prefit_err(:,5));
std_pre6 = std(prefit_err(:,6));

% post-fit error
std_post1 = std(postfit_err(:,1));
std_post2 = std(postfit_err(:,2));
std_post3 = std(postfit_err(:,3));
std_post4 = std(postfit_err(:,4));
std_post5 = std(postfit_err(:,5));
std_post6 = std(postfit_err(:,6));
