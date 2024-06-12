
close all; clear;

load("InertiaData.mat")
load("cygnss.mat")
load("measurementCov.mat")

%% Initial state - orbit, attitude, rates

%%%%%%%% INITIAL ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
a   = rE + 525;
e   = 0;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

oe          = [a; e; i; O; w; v];
rv_state    = OE2ECI(oe, muE);
RTNout      = rv2rtn(rv_state');
A_eci_rtn   = [RTNout(1:3)', RTNout(4:6)', RTNout(7:9)' ]';

% TODO tweak the initial attitude + rates slightly so there is some error for
% controller to correct

%%%%%%% INITIAL ATTITUDE AND RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth pointing attitude
R_eci_nadirRTN = [0 0 1; 0 1 0; -1 0 0] * A_eci_rtn;
q_0 = dcm2quat(R_eci_nadirRTN);
q_0 = q_0([2 3 4 1])';
w_0 = [0, -sqrt(muE/(a^3)),  0]';

% Inertial pointing attitude
% q_0 = [0, 0, 0, 1]';
% w_0 = deg2rad([0.5, 0.1, 1]');

%%%% SET SIMULATION TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
numPeriods  = 0.01;
dt          = 0.05; %seconds
tspan       = 0 : dt : T * numPeriods; % simulate once an minute?
% tspan      = 0:0.5:200;

%%%% INITIALIZE EPOCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];
% timeStep     = .5;
fractionDay  = dt/86400;
calday       = datetime(initialEpoch(1),initialEpoch(2), initialEpoch(3),0,0,dt);
gmst         = CAL2GMST(initialEpoch(1),initialEpoch(2), initialEpoch(3), fractionDay);


%% Attitude Determination: MEKF Initialization

% Initial State and Covariance
state_0        = [q_0; w_0];

cov_0          = [1e-3, 0,    0,    0,    0,    0;...
                   0,    1e-3, 0,    0,    0,    0;...
                   0,    0,    1e-3, 0,    0,    0;...
                   0,    0,    0,    1e-12,0,   0;...
                   0,    0,    0,    0,    1e-12,0;...
                   0,    0,    0,    0,    0,     1e-12];
cov_0 = cov_0 * 50;

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

% Control Matrix B:
B = [1/I_p(1,1), 0, 0; 0 1/I_p(2,2) 0; 0, 0, 1/I_p(3,3)];

%% Attitude Control: Spacecraft actuators

% reaction wheel triad
A = [1/sqrt(2), 1/sqrt(3),  1/sqrt(3);...
     0        , 1/sqrt(3), -1/sqrt(3);...
    -1/sqrt(2), 1/sqrt(3),  1/sqrt(3)];
Astar = pinv(A);

% Magnetometer
% TODO

Lw_0  = [0; 0; 0];

% TODO: anything needed to initialize sensors? error settings?

%% Simulation settings

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
gravityGrad     = 1;
magTorque       = 1;
aeroDrag        = 1;
SRP             = 1;
control         = 1;

%% Run the simulation
state_0 = [q_0; w_0; rv_state; Lw_0];
M_vec   = [0;0;0];

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

    timeStep = tspan(i);
    tspan_oneStep = [tspan(i) tspan(i+1)];

    %%%%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    reacWheelTorque_pa  = [0;0;0];
    magtorquerTorque_pa = [0;0;0];

    if ReacWheelOn
         % Calculate the actual reaction wheel torque that will be used for
         % control
         reacWheelTorque    = Control2ReacWheelTorque(Lw_0, I_p, R_i_pa, w, A, Astar);
         reacWheelTorque_pa = A * reacWheelTorque;

         % Now integrate to find the new angular momentum of the reaction
         % wheels if we were to apply that control torque given our current
         % state and reset Lw_0 so at the next time step we have our
         % updated reaction wheel angular momentum. 
         [t_prop, Lw_out] = ode113(@(t,Lw) Control2ReacWheelTorque(Lw,...
            I_p, R_i_pa, w, A, Astar),...
            tspan_oneStep, Lw_0, options);

         Lw_0 = Lw_out(end,:);
    end

    if MagnetorquerOn
        magtorquerTorque_pa = MagnetorquerWrapper(R_i_pa, rv, calday, gmst);
    end

    appliedControl = reacWheelTorque_pa + magtorquerTorque_pa; 

    %%%%%%%%%%%%% PROPAGATE GROUND TRUTH  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [gtTimeOut, gtState] = ode113(@(t,gtState)IntegratedODE(state, I_p, mu,...
        cygnss, timeStep, initialEpoch, appliedControl,...
        gravityGrad, magTorque, aeroDrag, SRP), tspan_oneStep, gtState_0, options);
    

    gtState_0 = gtState(end,:);


    %%%%%%%%%%%% NAVIGATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%% STEP 2 filter.Propagate() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

   %%%%%%%%%%%%% STEP 3 filter.TimeUpdate() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    obj.errorState(4:6) = obj.A22 * obj.errorState(4:6) + (B * Mrot * obj.dt);
    w_mekf(i,:)         = obj.errorState(4:6)';
            
    %%%% Update State Covariance
    obj.stateCov = obj.STM * obj.stateCov * obj.STM' + obj.processNoise;

    cov_mekf(i,1) = sqrt(obj.stateCov(1,1));
    cov_mekf(i,2) = sqrt(obj.stateCov(2,2));
    cov_mekf(i,3) = sqrt(obj.stateCov(3,3));

   % %%%%%%%%%%%%% STEP 4: filter.MeasurementUpdate() %%%%%%%%%%%%%%%%%%%%%%%

   %%% Get true measurements - 3n x 1 vector where n is the number of
   %%% measurements
   R_i_p    = quat2dcm(q_t1([4 1 2 3])'); % Ground truth attitude

   %%% Get modeled measurements - 3n x 1 vector where n is the number of
   %%% measurements
   quat_est = obj.state(1:4)';
   R_est    = quat2dcm(quat_est([4 1 2 3]));

   [sun_vec_i, ~]       = CalculateSunPositionECI(calday);
   sun_vec_i            = sun_vec_i/norm(sun_vec_i);
   [~, B_norm, B_vec_i] = CalculateMagneticTorque(w_t1,calday, gmst);
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

   %%%%%%%%%%%% STEP 5: filter.Reset() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
   [~,B_norm, B_vec_i] = CalculateMagneticTorque(w_t1,calday, gmst);
   B_vec_i             = B_vec_i/norm(B_vec_i);

   z1_post = R_est * sun_vec_i; % Sun sensor modeled measurement
   z2_post = R_est * B_vec_i;   % Magnetometer modeled measurement
   z_post  = [z1_post; z2_post];

   postfit_err(i,:) = y' - z_post';

   %%%%%%%%%%%% Update ODE State for next iteration %%%%%%%%%%%%%%%%%%%%%%%
   state_0 = [q_t1; w_t1; rv_t1; Lw_t1];
end

%% Plot results afterwards (anything we can't recover directly from state output)

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
ylim([-1, 1])
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
ylim([-1, 1])
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
