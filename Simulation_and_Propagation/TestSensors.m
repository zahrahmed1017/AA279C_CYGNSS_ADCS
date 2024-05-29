% testing the higher-fidelity sensor models

close all; clear; clc;

load("InertiaData.mat")
load("cygnss.mat")

%% Orbit + attitude prop
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

%% Initialize sun sensor

ss_ang_noise = deg2rad(0.1); % 1 sigma
ss_ang_quant = deg2rad(0.125);
ss = SunSensor(ss_ang_noise, ss_ang_quant, calday);

%% Initialize magnetometer

mag_v_bias = [0, 0, 0]'; % init as zero because we assume good calibration
mag_v_noise = 0.025; % volts
mag_F_matrix = 1e5 * eye(3); % V/T
mag = MagSensor(mag_v_bias, mag_v_noise, mag_F_matrix, calday, gmst);

%% Get sensor measurements for each time step and compare to ground truth

for i=1:length(t_out)
    
    [~,B_norm, B_vec_i] = CalculateMagneticTorque(state_out(i,8:10),calday, gmst);
    B_vec_i = B_vec_i/norm(B_vec_i);
   
    [sun_vec_i, ~] = CalculateSunPositionECI(calday);
    sun_vec_i = sun_vec_i/norm(sun_vec_i);

    R_i_p = quat2dcm(state_out(i,[4 1 2 3]));
    
    B_vec_p = R_i_p * B_vec_i; % ground truth body frame B vector
    sun_vec_p = R_i_p * sun_vec_i; % ground truth body frame sun vector

    B_meas = mag.get_measurement(R_i_p, state_out(i,8:10));

    sun_meas = ss.get_measurement(R_i_p);



end

