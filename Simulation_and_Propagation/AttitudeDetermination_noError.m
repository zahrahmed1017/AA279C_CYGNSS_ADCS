%% Testing attitude determination implementations
close all; clear; clc;

load("InertiaData.mat")
load("cygnss.mat")

%% Deterministic algorithm
% use just magnetometer and sun sensor (assuming full sun sensor access)

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
w_0         = [0, 0, deg2rad(1)]';


%%%% Initial Epoch for Testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];
timeStep     = 0; %seconds
fractionDay  = timeStep/86400;
calday       = datetime(initialEpoch(1),initialEpoch(2), initialEpoch(3),0,0,timeStep);
gmst         = CAL2GMST(initialEpoch(1),initialEpoch(2), initialEpoch(3), fractionDay);

T           = 2*pi*sqrt(a^3/muE); % period in seconds
T_days      = T/(24 * 60 * 60);
numPeriods  = 1;
tspan       = 0 : 5 : T * numPeriods; % simulate once an minute?

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


for i=1:length(t_out)

    % Get fake sensor measurements -- magnetic field vector and sun vector

    [~,B_norm, B_vec_i] = CalculateMagneticTorque(rv_state(1:3),calday, gmst);
    B_vec_i = B_vec_i/norm(B_vec_i);

    [sun_vec_i, ~] = CalculateSunPositionECI(calday);
    sun_vec_i = sun_vec_i/norm(sun_vec_i);

    % rotate into PA frame at this time step
    R_i_p = quat2dcm(state_out(i,1:4));
    B_vec_p = R_i_p * B_vec_i;
    sun_vec_p = R_i_p * sun_vec_i;
    star1_p =  R_i_p * star1_i;
    star2_p =  R_i_p * star2_i;

    % Deterministic method w/ 3 unit vectors
    V1 = [B_vec_i, sun_vec_i, star1_i];
    M1 = [B_vec_p, sun_vec_p, star1_p];
    R_est1 = M1 * inv(V1);
    DCM_est_det1 = cat(3, DCM_est_det1, R_est1);

    % Deterministic method w/ 2 unit vectors and "dummy"
    p_p = B_vec_p;
    q_p = cross(B_vec_p, sun_vec_p) / norm(cross(B_vec_p, sun_vec_p));
    r_p = cross(p_p, q_p);
    p_i = B_vec_i;
    q_i = cross(B_vec_i, sun_vec_i) / norm(cross(B_vec_i, sun_vec_i));
    r_i = cross(p_i, q_i);
    V2 = [p_i, q_i, r_i];
    M2 = [p_p, q_p, r_p];
    R_est2 = M2 * inv(V2);
    DCM_est_det2 = cat(3, DCM_est_det2, R_est2);

    % q-method
    weights = [1 1 2 2];
    weights = weights/norm(weights);
    W = [sqrt(weights); sqrt(weights); sqrt(weights)] .* [B_vec_p, sun_vec_p, star1_p, star2_p];
    U = [sqrt(weights); sqrt(weights); sqrt(weights)] .* [B_vec_i, sun_vec_i, star1_i, star2_i];
    B = W*U';
    S = B + B';
    Z = [B(2,3)-B(3,2), B(3,1)-B(1,3), B(1,2)-B(2,1)]';
    sigma = trace(B);
    K = [S-eye(3)*sigma, Z; 
         Z',       sigma ];
    [v,d] = eig(K);
    [~,ind] = max([d(1,1), d(2,2), d(3,3), d(4,4)]);
    q_est = v(:,ind);
    R_est3 = quat2dcm(q_est([4 1 2 3])');

    % rate measurement
    w_meas = 

end



