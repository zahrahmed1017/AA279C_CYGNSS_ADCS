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


for i=1:length(t_out)

    % Get fake sensor measurements -- magnetic field vector and sun vector

    [~,~, B_vec] = CalculateMagneticTorque(rv_state(1:3),calday, gmst);




end



