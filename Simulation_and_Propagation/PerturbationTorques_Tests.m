%% Gravity Gradient Torque

close all; clear;

load("InertiaData.mat")
load("cygnss.mat")

%% Testing Magnetic Torque Values for different altitudes 

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

%%%% Initial Epoch for Testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];
timeStep     = 0; %seconds
fractionDay  = timeStep/86400;
calday       = datetime(initialEpoch(1),initialEpoch(2), initialEpoch(3),0,0,timeStep);
gmst         = CAL2GMST(initialEpoch(1),initialEpoch(2), initialEpoch(3), fractionDay);

[magTorque_eci,earthMagField] = CalculateMagneticTorque(rv_state(1:3),calday, gmst);


% Testing for different altitudes
altitudes         = 500:100:30000;
magTorque_eci_all = zeros(length(altitudes),4);
magField_magnitude = zeros(length(altitudes),1);

for i = 1:length(altitudes)

    h        = altitudes(i);
    a        = rE + h;
    % oe       = [a; e; i; O; w; v];
    % rv_state = OE2ECI(oe, muE);

    rv_state = [a; 0; 0];

    [magTorque_eci,earthMagField] = CalculateMagneticTorque(rv_state(1:3), calday, gmst);

    magTorque_eci_all(i,1:3) = magTorque_eci;
    magTorque_eci_all(i,4)   = norm(magTorque_eci);
    magField_magnitude(i)    = earthMagField;
end

figure()
subplot(1,2,1)
plot(altitudes,magField_magnitude, 'LineWidth', 2)
title('Earth Magnetic Field Magnitude versus Altitude')
grid on;
fontsize(14,'points')
xlabel('Altitude [km]')
ylabel('Magnetic Field Magnitude [T]')

subplot(1,2,2)
plot(altitudes,magTorque_eci_all(:,4), 'LineWidth', 2)
title('Magnetic Torque Magnitude versus Altitude')
grid on;
fontsize(14,'points')
xlabel('Altitude [km]')
ylabel('Magnetic Torque Magnitude [Nm]')


%% Calculate Magnetic Torque Magnitude over a single orbit:

% Define Orbital Properties
rE  = 6378; % km
h   = 525;
a   = rE + h;
e   = 0;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

T      = 2*pi*sqrt(a^3/muE); % period in seconds
T_days = T/(24 * 60 * 60);

oe = [a; e; i; O; w; v];

% Convert to ECI position and velocity for orbit propagation 
state = OE2ECI(oe, muE);

% Propagate Orbit for multiple orbits
numPeriods  = 1;
tspan       = 0 : 60 : T * numPeriods; % simulate once an minute?
initial     = state;
options     = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[tout,Xout] = ode113(@(t,State) PropagateOrbit(State, muE),...
                     tspan, initial, options);


initialEpoch = [2016, 12, 16];
magTorque_1orbit_eci = zeros(length(tout),4);

for j = 1:length(tout)

    fractionDay = tout(j)/86400;
    calday      = datetime(initialEpoch(1),initialEpoch(2), initialEpoch(3),0,0,timeStep);
    gmst        = CAL2GMST(initialEpoch(1),initialEpoch(2), initialEpoch(3), fractionDay);

    [magTorque_eci,earthMagField] = CalculateMagneticTorque(Xout(j,1:3), calday, gmst);

    magTorque_1orbit_eci(j,1:3) = magTorque_eci;
    magTorque_1orbit_eci(j,4)   = norm(magTorque_eci);

end

totalPoints = length(tout);
step        = totalPoints/4;
markers     = 0 : step : totalPoints;



figure()
subplot(1,2,1)
plot(tout/3600,magTorque_1orbit_eci(:,4),'LineWidth',2)
title('Magnetic Torque Magnitude versus Time')
grid on;
hold on;
plot(tout(markers(2))/3600, magTorque_1orbit_eci(markers(2),4), 'o', 'MarkerSize',15,'MarkerFaceColor', 'r')
plot(tout(markers(3))/3600, magTorque_1orbit_eci(markers(3),4), 'o', 'MarkerSize',15,'MarkerFaceColor', 'b')
plot(tout(markers(4))/3600, magTorque_1orbit_eci(markers(4),4), 'o', 'MarkerSize',15,'MarkerFaceColor', 'g')
plot(tout(markers(5))/3600, magTorque_1orbit_eci(markers(5),4), 'o', 'MarkerSize',15,'MarkerFaceColor', 'c')
fontsize(14,'points')
xlabel('Time [hours]')
ylabel('Magnetic Torque Magnitude [Nm]')

subplot(1,2,2)
plot3(Xout(:,1),Xout(:,2),Xout(:,3),'Color','r','LineWidth',3)
title('ECI Position')
hold on;
grid on;
axis equal;
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
% title(['CYGNSS Orbit Propagation with No Perturbing Forces for ', num2str(numPeriods), ' orbits']);
fontsize(14,'points')

%plot Earth-sized sphere
rE = 6378;
[xE,yE,zE] = ellipsoid(0,0,0,rE,rE,rE,20);
surface(xE, yE , zE , 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha', 0.25);
view (3);

% Time markers
plot3(Xout(markers(2),1), Xout(markers(2),2), Xout(markers(2),3), 'o', 'MarkerSize',15,'MarkerFaceColor', 'r')
plot3(Xout(markers(3),1), Xout(markers(3),2), Xout(markers(3),3), 'o', 'MarkerSize',15,'MarkerFaceColor', 'b')
plot3(Xout(markers(4),1), Xout(markers(4),2), Xout(markers(4),3), 'o', 'MarkerSize',15,'MarkerFaceColor', 'g')
plot3(Xout(markers(5),1), Xout(markers(5),2), Xout(markers(5),3), 'o', 'MarkerSize',15,'MarkerFaceColor', 'c')


%% Testing Aerodynamic Torque 

%%%% INITIALIZE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
% a   = rE + 3000;     % km
a   = rE + 525;
e   = 0;
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
numPeriods  = 1;
tspan       = 0 : 60 : T * numPeriods; % simulate once an minute?

%%%% INITIALIZE ATTITUDE AND ANGULAR RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial angular velocity:
w_0         = [0, 0, deg2rad(5)]';

% Initial attitude:
e_vec = [1;1;1]; 
e     = e_vec/ norm(e_vec);
p     = 0; % for PS4-Q1
q_0   = [e(1)*sin(p/2);
         e(2)*sin(p/2);
         e(3)*sin(p/2);
         cos(p/2)];
dcm_0 = quat2dcm(q_0([4 1 2 3])');

aeroDragTorque_eci = CalculateDragTorque(cygnss, rv_state , w_0 , dcm_0);

%% INITIAL EPOCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];


%%% WITHOUT DRAG TORQUE:
state_0   = [q_0; w_0; rv_state];
gravityGrad = 0;
magTorque   = 0;
aeroDrag    = 0;
SRP         = 0;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t_noT, state_out_noT] = ...
    ode113(@(t,state) PropagateOrbit_and_Attitude_wPerturbations(state, I_p, muE, cygnss, t, initialEpoch, gravityGrad, magTorque, aeroDrag, SRP),...
    tspan, state_0, options);

%%% WITH DRAG TORQUE:

state_0   = [q_0; w_0; rv_state];
gravityGrad = 0;
magTorque   = 0;
aeroDrag    = 1;
SRP         = 0;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t, state_out] = ...
    ode113(@(t,state) PropagateOrbit_and_Attitude_wPerturbations(state, I_p, muE, cygnss, t, initialEpoch, gravityGrad, magTorque, aeroDrag, SRP),...
    tspan, state_0, options);

eulerAngs_noT = zeros(length(t_noT), 3);
eulerAngs     = zeros(length(t), 3);
aeroDrag      = zeros(length(t), 4);

for i = 1:length(t_noT)
    
    q_noT = state_out_noT(i,1:4);
    q_wT  = state_out(i,1:4);

    [yaw_noT, pitch_noT, roll_noT] = quat2angle(q_noT([4 1 2 3]), 'ZYX');
    [yaw, pitch, roll]             = quat2angle(q_wT([4 1 2 3]), 'ZYX');

    eulerAngs_noT(i, :) = [yaw_noT, pitch_noT, roll_noT];
    eulerAngs(i,:)      = [yaw, pitch, roll];

    % Calculate drag 
    dcm_i = quat2dcm(q_wT([4 1 2 3]));
    w_i   = state_out(i,5:7)';
    rv_i  = state_out(i,8:13)';

    drag_eci = CalculateDragTorque(cygnss, rv_i , w_i , dcm_i);
    drag_pa  = dcm_i * drag_eci;
    drag_norm = norm(drag_pa);

    aeroDrag(i,1:3) = drag_pa';
    aeroDrag(i,4)   = drag_norm;
    
end

figure();
subplot(1,3,1)
plot(t/3600, eulerAngs_noT(:,1), 'LineWidth', 2)
hold on;
plot(t/3600, eulerAngs(:,1), 'LineWidth', 2)
title('Yaw')
xlabel('Time [hr]')
ylabel('Angle [rad]')
legend('Without Drag Torque', 'With Drag Torque')
fontsize(14,'points')

subplot(1,3,2)
plot(t/3600, eulerAngs_noT(:,2), 'LineWidth', 2)
hold on;
plot(t/3600, eulerAngs(:,2), 'LineWidth', 2)
title('Pitch')
xlabel('Time [hr]')
ylabel('Angle [rad]')
legend('Without Drag Torque', 'With Drag Torque')
fontsize(14,'points')

subplot(1,3,3)
plot(t/3600, eulerAngs_noT(:,3), 'LineWidth', 2)
hold on;
plot(t/3600, eulerAngs(:,3), 'LineWidth', 2)
title('Roll')
xlabel('Time [hr]')
ylabel('Angle [rad]')
legend('Without Drag Torque', 'With Drag Torque')
fontsize(14,'points')

figure();
plot(t/3600, aeroDrag(:,1), 'LineWidth',2)
title('Aerodynamic Drag Components and Magnitude')
grid on;
hold on;
plot(t/3600, aeroDrag(:,2), 'LineWidth',2)
plot(t/3600, aeroDrag(:,3), 'LineWidth',2)
plot(t/3600, aeroDrag(:,4), 'LineWidth',2)
xlabel('Time [hr]')
ylabel('Magnitude [Nm]')
fontsize(14,'points')

legend('X Component', 'Y Component', 'Z Component', 'Magnitude')


%% Now Simulate Orbit and Attitude with Gravity Gradient, Magnetic Torque, and Aero Drag Together

%%%% INITIALIZE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rE  = 6378; % km
% a   = rE + 3000;     % km
a   = rE + 525;
e   = 0;
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
numPeriods  = 1;
tspan       = 0 : 60 : T * numPeriods; % simulate once an minute?

%%%% INITIALIZE ATTITUDE AND ANGULAR RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial angular velocity:
w_0         = [0, 0, deg2rad(5)]';

% Initial attitude:
e_vec = [1;1;1]; 
e     = e_vec/ norm(e_vec);
p     = 0; % for PS4-Q1
q_0   = [e(1)*sin(p/2);
         e(2)*sin(p/2);
         e(3)*sin(p/2);
         cos(p/2)];
dcm_0 = quat2dcm(q_0([4 1 2 3])');

%%%% INITIALIZE EPOCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialEpoch = [2016, 12, 16];

%%% RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_0   = [q_0; w_0; rv_state];
gravityGrad = 1;
magTorque   = 1;
aeroDrag    = 1;
SRP         = 0;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t, state_out] = ...
    ode113(@(t,state) PropagateOrbit_and_Attitude_wPerturbations(state, I_p, muE, cygnss, t, initialEpoch, gravityGrad, magTorque, aeroDrag, SRP),...
    tspan, state_0, options);

eulerAngs_noT = zeros(length(t), 3);
eulerAngs     = zeros(length(t), 3);
gravGrad      = zeros(length(t), 4);
magTorque     = zeros(length(t), 4);
aeroDrag      = zeros(length(t), 4);
totTorque     = zeros(length(t), 4);


for i = 1 : length(t)

    q_i   = state_out(i,1:4);
    w_i   = state_out(i,5:7);
    rv_i  = state_out(i,8:13);
    dcm_i = quat2dcm(q_i([4 1 2 3]));

    % Calculate Gravity Gradient:
    gg     = CalculateGravityGradient(q_i([4 1 2 3])', rv_i(1:3)' , muE , I_p); % 1 x 3
    gg_mag = norm(gg);

    gravGrad(i,:) = [gg, gg_mag];


    % Calculate Magnetic Torque:
    fractionDay = t(i)/86400;
    calday      = datetime(initialEpoch(1),initialEpoch(2), initialEpoch(3), 0, 0, t(i));
    gmst        = CAL2GMST(initialEpoch(1),initialEpoch(2), initialEpoch(3), fractionDay);

    [magTorque_eci,earthMagField] = CalculateMagneticTorque(rv_i(1:3), calday, gmst);

    magTorque_pa = dcm_i * magTorque_eci; % 3 x 1
    magTorque(i,1:3) = magTorque_pa;
    magTorque(i,4)   = norm(magTorque_pa);


    % Calculate Aerodynamic Drag:
    drag_eci  = CalculateDragTorque(cygnss, rv_i' , w_i' , dcm_i);
    drag_pa   = dcm_i * drag_eci;
    drag_norm = norm(drag_pa);

    aeroDrag(i,1:3) = drag_pa';
    aeroDrag(i,4)   = drag_norm;

    % Total Torques:
    tot_torque_i = gg + magTorque_pa' + drag_pa';
    tot_torque_norm = norm(tot_torque_i);

    totTorque(i,:) = [tot_torque_i, tot_torque_norm];

end

% % Calculate Max Theoretical Values:
% maxGravGrad = (3/2) * muE / (norm(rv_state(1:3))^3) * abs(I_p(3,3) - I_p(1,1));
% maxMag      = 2 * ((4*pi*1e-7) * 1 * 0.5 * 0.6 * 0.1) * 7.9e15 / (norm(rv_state(1:3)/1000)^3) ; 
% maxDrag     = (0.5) * 1e-11 * (0.5*0.6) * 2 ;


% Plot Components for Each Modeled Torque Over Time:
t_hr = t/3600;

figure()
% Gravity Gradient
subplot(3,1,1)
plot(t_hr, gravGrad(:,1), 'LineWidth',2)
hold on;
grid on;
plot(t_hr, gravGrad(:,2), 'LineWidth',2)
plot(t_hr, gravGrad(:,3), 'LineWidth', 2)
plot(t_hr, gravGrad(:,4), 'LineWidth', 2)
title('Gravity Gradient Components in Principal Frame')
xlabel('Time [hr]')
ylabel('Magnitude [Nm]')
legend('X Component', 'Y Component', 'Z Component', 'Magnitude')
fontsize(14,'points')

% Magnetic Torque
subplot(3,1,2)
plot(t_hr, magTorque(:,1), 'LineWidth',2)
hold on;
grid on;
plot(t_hr, magTorque(:,2), 'LineWidth',2)
plot(t_hr, magTorque(:,3), 'LineWidth', 2)
plot(t_hr, magTorque(:,4), 'LineWidth', 2)
title('Magnetic Torque Components in Principal Frame')
xlabel('Time [hr]')
ylabel('Magnitude [Nm]')
legend('X Component', 'Y Component', 'Z Component', 'Magnitude')
fontsize(14,'points')

% Aerodynamic Drag Torque
subplot(3,1,3)
plot(t_hr, aeroDrag(:,1), 'LineWidth',2)
hold on;
grid on;
plot(t_hr, aeroDrag(:,2), 'LineWidth',2)
plot(t_hr, aeroDrag(:,3), 'LineWidth', 2)
plot(t_hr, aeroDrag(:,4), 'LineWidth', 2)
title('Drag Torque Components in Principal Frame')
xlabel('Time [hr]')
ylabel('Magnitude [Nm]')
legend('X Component', 'Y Component', 'Z Component', 'Magnitude')
fontsize(14,'points')

% Total Torque
figure()
plot(t_hr, totTorque(:,1), 'LineWidth',2)
hold on;
grid on;
plot(t_hr, totTorque(:,2), 'LineWidth',2)
plot(t_hr, totTorque(:,3), 'LineWidth', 2)
plot(t_hr, totTorque(:,4), 'LineWidth', 2)
title('Total Torque Components in Principal Frame')
xlabel('Time [hr]')
ylabel('Magnitude [Nm]')
legend('X Component', 'Y Component', 'Z Component', 'Magnitude')
fontsize(14,'points')