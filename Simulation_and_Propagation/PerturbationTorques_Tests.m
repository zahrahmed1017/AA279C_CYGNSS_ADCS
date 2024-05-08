%% Gravity Gradient Torque

close all; clear;

load("InertiaData.mat")

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
magTorque_1orbit = zeros(length(tout),4);

for j = 1:length(tout)

    fractionDay = tout(j)/86400;
    calday      = datetime(initialEpoch(1),initialEpoch(2), initialEpoch(3),0,0,timeStep);
    gmst        = CAL2GMST(initialEpoch(1),initialEpoch(2), initialEpoch(3), fractionDay);

    [magTorque_eci,earthMagField] = CalculateMagneticTorque(Xout(j,1:3), calday, gmst);

    magTorque_1orbit(j,1:3) = magTorque_eci;
    magTorque_1orbit(j,4)   = norm(magTorque_eci);

end

totalPoints = length(tout);
step        = totalPoints/4;
markers     = 0 : step : totalPoints;



figure()
subplot(1,2,1)
plot(tout/3600,magTorque_1orbit(:,4),'LineWidth',2)
title('Magnetic Torque Magnitude versus Time')
grid on;
hold on;
plot(tout(markers(2))/3600, magTorque_1orbit(markers(2),4), 'o', 'MarkerSize',15,'MarkerFaceColor', 'r')
plot(tout(markers(3))/3600, magTorque_1orbit(markers(3),4), 'o', 'MarkerSize',15,'MarkerFaceColor', 'b')
plot(tout(markers(4))/3600, magTorque_1orbit(markers(4),4), 'o', 'MarkerSize',15,'MarkerFaceColor', 'g')
plot(tout(markers(5))/3600, magTorque_1orbit(markers(5),4), 'o', 'MarkerSize',15,'MarkerFaceColor', 'c')
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


