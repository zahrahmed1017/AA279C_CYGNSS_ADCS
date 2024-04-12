%% CYGNSS Orbital Properties:

% Define Orbital Properties
a   = 6903;         % km
e   = 0.00162;
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
numPeriods  = 3;
tspan       = 0 : 60 : T * numPeriods; % simulate once an hour?
initial     = state;
options     = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[tout,Xout] = ode113(@(t,State) PropagateOrbit(State, muE),...
                     tspan, initial, options);

figure()
plot3(Xout(:,1),Xout(:,2),Xout(:,3),'Color','r','LineWidth',3)
hold on;
grid on;
axis equal;
xlabel('X')
ylabel('Y')
zlabel('Z')
% title(['CYGNSS Orbit Propagation with No Perturbing Forces for ', num2str(numPeriods), ' orbits']);
fontsize(14,'points')

%plot Earth-sized sphere
rE = 6378;
[xE,yE,zE] = ellipsoid(0,0,0,rE,rE,rE,20);
surface(xE, yE , zE , 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha', 0.25);
view (3);
