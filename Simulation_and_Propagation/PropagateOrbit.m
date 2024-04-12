function newState = PropagateOrbit(state, mu)

% No perturbing forces included, just the fundamental orbital differential
% equation

% Parse state vector:
x  = state(1);
y  = state(2);
z  = state(3);
vx = state(4);
vy = state(5);
vz = state(6);

% Calculate Acceleration:
r_mag = norm([x,y,z]);
ax = -mu * x/r_mag^3;
ay = -mu * y/r_mag^3;
az = -mu * z/r_mag^3;

% Package derivatives
newState = zeros(6,1);
newState(1) = vx;
newState(2) = vy;
newState(3) = vz;
newState(4) = ax;
newState(5) = ay;
newState(6) = az;


end