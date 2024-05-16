function [magTorque_eci,B_norm, B] = CalculateMagneticTorque(eci, caldate, gmst)
% Calculating the magnetic torque as a dipole in the inertial frame!!
% T_m = m x B where m is the satellites magnetic moment vector and B is the
% Earth's magnetic field vector.

% Calculating Years since Epoch 1975 
epoch          = datetime(1975,1,1,0,0,0);
timeSinceEpoch = caldate - epoch;
deltaYear      = years(timeSinceEpoch);

% Gaussian Coefficients: [Wertz H-17]
g01 = -30186 + 25.6 * deltaYear; % nT
g11 = -2036 + 10 * deltaYear; % nT
h11 = 5735 - 10.2 * deltaYear; % n

% Dipole strength, coelevation, east longitude [Wertz H-18 through
% H-20]

H0            = sqrt(g01^2 + g11^2 + h11^2); % nT 
theta_prime_m = acos(g01/H0);
phi_prime_m   = atan2(h11,g11);

% Earth vector dipole m_hat [Wertz H-23]
alpha_m = gmst + phi_prime_m; 
m_hat   = [sin(theta_prime_m)*cos(alpha_m); sin(theta_prime_m)*sin(alpha_m); cos(theta_prime_m)];


% Earth Magnetic Field B [Wertz H-22]
if size(eci,2) ~= 1
    eci = eci';
end
R_E    = 6378; %6371.2; % km
R      = norm(eci); % km
R_hat  = eci / R; 
B      = (R_E^3 * H0 * 1e-9/ R^3) * ((3 * dot(m_hat, R_hat) * R_hat) - m_hat); % T 
B_norm = norm(B);


% Satellite Magnetic Moment (hard-coding this for now)
Ssat  = 0.6 * 0.5; % m^2 area of the top face of the bus (ignoring tapering and solar panels)
m_max = (4*pi*1e-7) * 1 * Ssat * 0.1; % Nm/T
% m_max = 0.1;
m     = m_max * [0;0;1];


% Finally, calculate the magnetic torque
magTorque_eci = cross(m,B); % Nm


end
