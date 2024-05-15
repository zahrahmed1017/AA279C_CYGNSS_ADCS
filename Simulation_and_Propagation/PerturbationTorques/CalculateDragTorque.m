function aeroDragTorque_eci = CalculateDragTorque(cygnss, eci , w_sat_pa , R_i_pa)

% Constants
r_E   = 6378137; % [m]
H     = 7500; % [m], characteristic height
rho_0 = 1.225 ; %kg/m^3
Cd    = 2; 
w_E   = [0; 0; 7.292e-5]; % rad/s rotation of earth in ECI

% ECI Position and velocity of the satellite 
satPos_eci = eci(1:3); % km
satVel_eci = eci(4:6); % km/s
r_mag      = norm(satPos_eci); % km

% Convert angular velocity of satellite from PA to inertial:
w_sat_eci = R_i_pa' * w_sat_pa; % rad/s

% Velocity of COM wrt Earth
satVel_wrtE = satVel_eci - cross(w_E, satPos_eci); % km/s

% Atmospheric density function:
% rho = rho_0 * exp(-((r_mag*1000) - r_E)/H); % kg / m^3
rho = 1e-11;

% Component Names:
components = {'sp','bus','sensor'};

% Initialize Drag torque:
aeroDragTorque_eci = [0;0;0];

% For each component
for i = 1 : length(components)

    component = components{i};
    faces = fieldnames(cygnss.(component).coord);

    % For each face of each component
    for j = 1 : length(faces)
        face = faces{j};

        facePos_pa    = cygnss.(component).barycenter.(face)'; % m
        faceNormal_pa = cygnss.(component).normal.(face)'; 
        faceArea      = cygnss.(component).area.(face) ; % m^2

        facePos_eci    = R_i_pa' * facePos_pa;    % m
        faceNormal_eci = R_i_pa' * faceNormal_pa; % m 
        
        % Velocity of surface wrt Earth (assuming no winds)
        surfVel_wrtE = satVel_wrtE + cross(w_sat_eci, facePos_eci/1000); % km/s
        surfVel_norm = norm(surfVel_wrtE);                               % km/s 
        surfVel_unit = surfVel_wrtE / surfVel_norm; 

        % Calculate Drag force for this face:
        Fdrag_surf   = -0.5 * rho * (surfVel_norm * 1000)^2 * faceArea * Cd * dot(surfVel_unit, faceNormal_eci) * surfVel_unit; % Newtons

        % Calculate Drag torque 
        Maero_surf   = cross(facePos_eci, Fdrag_surf);

        aeroDragTorque_eci = aeroDragTorque_eci + Maero_surf;

    end
    
end