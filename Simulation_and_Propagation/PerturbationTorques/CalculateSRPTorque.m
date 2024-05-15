function [srpTorque_eci, earth2sun_eci] = CalculateSRPTorque(cygnss, eci , R_i_pa, caldate)

% SRP Force Constants 
P     = 4.57e-6; % [N/m^2] from AA279A Lecture 10 "Extra" Notes 
Ca    = 0.1; % absorption coefficient
Cs    = 0.2; % specular reflection coefficient
Cdiff = 1 - Ca - Cs; % diffuse reflection coefficient
rE    = 6378; % [km] 

% Calculate the unit vector from Earth to Sun in ECI
[earth2sun_unit_eci, earth2sun_eci] = CalculateSunPositionECI(caldate);

% Calculate the unit vector from the spacecraft COM to the Sun in ECI
if size(eci,2) ~= 1
    eci = eci';
end

com2sun_eci        = earth2sun_eci - eci(1:3);
com2sun_unit_eci   = com2sun_eci / norm(com2sun_eci);

% Component Names:
components = {'sp','bus','sensor'};

% Initialize SRP torque:
srpTorque_eci = [0;0;0];

% For each component
for i = 1 : length(components)

    component = components{i};
    faces = fieldnames(cygnss.(component).coord);

    % For each face of each component
    for j = 1 : length(faces)
        face = faces{j};

        facePos_pa     = cygnss.(component).barycenter.(face)'; % m
        faceNormal_pa  = cygnss.(component).normal.(face)'; 
        facePos_eci    = R_i_pa' * facePos_pa;    % m
        faceNormal_eci = R_i_pa' * faceNormal_pa; % m 
        faceArea       = cygnss.(component).area.(face) ; % m^2


        % Check if the surface is illuminated
        e = dot(faceNormal_eci, com2sun_unit_eci);

        % Check if the surface in Earth's shadow
        % Decompose satellite eci position into parallel and
        % perpendicular to sun vector
        r_pll  = dot(eci(1:3), earth2sun_unit_eci) * earth2sun_unit_eci;
        r_perp = eci(1:3) - r_pll; 

        if norm(r_perp) <= rE && dot(r_pll, earth2sun_unit_eci) < 0
            a = 1; % Spacecraft is in the shadow
        else
            a = 0; % Spacecraft not in shadow
        end

        if e > 0 && a == 0 % Surface is illuminated and spacecraft is not in Earth's shadow
            % Calculate SRP force for this face:
            theta = acos(e) / (norm(faceNormal_eci) * norm(com2sun_unit_eci));

            term1     = (1 - Cs) * com2sun_unit_eci; % unitless
            term2     = 2 * faceNormal_eci * (Cs * cos(theta) + (1/3) * Cdiff); % unitless
            Fsrp_surf = -P * cos(theta) * faceArea *...
                         (term1 + term2); % Newtons
    
            % Calculate SRP torque 
            Msrp_surf     = cross(facePos_eci, Fsrp_surf);
            srpTorque_eci = srpTorque_eci + Msrp_surf;
        end

    end


end 