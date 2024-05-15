function gravityGrad = CalculateGravityGradient(q, r_eci, mu, I_p)

% Inputs:
    % q   - quaternion ([q1; q2; q3; q0], scalar first) from inertial to
    %      principal
    % r   - inertial (ECI) position as a column vector (km)
    % mu  - central body gravitational constant (km^3/s^2)
    % I_p - Principal Moments of Inertia (kg * m^2)

    % OUTPUT GG in principal coordinates IN N-m !!!!

    % Convert mu to m^3/s^2
    mu_m    = mu * 1e9; % m^3/s^2
    r_eci_m = r_eci * 1000; 


    % Calculate position vector in principal coordinates:
    R_i_pa = quat2dcm(q([4 1 2 3])');
    r_pa   = R_i_pa * r_eci_m;
    c_pa   = r_pa / norm(r_pa); % Unit Vector

% Calculate gravity gradient torque:
    r_norm = norm(r_pa);
    GG_x = (3 * mu_m /r_norm^3) * (I_p(3,3) - I_p(2,2)) * c_pa(2) * c_pa(3);
    GG_y = (3 * mu_m /r_norm^3) * (I_p(1,1) - I_p(3,3)) * c_pa(3) * c_pa(1);
    GG_z = (3 * mu_m /r_norm^3) * (I_p(2,2) - I_p(1,1)) * c_pa(1) * c_pa(2);
    gravityGrad = [GG_x, GG_y, GG_z];

end
