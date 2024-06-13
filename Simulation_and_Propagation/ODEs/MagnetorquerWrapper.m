function M_mag = MagnetorquerWrapper(R_i_pa, I_p, rv, w, calday, gmst)

    % Given the current ECI position and time, computes the
    % magnetorquer torque

    % calculate mean motion n (hacky)
    muE = 398600;       % [km^3/s^2]
    a = sqrt(rv(1)^2 + rv(2)^2 + rv(3)^2);
    n = sqrt(muE/(a^3));

    kmag = 7e-13;
    L_i = I_p * w; % TODO does this need to have cross product term in it?
    L_tar = I_p * [0; -n; 0];
    dH = L_i - L_tar;
    [~,~, B_vec_i] = CalculateMagneticTorque(rv(1:3) ,calday, gmst);
    B_p = R_i_pa * B_vec_i;
    M_mag = ComputeMagnetorquerTorque(dH, B_p, kmag);
    
end 