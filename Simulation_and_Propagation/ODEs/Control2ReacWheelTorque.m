function Lw_d = Control2ReacWheelTorque(Lw, I_p, R_i_pa, w, A, Astar, rv)
    
    earthPointing = false;

    % state     = Lw (angular momentum of the reaction wheel triad)
    % state_dot = reaction wheel torque FOR EACH REACTION WHEEL, NOTTTTT in
    % pa frame

    % Calculate desired input control torque
    if earthPointing

        % solve for n from rv
        muE = 398600;       % [km^3/s^2]
        a_orbit = norm(rv(1:3)) ; % km 
        n = sqrt(muE/(a_orbit^3));

        % get RTN frame
        RTNout      = rv2rtn(rv');
        R_eci_rtn   = [RTNout(1:3)', RTNout(4:6)', RTNout(7:9)' ]';

        R_error     = R_i_pa * R_eci_rtn';
        w_error     = w - [n; 0; 0];
        % w_error     = w;

        control_angs = [R_error(2,3), -R_error(1,3), R_error(1,2)]';

        Mcontrol = controlTorque_earthPointing(I_p, control_angs, w_error, n);


    else
        control_angs = [R_i_pa(2,3), -R_i_pa(1,3), R_i_pa(1,2)]'; % assume a 3-2-1 rotation
        Mcontrol = controlTorque_inertial_discrete(I_p, control_angs, w);
    end
    
    % Calculate actuator torque 
    Lw_d = ComputeActuatorTorque(Lw, Mcontrol, w, A, Astar);

end