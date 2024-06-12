function Lw_d = Control2ReacWheelTorque(Lw, I_p, R_i_pa, w, A, Astar)
    
    % state     = Lw (angular momentum of the reaction wheel triad)
    % state_dot = reaction wheel torque FOR EACH REACTION WHEEL, NOTTTTT in
    % pa frame

    % Calculate desired input control torque
    control_angs = [R_i_pa(2,3), -R_i_pa(1,3), R_i_pa(1,2)]'; % assume a 3-2-1 rotation
    Mcontrol = controlTorque_inertial(I_p, control_angs, w);
    
    % Calculate actuator torque 
    Lw_d = ComputeActuatorTorque(Lw, Mcontrol, w, A, Astar);

end