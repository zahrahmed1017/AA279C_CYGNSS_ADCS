function [state_d, gravGrad_pa, magTorque_pa, dragTorque_pa, srpTorque_pa] =...
    PropagateOrbit_Attitude_wPert_wControl(state, I_p, mu, cygnss, A, Astar, timeStep,...
    initialEpoch, gravityGrad, magTorque, aeroDrag, SRP, control)
% quaternion is [q1 q2 q3 q0] (scalar last)

% Unpack state:
q  = state(1:4);
w  = state(5:7);
rv = state(8:13);
Lw = state(14:16);

% calculate mean motion n (hacky)
muE = 398600;       % [km^3/s^2]
a = sqrt(rv(1)^2 + rv(2)^2 + rv(3)^2);
n = sqrt(muE/(a^3));

% Calculate rotation matrix from inertial to principal
R_i_pa = quat2dcm(q([4 1 2 3])');

% Calculate current time:
fractionDay = timeStep/86400;
caldate     = datetime(initialEpoch(1), initialEpoch(2), initialEpoch(3), 0, 0, timeStep);
gmst        = CAL2GMST(initialEpoch(1), initialEpoch(2), initialEpoch(3), fractionDay);

% Calculate control torque

% small angle approx. of angle error -- assume target attitude is [0 0 0 1]
% this means the current attitude IS the error!
% for now, assume perfect attitude knowledge
control_angs = [R_i_pa(2,3), -R_i_pa(1,3), R_i_pa(1,2)]'; % assume a 3-2-1 rotation
Mcontrol = controlTorque_inertial(I_p, control_angs, w);

% Calculate Actuator Torque
actuatorTorque_pa = [0;0;0];
Lw_d = [0;0;0];
M_mag = [0;0;0];
if control
    Lw_d = ComputeActuatorTorque(Lw, Mcontrol, w, A, Astar);
    actuatorTorque_pa = A * Lw_d;

    kmag = 7e-13;
    L_i = I_p * w; % TODO does this need to have cross product term in it?
    L_tar = I_p * [0; -n; 0];
    dH = L_i - L_tar;
    [~,B_norm, B_vec_i] = CalculateMagneticTorque(rv(1:3) ,caldate, gmst);
    B_p = R_i_pa * B_vec_i;
    M_mag = ComputeMagnetorquerTorque(dH, B_p, kmag);
end

% Calculate Gravity Gradient Torque (already in principal coordinates)
gravGrad_pa = [0;0;0];
if gravityGrad
    gravGrad_pa = CalculateGravityGradient(q, rv(1:3), mu, I_p);
end

% Calculate Magnetic Torque
magTorque_pa = [0;0;0];
if magTorque
    [magTorque_eci,~] = CalculateMagneticTorque(rv(1:3),caldate,gmst);
    magTorque_pa      = R_i_pa * magTorque_eci;
end

% Calculate Aerodynamic Drag Torque
dragTorque_pa = [0;0;0];
if aeroDrag
    dragTorque_eci = CalculateDragTorque(cygnss, rv, w, R_i_pa);
    dragTorque_pa  = R_i_pa * dragTorque_eci;
end

% Calculate SRP Torque
srpTorque_pa = [0;0;0];
if SRP
    srpTorque_eci = CalculateSRPTorque(cygnss, rv , R_i_pa, caldate);
    srpTorque_pa  = R_i_pa * srpTorque_eci; 
end

% Total External Torque:
M_ext = actuatorTorque_pa + M_mag + gravGrad_pa + magTorque_pa + dragTorque_pa + srpTorque_pa; 
% M_ext = -Mcontrol + gravGrad_pa + magTorque_pa + dragTorque_pa + srpTorque_pa;

% Propagate the orbit:
orbit_dot = PropagateOrbit(rv, mu);

% Propagate Kinematics: 
q_d = QuaternionKinematics(q,w);

% Propagate Eulers: 
w_d = PropagateEuler(w, M_ext, I_p);

state_d = [q_d; w_d; orbit_dot; Lw_d];

end