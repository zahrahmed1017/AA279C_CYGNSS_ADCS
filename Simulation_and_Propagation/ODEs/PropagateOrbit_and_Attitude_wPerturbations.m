function state_d = PropagateOrbit_and_Attitude_wPerturbations(state, I_p, mu, timeStep, initialEpoch, magTorque, aeroDrag, SRP, gravityGrad)
% quaternion is [q1 q2 q3 q0] (scalar last)

% Unpack state:
q  = state(1:4);
w  = state(5:7);
rv = state(8:13);

% Calculate current time:
fractionDay = timeStep/86400;
caldate     = datetime(initialEpoch(1), initialEpoch(2), initialEpoch(3), 0, 0, timeStep);
gmst        = CAL2GMST(initialEpoch(1),initialEpoch(2),initialEpoch(3), fractionDay);

% Initialize external torques:
M_ext = [0;0;0];

% Calculate Gravity Gradient Torque:
if gravityGrad
    GG_vec = CalculateGravityGradient(q, rv(1:3), mu, I_p);
    M_ext  = M_ext + GG_vec;
end

% Calculate Magnetic Torque
if magTorque
    [magTorque_eci,~] = CalculateMagneticTorque(rv(1:3),caldate,gmst);
end

% Calculate Aerodynamic Drag Torque
if aeroDrag
end

% Calculate SRP Torque
if SRP
end

% Propagate the orbit:
orbit_dot = PropagateOrbit(rv, mu);

% Propagate Kinematics: 
q_d = QuaternionKinematics(q,w);

% Propagate Eulers: 
w_d = PropagateEuler(w, M_ext, I_p);

state_d = [q_d; w_d; orbit_dot];

end