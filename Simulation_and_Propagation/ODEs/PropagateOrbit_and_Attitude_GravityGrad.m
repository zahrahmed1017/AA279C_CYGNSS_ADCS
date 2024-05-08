function state_d = PropagateOrbit_and_Attitude_GravityGrad(state, I_p, mu)
% quaternion is [q1 q2 q3 q0] (scalar last)

% Unpack state:
q  = state(1:4);
w  = state(5:7);
rv = state(8:13);

% Calculate Gravity Gradient Torque:
GG_vec = CalculateGravityGradient(q, rv(1:3), mu, I_p);

% Propagate the orbit:
orbit_dot = PropagateOrbit(rv, mu);

% Propagate Kinematics: 
q_d = QuaternionKinematics(q,w);

% Propagate Eulers: 
w_d = PropagateEuler(w, GG_vec, I_p);

state_d = [q_d; w_d; orbit_dot];

end