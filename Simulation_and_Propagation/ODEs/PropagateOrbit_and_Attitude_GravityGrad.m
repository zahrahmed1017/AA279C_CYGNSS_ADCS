function state_d = PropagateOrbit_and_Attitude_GravityGrad(state, I_p, mu)
% quaternion is [q1 q2 q3 q0] (scalar last)

% Unpack state:
q  = state(1:4);
w  = state(5:7);
rv = state(8:13);

% % Calculate position vector in principal coordinates:
% R_i_pa = quat2dcm(q([4 1 2 3])');
% r_pa   = R_i_pa * rv(1:3);
% 
% % Calculate gravity gradient torque:
% r_norm = norm(r_pa);
% GG_x = (3*mu/r_norm) * (I_p(3,3) - I_p(2,2)) * r_pa(2) * r_pa(3);
% GG_y = (3*mu/r_norm) * (I_p(1,1) - I_p(3,3)) * r_pa(3) * r_pa(1);
% GG_z = (3*mu/r_norm) * (I_p(2,2) - I_p(1,1)) * r_pa(1) * r_pa(2);
% GG_vec = [GG_x, GG_y, GG_z];

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