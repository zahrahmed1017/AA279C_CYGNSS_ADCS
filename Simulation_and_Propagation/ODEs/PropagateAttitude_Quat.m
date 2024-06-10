function qw_d = PropagateAttitude_Quat(qw, M_vec, I_p)

% qw --> [q1 q2 q3 q4 w1 w2 w3]

q = qw(1:4);
w = qw(5:7);

% M_contrl calc

% M_actuator torque M_vec = actuator torque

q_d = QuaternionKinematics(q,w);

w_d = PropagateEuler(w, M_vec, I_p);

qw_d = [q_d; w_d];

end

















