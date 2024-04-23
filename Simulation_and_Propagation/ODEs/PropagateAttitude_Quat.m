function qw_d = PropagateAttitude_Quat(qw, M_vec, I_p)

% qw --> [q1 q2 q3 q4 w1 w2 w3]

% Not sure if I can actually do it this way.....was trying to just make one
% equation to propagate the dynamics (Euler) and then kinematics
% (w quaternions) in one function. But we need to integrate the dynamics at
% every time step and pass that to kinematics?

q = qw(1:4);
w = qw(5:7);

q_d = QuaternionKinematics(q,w);

w_d = PropagateEuler(w, M_vec, I_p);

qw_d = [q_d; w_d];

end

















