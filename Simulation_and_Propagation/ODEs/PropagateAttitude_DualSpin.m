function qw_d = PropagateAttitude_DualSpin(qw, M_vec, I_p, M_r, I_r, r_rot)

% qw --> [q1 q2 q3 q4 w1 w2 w3 wr]

q = qw(1:4);
w = qw(5:8);

q_d = QuaternionKinematics(q,w);

w_d = PropagateEulerDualSpin(w, M_vec, I_p, M_r, I_r, r_rot);

if size(w_d, 2) ~= 1;
    w_d = w_d';
end

if size(q_d, 2) ~= 1;
    q_d = q_d';
end

qw_d = [q_d; w_d];

end