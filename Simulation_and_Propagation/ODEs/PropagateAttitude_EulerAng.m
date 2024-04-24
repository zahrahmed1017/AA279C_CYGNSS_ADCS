function ang_w_d = PropagateAttitude_EulerAng(ang_w, M_vec, I_p)

% ang_w --> [phi theta psi w1 w2 w3]
% euler angles correspond to 3-1-3 sequence, in radians

angs = ang_w(1:3);
w    = ang_w(4:6);

% if abs(ang_w(3)) < 0.0001
%     angs_d = [NaN; NaN; NaN];
% else
%     angs_d = EulerAngKinematics(angs,w);
% end
angs_d = EulerAngKinematics(angs,w);
w_d = PropagateEuler(w, M_vec, I_p);

ang_w_d = [angs_d; w_d];

end

















