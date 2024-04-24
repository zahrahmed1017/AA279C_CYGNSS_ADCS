function dcm_w_d = PropagateAttitude_DCM(dcm_w, M_vec, I_p)

% dcm_w = [dcm(:,1); dcm(:,2); dcm(:,3); w1; w2; w3]

dcm = [dcm_w(1:3), dcm_w(4:6), dcm_w(7:9)];
w   = dcm_w(10:12);

dcm_d = DCMKinematics(dcm,w);
w_d   = PropagateEuler(w, M_vec, I_p);

dcm_w_d = [dcm_d(:,1); dcm_d(:,2); dcm_d(:,3); w_d];

end