function angs_d = EulerAngKinematics(angs, w)
% 3-1-3 Euler angle kinematics
% angs -> [phi, theta, psi] in radians
% phi is the initial rotation about axis 3
% theta is the second rotation about axis 1'
% psi is the third rotation about axis z

phi = angs(1);
theta = angs(2);
psi = angs(3);

wx = w(1);
wy = w(2);
wz = w(3);

% kinematics euqaiotns from slides
phi_d = (wx * sin(psi) + wy * cos(psi)) / sin(theta);
theta_d = wx * cos(psi) - wy * sin(psi);
psi_d = wz - (wx * sin(psi) + wy * cos(psi)) / cot(theta);

angs_d = [phi_d; theta_d; psi_d];

end