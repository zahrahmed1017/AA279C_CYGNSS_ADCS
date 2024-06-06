function Mc = controlTorque_inertial_ver1(I_p, a, a_dot)
% Control torque for inertial pointing, assuming small displacement from
% nominal
% Assume RW triad

% Mc = I * w_dot;

f = 1/10;

kp = (f^2) / I_p(3,3);
kd = 2*sqrt(I_p(3,3) * kp); % zeta = 1

Mc =  5*kp*a + kd*a_dot;

end
