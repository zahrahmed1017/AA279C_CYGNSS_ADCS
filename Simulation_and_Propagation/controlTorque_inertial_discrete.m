function Mc = controlTorque_inertial_discrete(I_p, a, a_dot)
% Control torque for inertial pointing, assuming small displacement from
% nominal
% Assume RW triad

% Mc = I * w_dot;

% f = 1/20;
f = 1/20;

k_p = zeros(3,1);
k_d = zeros(3,1);

for i=1:3
    k_p(i) = (f^2) / I_p(i,i);
    k_d(i) = 2*sqrt(I_p(i,i) * k_p(i)); % zeta = 1
end

Mc =  k_p.*a + k_d.*a_dot;



end
