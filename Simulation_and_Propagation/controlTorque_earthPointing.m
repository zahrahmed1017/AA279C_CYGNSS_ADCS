function Mc = controlTorque_earthPointing(I_p, a, a_dot, n)
% Control torque for earth pointing, assuming small displacement from
% nominal
% Assume we want z aligned w/ R, x aligned w/ T
% n is mean motion

% Mc = I * w_dot;

% f = 1/20;
f = 1/20;

k_p = zeros(3,1);
k_d = zeros(3,1);

% x
k_p(1) = (f^2)/I_p(1,1) - (n^2)*(I_p(2,2) - I_p(3,3));
k_d(1) = 2*sqrt(I_p(1,1) * k_p(1));

% y
k_p(2) = (f^2)/I_p(2,2);
k_d(2) = 2*sqrt(I_p(2,2) * k_p(2));

% z
k_p(3) =  (f^2)/I_p(3,3) - (n^2)*(I_p(2,2) - I_p(1,1));
k_d(3) = 2*sqrt(I_p(3,3) * k_p(3));


Mc =  k_p.*a + k_d.*a_dot;

end
