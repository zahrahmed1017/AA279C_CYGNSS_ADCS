function A = quat2dcm(quat)
% quat --> [q1 q2 q3 q4], where q4 is the scalar

q4 = quat(4);
q = quat(1:3);

qx = crossMatrix(q);

A = (q4^2 - q' * q) * eye(3) + 2*q*q' - 2*q4*qx;

end