function q_dot = QuaternionKinematics(q, w)

q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4); % SCALAR

% Skew symmetric matrix
wx = w(1);
wy = w(2);
wz = w(3);

Om = [0,  wz, -wy,  wx;...
     -wz, 0,   wx,  wy;...
      wy, -wx, 0,   wz;...
     -wx, -wy, -wz, 0];

q_dot_temp = 0.5 .* Om * [q1;q2;q3;q4];

% Normalize the quaternion to unit length
q_dot = q_dot_temp ./ norm(q_dot_temp);

end