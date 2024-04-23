function q_dot = QuaternionKinematics(q, w)

q_norm = q ./ (norm(q));

q1 = q_norm(1);
q2 = q_norm(2);
q3 = q_norm(3);
q4 = q_norm(4); % SCALAR

% Skew symmetric matrix
wx = w(1);
wy = w(2);
wz = w(3);

Om = [0,  wz, -wy,  wx;...
     -wz, 0,   wx,  wy;...
      wy, -wx, 0,   wz;...
     -wx, -wy, -wz, 0];

q_dot = 0.5 .* Om * [q1;q2;q3;q4];


end