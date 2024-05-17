function R = deterministicAtt2(m1, m2, v1, v2)
% Determinstic attitude determination algorithm using two measurements
% all v and m must be 3-element column vectors
% m -> measurement (in a body-fixed frame)
% v -> reference vector (in inertial frame)

p_p = m1;
q_p = cross(m1, m2) / norm(cross(m1, m2));
r_p = cross(p_p, q_p);

p_i = v1;
q_i = cross(v1, v2) / norm(cross(v1, v2));
r_i = cross(p_i, q_i);

V2 = [p_i, q_i, r_i];
M2 = [p_p, q_p, r_p];

R = M2 * inv(V2);

end