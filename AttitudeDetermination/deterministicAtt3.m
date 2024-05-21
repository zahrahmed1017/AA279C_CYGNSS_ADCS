function R = deterministicAtt3(m1, m2, m3, v1, v2, v3)
% Determinstic attitude determination algorithm using three measurements
% all v and m must be 3-element column vectors
% m -> measurement (in a body-fixed frame)
% v -> reference vector (in inertial frame)

M1 = [m1, m2, m3];
V1 = [v1, v2, v3];

R = M1 * inv(V1);

end