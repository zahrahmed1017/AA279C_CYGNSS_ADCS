function angs = dcm2eulerAng(A)
% angs -> [phi, theta, psi] in radians for 3-1-3 rotation

phi = -atan2(A(3,1), A(3,2)); % not sure if atan2 is needed, but probably can't hurt

theta = acos(A(3,3));

psi = atan2(A(1,3), A(2,3));

angs = [phi; theta; psi];

end