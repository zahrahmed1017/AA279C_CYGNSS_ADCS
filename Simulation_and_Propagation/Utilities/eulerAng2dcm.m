function A = eulerAng2dcm(angs)
% angs -> [phi, theta, psi] in radians for 3-1-3 rotation

phi = angs(1); 
theta = angs(2);
psi = angs(3);

A = [cos(psi)*cos(phi) - cos(theta)*sin(psi)*sin(phi), ...
    cos(psi)*sin(phi) + cos(theta)*sin(psi)*cos(phi), ...
    sin(theta)*sin(psi); % row 1
    -sin(psi)*cos(phi) - cos(theta)*cos(psi)*sin(phi), ...
    -sin(psi)*sin(phi) + cos(theta)*cos(psi)*cos(phi), ...
    sin(theta)*cos(psi); % row 2
    sin(theta)*sin(phi), ...
    -sin(theta)*cos(phi), ...
    cos(theta)]; % row 3

end