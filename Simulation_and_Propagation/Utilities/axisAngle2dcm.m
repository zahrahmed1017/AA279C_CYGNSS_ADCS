function A = axisAngle2dcm(axis, angle)
% angle must be in radians

A = cos(angle)*eye(3) + (1 - cos(angle))*axis*axis' - sin(angle)*crossMatrix(axis) ;

end