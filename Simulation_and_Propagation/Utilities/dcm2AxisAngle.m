function [axis, angle] = dcm2AxisAngle(A)

    angle = acos( .5 * (trace(A) - 1) );
    
    e1 = (A(2,3) - A(3,2)) / (2*sin(angle));
    e2 = (A(3,1) - A(1,3)) / (2*sin(angle));
    e3 = (A(1,2) - A(2,1)) / (2*sin(angle));

    axis = [e1; e2; e3];

end