function quaternion = quatMul(q1, q2)

    % Quaternion multiplication, assumes scalar last!!!!
    a1 = q1(4);
    b1 = q1(1);
    c1 = q1(2);
    d1 = q1(3);
    
    a2 = q2(4);
    b2 = q2(1);
    c2 = q2(2);
    d2 = q2(3);
    
    a = a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2;
    b = a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2;
    c = a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2;
    d = a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2;
    
    quaternion = [b; c; d; a]; 

end