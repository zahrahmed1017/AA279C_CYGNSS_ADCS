function w_dot = PropagateEulerDualSpin(w, M_vec, I_p, M_r, I_r, r_rot)

%{
    Description: 
        This function describes the Euler equations for a dual-spin
        satellite (i.e. with a rotor)
    Inputs: 
        - w    : angular velocity vector [wx, wy, wz, wr] where wr is the
                 rotor angular velocity
        - M_vec: vector of applied torques (constant) [Mx, My, Mz];
        - I_p  : Principal moments of inertia [Ix, Iy, Iz]; 
        - M_r  : applied rotor torque
        - I_r  : moment of inertia of the rotor
        - r_rot: unit vector for the rotation axis of the rotor in
                 principal coordinates

%}
    % Angular velocity
    wx = w(1);
    wy = w(2);
    wz = w(3);
    wr = w(4);
    % Moment of inertia
    Ix = I_p(1,1);
    Iy = I_p(2,2);
    Iz = I_p(3,3);
    Ir = I_r;
    % External Torque
    Mx = M_vec(1);
    My = M_vec(2);
    Mz = M_vec(3);
    Mr = M_r;
    % Rotation axis components for rotor
    rx = r_rot(1);
    ry = r_rot(2);
    rz = r_rot(3);

    wr_d = Mr / Ir;
    wx_d = (Mx - (Ir * wr_d * rx) - ((Iz -Iy)*wy*wz) - (Ir * wr * (wy*rz - wz*ry))) / Ix;
    wy_d = (My - (Ir * wr_d * ry) - ((Ix -Iz)*wz*wx) - (Ir * wr * (wz*rx - wx*rz))) / Iy;
    wz_d = (Mz - (Ir * wr_d * rz) - ((Iy -Ix)*wx*wy) - (Ir * wr * (wx*ry - wy*rx))) / Iz;

    w_dot = [wx_d, wy_d, wz_d, wr_d];

end