function angularVelocity = Analytical_Euler_TF_AS(w0,t, Ixy, Iz)

% Analytical solution to torque free, axially symmetric Euler Eqns
%{
    Inputs:
        - w0  - 3 x 1 initial angular velocity vector
        - t   - m x 1 array of time in sections
        - Ixy - moment of inertia about x and y axis. Note Ixy = Ix = Iy
                for axial-symmetric object
        - Iz  - moment of inertia about z axis.
%}

    wx0 = w0(1);
    wy0 = w0(2);
    wz0 = w0(3);

    wz = wz0;
    A  = ((Ixy - Iz)/Ixy) * wz;

    angularVelocity = zeros(length(t), 3);

    for i = 1:length(t)
        time = t(i);
        wx = (wx0 * cos(A*time)) + (wy0 * sin(A*time));
        wy = (wy0 * cos(A*time)) - (wx0 * sin(A*time));

        angularVelocity(i,:) = [wx, wy, wz];
    end
    
end