function w_dot = PropagateAttitude(w, M_vec, I_p)
% w -> initial angular velocity in PA frame. Units: rad/s
% M_vec -> vector of constant torque values in PA frame. Set to [0, 0, 0]
% for toque-free motion. Units: N-m
% I_p -> Inertia tensor of spacecraft in PA frame (diagonal matrix). Units:
% kg*m^2

% w_dot -> rate of change of angular velocities in PA frame

wx = w(1);
wy = w(2);
wz = w(3);
Ix = I_p(1,1);
Iy = I_p(2,2);
Iz = I_p(3,3);
Mx = M_vec(1);
My = M_vec(2);
Mz = M_vec(3);

wx_d = (Mx - (Iz - Iy)*wz*wy)/Ix;
wy_d = (My - (Ix - Iz)*wx*wz)/Iy;
wz_d = (Mz - (Iy - Ix)*wy*wx)/Iz;


w_dot = [wx_d, wy_d, wz_d]';


end