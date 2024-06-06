function state_d = PropagateAttitude_Quat_wControl(state, t, I_p, A, Astar)

% state --> [q1 q2 q3 q4 w1 w2 w3 Lw1 Lw2 Lw3]

q  = state(1:4);
w  = state(5:7);
Lw = state(8:10);

% M_contrl calc
M1 = sin(0.1*t);
M2 = sin(0.2*t);
M3 = cos(0.2*t);
Mc = [M1; M2; M3];

% M_actuator torque M_vec = actuator torque
Lw_d = ComputeActuatorTorque(Lw, Mc, w, A, Astar);
Lw_pa = A * Lw_d;

% Propagate Kinematics
q_d = QuaternionKinematics(q,w);

% Propagate Dynamics
Mvec = [0; 0; 0];
w_d = PropagateEuler(w, Mvec, I_p);

% Output State
state_d = [q_d; w_d; Lw_d];

end