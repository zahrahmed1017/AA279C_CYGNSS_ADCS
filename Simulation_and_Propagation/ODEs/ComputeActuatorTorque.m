function Lw_d = ComputeActuatorTorque(Lw, Mc, w, A, Astar)

Lw_d = Astar * (-Mc - cross(w, A*Lw));

end