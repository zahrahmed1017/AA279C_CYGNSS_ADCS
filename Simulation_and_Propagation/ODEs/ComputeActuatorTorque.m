function Lw_d = ComputeActuatorTorque(Lw, Mc, w, A, Astar)

    if size(w,2) ~= 1
        w = w';
    end
    
    if size(Lw, 2) ~= 1
        Lw = Lw';
    end


Lw_d = Astar * (-Mc - cross(w, A*Lw));

end