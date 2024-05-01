function rtn_series = rv2rtn(rv_series)
% takes in a time series of position, velocity info as n x 6 matrix. 
% Outputs a time series of R, T, N unit vector positions as n x 9 matrix.
[n, ~] = size(rv_series);

% Calculating the R-T-N unit vectors:
R_unitvec = zeros(n,3);
T_unitvec = zeros(n,3);
N_unitvec = zeros(n,3);

for i = 1:n
    rx = rv_series(i,1);
    ry = rv_series(i,2);
    rz = rv_series(i,3);
    vx = rv_series(i,4);
    vy = rv_series(i,5);
    vz = rv_series(i,6);
    r_norm = sqrt(rx^2 + ry^2 + rz^2);
    v_norm = sqrt(vx^2 + vy^2 + vz^2);
    h      = cross([rx,ry,rz], [vx, vy, vz]);
    h_norm = sqrt(h(1)^2 + h(2)^2 + h(3)^2);

    R_i = [rx/r_norm, ry/r_norm, rz/r_norm];
    N_i = [h(1)/h_norm, h(2)/h_norm, h(3)/h_norm];
    T_i = cross(N_i, R_i);

    % T_i = [vx/v_norm, vy/v_norm, vz/v_norm];
    % N_i = cross(R_i,T_i);

    R_unitvec(i,:) = R_i;
    T_unitvec(i,:) = T_i;
    N_unitvec(i,:) = N_i;
    
end


rtn_series = [R_unitvec, T_unitvec, N_unitvec]; % concatenate horizontally -> n x 9




end