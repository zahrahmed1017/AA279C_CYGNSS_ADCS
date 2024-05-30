function quat = mrp2quat(mrp)

%%%% SCALAR LAST

sigma = norm(mrp);

quat = [2 * mrp(1); 2 * mrp(2); 2 * mrp(3) ; 1 - sigma^2];
quat = quat / (1 + sigma^2);




end 