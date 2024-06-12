function M_mag = ComputeMagnetorquerTorque(dH, B, kmag)
% Lsat --> momentum of whole satellite
% B --> magnetic field vector
% dH --> momentum error, i.e. what you want to correct (vector)
% kmag --> magnetorquer gain (not sure if this should be here)

b = B ./ norm(B);

m = kmag * (cross(dH, b));

M_mag = cross(m, B);

end