function R = qMethod(weights, all_meas, all_ref)
% weights must be unit norm and dimensions 1xn, where n is the number of measurements
% all_meas must be 3xn (concatenation of all measurements)
% "  " for reference directions

W = [sqrt(weights); sqrt(weights); sqrt(weights)] .* all_meas;
U = [sqrt(weights); sqrt(weights); sqrt(weights)] .* all_ref;

B = W*U';
S = B + B';
Z = [B(2,3)-B(3,2), B(3,1)-B(1,3), B(1,2)-B(2,1)]';
sigma = trace(B);
K = [S-eye(3)*sigma, Z; 
     Z',       sigma ];
[v,d] = eig(K);
[~,ind] = max([d(1,1), d(2,2), d(3,3), d(4,4)]);
q_est = v(:,ind);

R = quat2dcm(q_est([4 1 2 3])');