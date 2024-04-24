function quat = dcm2quaternion(A)

quat_temp = dcm2quat(A);
quat      = quat_temp([2;3;4;1])';


% based on backup slides from lecture 5
% 4th quaternion element is the scalar

% we determine 3 other quat entries based on whichever entry is largest.
% See: lecture slides, Sheppard method
% https://arc.aiaa.org/doi/epdf/10.2514/3.55767b 

% T = trace(A);
% 
% p_1_2 = 1 + 2*A(1,1) - T; % p_1^2
% p_2_2 = 1 + 2*A(2,2) - T;
% p_3_2 = 1 + 2*A(3,3) - T;
% p_4_2 = 1 + T;
% 
% if p_1_2 == max([p_1_2, p_2_2, p_3_2, p_4_2])
% 
%     p_1 = sqrt(p_1_2);
%     p_2 = (A(2,1) + A(1,2))/p_1;
%     p_3 = (A(1,3) + A(3,1))/p_1;
%     p_4 = (A(3,2) - A(2,3))/p_1;
% 
% elseif p_2_2 == max([p_1_2, p_2_2, p_3_2, p_4_2])
% 
%     p_2 = sqrt(p_2_2);
%     p_1 = (A(2,1) + A(1,2)) / p_2;
%     p_3 = (A(3,2) + A(2,3)) / p_2;
%     p_4 = (A(1,3) - A(3,1)) / p_2;
% 
% elseif p_3_2 == max([p_1_2, p_2_2, p_3_2, p_4_2])
% 
%     p_3 = sqrt(p_3_2);
%     p_1 = (A(1,3) + A(3,1)) / p_3;
%     p_2 = (A(3,2) + A(2,3)) / p_3;
%     p_4 = (A(2,1) - A(1,2)) / p_3;
% 
% else % p_4_2 == max([p_1_2, p_2_2, p_3_2, p_4_2])
% 
%     p_4 = sqrt(p_4_2);
%     p_1 = (A(3,2) - A(2,3)) / p_4;
%     p_2 = (A(1,3) - A(3,1)) / p_4;
%     p_3 = (A(2,1) - A(1,2)) / p_4;
% 
% end
% 
% quat = [p_1; p_2; p_3; p_4] / 2;

end