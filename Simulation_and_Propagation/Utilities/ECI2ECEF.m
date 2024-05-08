function ecef = ECI2ECEF(eci, GMST)
% Convert ECI position to ECEF position
%
% Inputs:
%           eci - 3xN ECI position vectors (km)
%           GMST - 1xN list of GMST (radians)
% Outputs:
%           ecef - 3xN ECEF position vectors (km)

% Number of input vectors
N = size(eci, 2);

% Check for correct input dimension
if (size(eci, 1) == 3) && (length(GMST) == N)
    
    % Size output
    ecef = zeros(3, N);
    
    for ii = 1:N
        % Compute rotation from ECI to ECEF
        R_eci2ecef = [cos(GMST(ii)),  sin(GMST(ii)), 0; ...
                      -sin(GMST(ii)), cos(GMST(ii)), 0; ...
                      0,              0,             1];

        % Compute ECEF position
        ecef(:, ii) = R_eci2ecef * eci(:, ii);
    end
    
else
    ecef = nan(3, N);
    disp('Wrong input dimension for ECI2ECEF');
end
    
