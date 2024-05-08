function geoc = ECEF2GEOC(ecef, RP)
% Convert ECEF position to geocentric coordinates
%
% Inputs:
%           ecef - 3xN ECEF position vectors (km)
%           RP - radius of planet (km)
% Outputs:
%           geoc - 3xN geocentric position (latitude (deg), longitude (deg), height (km))

% Number of input vectors
N = size(ecef, 2);

% Check for correct input dimension
if size(ecef, 1) == 3
    
    % Size output matrix
    geoc = zeros(3, N);
    
    for ii = 1:N
        % Compute latitude and longitude and height
        geoc(1, ii) = asind(ecef(3, ii) / norm(ecef(:, ii)));
        geoc(2, ii) = atan2d(ecef(2, ii), ecef(1, ii));
        geoc(3, ii) = norm(ecef(:, ii)) - RP;
    end
    
else
    geoc = nan(3, N);
    disp('Wrong input dimension for ECEF2GEOC')
end