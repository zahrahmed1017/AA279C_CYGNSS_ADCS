function M = E2M(E, e)
% Converts eccentric anomaly to mean anomaly
%
% Inputs:
%           E - eccentric anomaly (radians)
%           e - eccentricity 
% 
% Outputs:
%           M - mean anomaly (radians)

% Check for valid eccentricity
if (e >= 0 && e <= 1)
    E = wrapTo2Pi(E);
    M = E - e * sin(E);
    M = wrapTo2Pi(M);
else
    M = E;
    disp('Invalid eccentricity for E2M');
end