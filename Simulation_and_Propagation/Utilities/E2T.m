function v = E2T(E, e)
% Convert eccentric anomaly to true anomaly
%
% Inputs: 
%           E - eccentric anomaly (radians)
%           e - eccentricity
% Outputs:
%           v - true anomaly (radians)

E = wrapTo2Pi(E);

% Check for valid eccentricity
if (e >= 0 && e <= 1)

    % Compute true anomaly
    sinv = sin(E) * sqrt(1 - e^2) / (1 - e * cos(E));
    cosv = (cos(E) - e) / (1 - e * cos(E));
    v = atan2(real(sinv), real(cosv));
    v = wrapTo2Pi(v);

else
    v = E;
    disp('Invalid eccentricity for E2T')
end