function GMST = CAL2GMST(Y1, M1, D1, D)
% Convert calendar time to GMST time
% Inputs:
%           Y1 - year
%           M1 - month
%           D1 - day
%           D - fractional day
%   Outputs:
%           GMST - GMST time (radians)

% Compute modified month and year
if M1 <= 2
    Y2 = Y1 - 1;
    M2 = M1 + 12;
else
    Y2 = Y1;
    M2 = M1;
end
B = Y1 / 400 - Y1 / 100 + Y1 / 4;

% Decimal days
D2 = D1 + D;
% Modified Julian Date
MJD = 365 * Y2 - 679004 + floor(B) + floor(30.6001 * (M2 + 1)) + D2;
d = MJD - 51544.5;

% GMST in degrees
GMST = 280.4606 + 360.9856473 * d;
% Convert to radians
GMST = GMST * pi / 180;
GMST = wrapTo2Pi(GMST);

end