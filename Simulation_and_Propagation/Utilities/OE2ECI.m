function eci = OE2ECI(oe, MU)
% Take orbital elements and compute position and velocity in ECI frame
%
% Inputs: 
%           oe - 6xN orbit element vectors consisting of
%               a  - semi-major axis (km)
%               e  - eccentricity
%               i  - inclination (radians)
%               O  - right ascension of ascending node (radians)
%               w  - argument of periapsis (radians)
%               va - true anomaly (radians)
%           MU - gravitational parameter of planet (km^3/s^2)
% Outputs:
%           eci - 6xN ECI state vectors consisting of three position components and
%                 three velocity components (km, km, km, km/s, km/s, km/s)

% Number of input vectors
N = size(oe, 2);
% Tolerance for integer rounding for e and i
eps = 1E-5;

% Check for valid inputs
if size(oe, 1) == 6 && MU > 0

    % Size output
    eci = zeros(6, N);

    for ii = 1:N

        % Get individual orbit elements
        a  = oe(1, ii);
        e  = oe(2, ii);
        i  = wrapToPi(oe(3, ii));
        O  = wrapTo2Pi(oe(4, ii));
        w  = wrapTo2Pi(oe(5, ii));
        va = wrapTo2Pi(oe(6, ii));

        % Test for valid inputs
        if (e >= 0 && e <= 1 && a >= 0)

            % Test for edge cases
            if (e < eps)
                % Circular orbit
                if (i < eps) || (abs(i - pi) < eps)
                    % Circular equatorial orbit
                    w = 0;
                    O = 0;
                    va = O + w + va;
                else
                    % Circular inclined orbit
                    w = 0;
                    va = w + va;
                end
            else
                if (i < eps) || (abs(i - pi) < eps)
                    % Elliptical equatorial orbit
                    O = 0;
                    w = O + w; 
                end
            end

            % Test for parabolic orbit
            if (abs(e - 1) < eps)
                disp('Parabolic orbit: please input rp in place of a.')
                p = 2 * a;
            else
                p = a * (1 - e^2);
            end
            % Compute radius and velocity of orbit
            r = p / (1 + e * cos(va));
            v = sqrt(MU / p);

            % Find position and velocity in perifocal coordinates
            rpqw = r * [cos(va); sin(va); 0];
            vpqw = v * [-sin(va); e + cos(va); 0];

            % Define rotation matrix from perifocal to ECI (passive rotations)
            % R = Rz(-rasc)Rx(-inc)Rz(-argp)
            R1 = [cos(-O) sin(-O) 0;...
                 -sin(-O) cos(-O) 0;...
                  0       0       1];
            R2 = [1  0        0;...
                  0  cos(-i) sin(-i);...
                  0 -sin(-i) cos(-i)];
            R3 = [cos(-w) sin(-w) 0;...
                 -sin(-w) cos(-w) 0;...
                  0        0        1];
            R = R1 * R2 * R3;
            
            % Find position and velocity in ECI coordinates
            eci(:, ii) = [R * rpqw; R * vpqw];

        else
            eci(:, ii) = nan(6, 1);
            disp('Invalid orbit elements in OE2ECI');
        end
    
    end
else
    eci = nan(6, N);
    disp('Invalid inputs to OE2ECI');
end