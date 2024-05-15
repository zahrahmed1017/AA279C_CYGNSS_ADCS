function E = M2E(M, e, eps)
% Convert mean anomaly in radians to eccentric anomaly in radians
% Inputs:
%           M - mean anomaly (radians)
%           e - eccentricity
%           eps - error tolerance
% Outputs:
%           E - eccentric anomaly (radians)

% Wrap input to pi
M = wrapTo2Pi(M);

% Initial guess for E
if (e >= 0.8)
	Eold = pi;
	E = pi;
else
	Eold = M;
	E = M;
end

% Initial error
error = eps + 1;
% Counter
count = 1;

% Check for valid e
if (e >= 0 && e <= 1 && eps > 0)

    % If M = 0 or pi, we know E = 0 or pi also
    if (M < eps) || (abs(M - 2 * pi) < eps)
        E = 0;
    elseif (abs(M - pi) < eps)
        E = pi;
    else
        % Otherwise loop the algorithm until tolerance is met
        while (error > eps) && (count < 100)
            % Compute new E value
            Enew = Eold - (Eold - e * sin(Eold) - M) / (1 - e * cos(Eold));
            % Compute change in E
            error = abs(Enew - Eold);
            % Remember computed E
            Eold = Enew;
            % Increase counter
            count = count + 1;
        end
        % Get final E for output
        E = Enew;
		% Check for convergence
		if count == 100
			disp('Convergence problem in M2E');
		end
    end
    
else
    disp('Invalid inputs for M2E')
end

        
        
