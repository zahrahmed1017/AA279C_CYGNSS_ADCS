classdef SunSensor
    
    properties
        angNoise
        angQuant
        calDay
        % gmst % TODO: see if we need to increment this (does it affect our mag model?)
    end

    methods
        function obj = SunSensor(angNoise, angQuant, calDay)
            obj.angNoise =  angNoise;
            obj.angQuant = angQuant;
            obj.calDay = calDay;
        end

        function meas = get_measurement(obj, R_i_p)
            % pos is 1x3 vector of s/c postion in ECI
            % R_i_p is the current attitude matrix, i.e. the rotation from
            % inertial to prinicpal/body frame
            
            % original vector measurement in inertial frame
            [sun_vec_i, ~] = CalculateSunPositionECI(obj.calDay);

            % rotate into body-fixed frame
            sun_vec_p = R_i_p * sun_vec_i;


            % 2-axis sun sensor model from Markely and Crassidis
            n = 1.4553;
            h = 0.448;

            z = (h*sun_vec_p(3)) / sqrt(n^2 - sun_vec_p(1)^2  - sun_vec_p(3)^2);
            % need to check the sign -- should match sign of z? TODO

            x = z * sun_vec_p(1) / sun_vec_p(3);

            phi = atan2(x, z);
            theta = atan2( (n*sqrt(x^2 + z^2)) , sqrt(h^2 - (n^2-1)*(x^2+z^2) )  );
          
            % add some noise, then quantize
            phi = phi + randn(1)*obj.angNoise;
            theta = theta + randn(1)*obj.angNoise;

            phi = round(phi/obj.angQuant) * obj.angQuant;
            theta = round(theta/obj.angQuant) * obj.angQuant;

            tan_a = tan(theta)*sin(phi);
            tan_b = tan(theta)*cos(phi);

            meas = (1 / sqrt(1 + tan_a^2 + tan_b^2)) * [tan_a; 1; tan_b];

            % check that the sign is right
            % if dot(meas, sun_vec_p) < 0
            %     meas = - meas;
            % end

            % check sign for each component
            % assume low enough error that this is ok
            for comp=1:3
                if sign(meas(comp)) ~= sign(sun_vec_p(comp))
                    meas(comp) = - meas(comp);
                end
            end


        end

    end

end

