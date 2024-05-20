classdef BasicMagSensor
    properties
        measNoise
        measBias
        calDay
        gmst % TODO: see if we need to increment this (does it affect our mag model?)
    end

    methods
        function obj = BasicSunSensor(measNoise, measBias, calDay)
            obj.measNoise =  measNoise;
            obj.measBias = measBias;
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

            % for now: model bias as a rotation about body +z direction
            % (arbitrary)
            R_bias = [ cos(obj.measBias), sin(obj.measBias), 0; 
                      -sin(obj.measBias), cos(obj.measBias), 0;
                                       0,                 0, 1];
            sun_vec_p = R_bias * sun_vec_p;

            % add Gaussian noise

            % Old method (didn't work; nonlinear combo of angles --> too
            % much error)
            % comp_noise = sqrt( (self.measNoise^2) / 3) ;
            % noise_this_step = randn(3,1) * comp_noise; % euler angles
            % % R_lin_noise = [1,                   noise_this_step(1), -noise_this_step(2);
            % %                -noise_this_step(1),                  1,  noise_this_step(3);
            % %                 noise_this_step(2), -noise_this_step(3),                 1];
            % R_lin_noise = eulerAng2dcm(noise_this_step); % fyi this func does 3-1-3 rotation
            % 

            % new method: choose an arbitrary axis, rotate around it by an
            % amount sampled from dist. parametrized by measNoise
            rot_ax = rand(3,1);
            rot_ax = rot_ax ./ norm(rot_ax);

            R_noise = axisAngle2dcm(rot_ax, randn(1)*obj.measNoise);

            sun_vec_p = R_noise * sun_vec_p;

            % sanity check: what's the total rotation error now?
            [~, total_noise_error] = dcm2AxisAngle(R_noise);

            meas = sun_vec_p;
        end

    end

end

