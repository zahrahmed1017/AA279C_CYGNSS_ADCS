classdef BasicMagSensor
    properties
        measNoise
        measBias
        calDay
        gmst % TODO: see if we need to increment this (does it affect our mag model?)
    end

    methods
        function obj = BasicMagSensor(measNoise, measBias, calDay, gmst)
            obj.measNoise =  measNoise;
            obj.measBias = measBias;
            obj.calDay = calDay;
            obj.gmst = gmst;
        end

        function meas = get_measurement(obj, R_i_p, pos)
            % pos is 1x3 vector of s/c postion in ECI
            % R_i_p is the current attitude matrix, i.e. the rotation from
            % inertial to prinicpal/body frame
            
            % original vector measurement in inertial frame
            [~,B_norm, B_vec_i] = CalculateMagneticTorque(pos, obj.calDay, obj.gmst);
            B_vec_i = B_vec_i/norm(B_vec_i);

            % rotate into body-fixed frame
            B_vec_p = R_i_p * B_vec_i;

            % for now: model bias as a rotation about body +z direction
            % (arbitrary)
            R_bias = [ cos(obj.measBias), sin(obj.measBias), 0; 
                      -sin(obj.measBias), cos(obj.measBias), 0;
                                       0,                 0, 1];
            B_vec_p = R_bias * B_vec_p;

            % add Gaussian noise

            % new method: choose an arbitrary axis, rotate around it by an
            % amount sampled from dist. parametrized by measNoise
            rot_ax = rand(3,1);
            rot_ax = rot_ax ./ norm(rot_ax);

            R_noise = axisAngle2dcm(rot_ax, randn(1)*obj.measNoise);

            B_vec_p = R_noise * B_vec_p;

            % sanity check: what's the total rotation error now?
            [~, total_noise_error] = dcm2AxisAngle(R_noise);

            meas = B_vec_p;

            % TODO: increment gmst?
        end

    end

end

