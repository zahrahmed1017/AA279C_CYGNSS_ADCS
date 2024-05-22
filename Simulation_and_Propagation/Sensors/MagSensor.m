classdef MagSensor
    % UNDER CONSTRUCTION
    properties
        % for V = F*B + V0
        voltageBias % V0
        measMatrix  % F
        calDay
        gmst % TODO: see if we need to increment this (does it affect our mag model?)
    end

    methods
        function obj = MagSensor(voltageBias, measMatrix, calDay, gmst)
            obj.voltageBias =  voltageBias;
            obj.measMatrix = measMatrix;
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

            % predict onboard voltage
            V = obj.measMatrix * B_vec_p + obj.voltageBias;


            %TODO: ADC simulation

            meas = B_vec_p; % for now

            % TODO: increment gmst?
        end

    end

end

