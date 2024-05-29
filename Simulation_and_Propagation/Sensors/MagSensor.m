classdef MagSensor
    % UNDER CONSTRUCTION
    properties
        % for V = F*B + V0
        voltageBias % V0
        measMatrix  % F
        voltageNoise % W
        calDay
        gmst % TODO: see if we need to increment this (does it affect our mag model?)
    end

    methods
        function obj = MagSensor(voltageBias, voltageNoise, measMatrix, calDay, gmst)
            obj.voltageBias =  voltageBias;
            obj.measMatrix = measMatrix;
            obj.voltageNoise = voltageNoise;
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
            V = obj.measMatrix * B_vec_p ...
                + obj.voltageBias ...
                + obj.voltageNoise * randn(3,1);

            % round V to nearest 100mV
            V = round(V*100*1e6)/(100*1e6);

            % Convert to B vector estimate
            % we subtract out voltagebias becasue we know it; we can't 
            % subtract noise bc we dont' know it
            meas = inv(self.measMatrix)*(V - obj.voltageBias);


        end

    end

end

