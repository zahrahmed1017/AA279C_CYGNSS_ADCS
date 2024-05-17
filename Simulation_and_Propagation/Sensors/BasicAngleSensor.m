classdef BasicAngleSensor
    % In progress, not sure I'll finish
    properties
        measErr
        measBias
    end

    methods
        function obj = BasicAngleSensor(measErr, measBias)
            obj.measErr =  measErr;
            obj.measBias = measBias;
        end

        function meas = get_measurement()
        end

    end

end
   
            