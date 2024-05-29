classdef MEKF

    properties (Access = public)
        % properties
        state           % State Vector [q0 q1 q2 q3 w1 w2 w3]^T !!!! scalar first in the quat
        errorState      % Error state using MRP [p1 p2 p3 w1 w2 w3 w4]^T
        stateCov        % State covariance
        processNoise    % Process noise
        quatKin         % Linearized, discrete-time quaternion kinematics

        % Other parameters
        dt              % Discrete time step
        I               % Principal Inertia Tensor

        % STM
        A11             % dp_dot/dp
        A12             % dp_dot/dw
        A22             % dw_dot/dw
        STM
    end


    methods 
        %%% Constructor
        function obj = MEKF(initialState, initialStateCov, initialProcessNoise, dt, I)

            if size(initialState,2) ~= 1
                initialState = initialState';
            end

            obj.state         = initialState;
            obj.errorState    = [0; 0; 0; initialState(5:7)];
            obj.stateCov      = initialStateCov;
            obj.processNoise  = initialProcessNoise; 
            obj.quatKin       = [];
            obj.dt            = dt;
            obj.I             = I;




        end

        %%% MEKF Steps:
        function obj = Propagate(obj)
            
            obj.QuatKin_Linear_Discrete();
            obj.PropagateQuat();
            obj.errorState(1:3) = [0; 0; 0];

        end

        function obj = TimeUpdate(obj)

            obj.UpdateSTM();
            obj.UpdateAngularVelocity();
            obj.UpdateCovariance()

        end

        function obj = MeasurementUpdate(obj)


        end

        function obj = Reset(obj)
            % Push error state to absolute state
            % MRP -> Quat
            % angular velocity in angular state transfers to absolute state
            dcm_old = quat2dcm(obj.state(1:4)');
            dcm_err = rod2dcm(obj.errorState(1:3)');
            dcm_new = dcm_err * dcm_old;
            q_new   = dcm2quat(dcm_new);

            obj.state(1:4) = q_new';
            obj.state(5:7) = obj.errorState(4:6);
            

        end

        %%% Subfunctions:
        function obj = QuatKin_Linear_Discrete(obj)

            w = obj.errorState(4:6);
            w1 = w(1);
            w2 = w(2);
            w3 = w(3);
            skewSym      = [0,  w3, -w2, w1;...
                           -w3, 0,   w1, w2;...
                            w2, -w1, 0,  w3;...
                           -w1, -w2, -w3, 0];
            obj.quatKin = eye(4,4) + (0.5).*skewSym.*obj.dt;

        end

        function obj = PropagateQuat(obj)
            
            Omega = obj.quatKin;
            qt_1  = obj.state(1:4);
            qt    = Omega * qt_1;

            obj.state(1:4) = qt;

        end

        function obj = UpdateAngularVelocity(obj)

            obj.errorState(4:6) = obj.A22 * obj.errorState(4:6);
            
        end

        function obj = UpdateCovariance(obj)
        
            obj.stateCov = obj.STM * obj.stateCov * obj.STM' + obj.processNoise; 

        end

        %%% Get Methods:
        function stm = get.STM(obj)
        
            stm = [obj.A11, obj.A12; ...
                       zeros(3,3), obj.A22];

        end 

        function a11 = get.A11(obj)
            
            if size(obj.errorState,2) ~= 1
                obj.errorState = obj.errorState';
            end

            wx = crossMatrix(obj.errorState(4:6));
            a11 = eye(3,3) + (1/2) * obj.dt .* wx ; 

        end

        function a12 = get.A12(obj)
            
            a12 = (1/4) * obj.dt .* eye(3,3);

        end

        function a22 = get.A22(obj)
            
            Ix = obj.I(1,1);
            Iy = obj.I(2,2);
            Iz = obj.I(3,3);
            wx = obj.errorState(4);
            wy = obj.errorState(5);
            wz = obj.errorState(6);

            eulEqLin = [ 0,                 (Iy - Iz)/Ix * wx, (Iy - Iz)/Ix * wy;...
                        (Iz - Ix)/Iy * wz,   0,                (Iz - Ix)/Iy * wx;...
                        (Ix - Iy)/Iz * wy,   (Ix - Iy)/Iz * wx, 0];
            a22  = eye(3,3) + obj.dt .* eulEqLin; 

        end

        %%% Helper Functions
        function vx = crossMatrix(v)
            % v is a 3x1 column matrix
            % outputs the [vx] matrix
            
            vx = [    0,  -v(3),  v(2);
                   v(3),      0, -v(1);
                  -v(2),   v(1),     0];

        end


    end



end
        
