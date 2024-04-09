%%% Moment of Inertia Calculations %%%%%%%%

clear;
%% Define Geometry:

% Defining (0, 0, 0) to be at the center of the top face of solar panel.
% Calling this frame A.
% length = along y axis
% width  = along x axis
% height = along z axis
sp_l = 1.6;
sp_w = 0.5;
sp_h = 0.01; 

% Satellite bus is a rectangular prism minus 2 triangular prisms
busRec_l = 0.6;
busRec_w = 0.5;
busRec_h = 0.25;
busTri_l  = (busRec_l - 0.15)/2;
busTri_w  = 0.5;
busTri_h  = 0.25 - 0.15;

% Sensor
sensor_l = 0.15;
sensor_w = 0.3;
sensor_h = 0.01;

%% Calculate C.O.M. coordinates 

% % The y-coordinate of CM is 0 from symmetry. The x- and z- coordinates will
% % be shifted.

% Defining wrt to solar panel frame.
v_sp  = sp_l * sp_w * sp_h;
v_bus = (busRec_l * busRec_w * busRec_h) - (2 * 0.5 * busTri_l * busTri_w * busTri_h);
v_sen = (sensor_l * sensor_w * sensor_h);
v_tot = v_sp + v_bus + v_sen;

m_tot = 25; %kg
m_sp  = m_tot * (v_sp / v_tot);
m_bus = m_tot * (v_bus / v_tot);
m_sen = m_tot * (v_sen / v_tot);


zCM = ((v_sp * sp_h/2) + (v_bus * (sp_h + busRec_h/2)) + ...
      (v_sen * (sp_h + busRec_h + sensor_h/2))) / v_tot;
xCM = (v_sen * (busTri_w/2 - (sensor_w/2))) / v_tot;
yCM = 0;

CM = [xCM, yCM, zCM]; % This is defined wrt the solar panel frame

%% Inertia Tensor calculations - Inertia tensor is with respect to body frame
% Body frame is the same as solar panel frame

% Using parallel axis theorem to find moments of inertia
% Solar Panels
Ixx_sp = ((1/12) * m_sp * (sp_l^2 + sp_h^2)) + (m_sp * xCM^2);
Iyy_sp = (1/12) * m_sp * (sp_w^2 + sp_h^2);
Izz_sp = ((1/12) * m_sp * (sp_l^2 + sp_w^2)) + (m_sp * (zCM - sp_h/2)^2);
% Satellite Bus
Ixx_prime_bus = ((1/12) * m_bus * (busRec_l^2 + busRec_h^2));
Izz_prime_bus = (1/12) * m_bus * (busRec_l^2 + busRec_w^2);
Ixx_bus       = Ixx_prime_bus + (m_bus * xCM^2);
Iyy_bus       = (1/12) * m_bus * (busRec_w^2 + busRec_h^2);
Izz_bus       = Izz_prime_bus + (m_bus * (((busRec_h/2) + sp_h) - zCM)^2);
% Sensor:
Ixx_sen = ((1/12) * m_sen * (sensor_l^2 + sensor_h^2)) + (m_sen * (0.1 - xCM)^2);
Iyy_sen = (1/12) * m_sen * (sensor_w^2 + sensor_h^2);
Izz_sen = ((1/12) * m_sen * (sensor_l^2 + sensor_w^2)) + (m_sen * ((sp_h + busRec_h + sensor_h/2) - zCM)^2);
% Total 
Ixx = Ixx_sp + Ixx_bus + Ixx_sen;
Iyy = Iyy_sp + Iyy_bus + Iyy_sen;
Izz = Izz_sp + Izz_bus + Izz_sen;

% Product of inertia - due to symmetry, the only product of inertia term
% will be a Ixz term.
Ixz_sp  = zCM * xCM * (sp_w * sp_h);
Ixz_bus = ((busRec_h/2) + sp_h - zCM) * (xCM) * (busRec_w * busRec_h);
Ixz_sen = ((sp_h + busRec_h + sensor_h/2) - zCM) * (0.1 - xCM) * (sensor_h * sensor_w);

Ixz = Ixz_sp + Ixz_bus + Ixz_sen;

% Other product of inertia's are zero:
Ixy = 0;
Iyz = 0;

%% Inertia Tensor:

I = [Ixx, Ixy, Ixz;...
     Ixy, Iyy, Iyz;...
     Ixz, Iyz, Izz];

% Checking with CAD output (which is in lb*in^2)
I_lbin = I .* 3410;




