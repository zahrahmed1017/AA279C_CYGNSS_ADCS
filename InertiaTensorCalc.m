%%% Moment of Inertia Calculations %%%%%%%%

clear; clc;
%% Define Geometry:

% Defining (0, 0, 0) to be at the center of the top face of solar panels
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

sp_area_yz     = (sp_l * sp_h);
busRec_area_yz = (busRec_l * busRec_h);
busTri_area_yz = (0.5 * busTri_l * busTri_h);
sensor_area_yz = (sensor_l * sensor_h);
total_area_yz  = sp_area_yz + busRec_area_yz + sensor_area_yz - (2 * busTri_area_yz);

zCM   = ((sp_area_yz * sp_h/2)...
    + (busRec_area_yz * (sp_h + busRec_h/2))...
    + (sensor_area_yz * (sp_h + busRec_h + sensor_h/2))...
    - (2 * busTri_area_yz * (sp_h + 0.15 + (2/3 * busTri_h))))/total_area_yz;


sp_area_xy      = (sp_l * sp_w);
busRec_area_xy  = (busRec_l * busRec_w);
busTri_area_xy  = (busTri_l * busTri_w);
sensor_area_xy  = (sensor_l * sensor_w);
total_area_xy   = sp_area_xy + busRec_area_xy + sensor_area_xy - (2 * busTri_area_xy);

xCM = (sensor_area_xy * (busTri_w/2 - (sensor_w/2)))/total_area_xy;
yCM = 0;

CM = [xCM, yCM, zCM];

%% Inertia Tensor calculations

% Masses of components:
m_tot = 25; %kg
v_sp  = sp_l * sp_w * sp_h;
v_bus = (busRec_l * busRec_w * busRec_h) - (2 * 0.5 * busTri_l * busTri_w * busTri_h);
v_sen = (sensor_l * sensor_w * sensor_h);
v_tot = v_sp + v_bus + v_sen;

m_sp  = m_tot * (v_sp / v_tot);
m_bus = m_tot * (v_bus / v_tot);
m_sen = m_tot * (v_sen / v_tot);

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

% Checking with CAD output (which is in lb*in^2)
Ixx_lbin = Ixx * 3410;
Iyy_lbin = Iyy * 3410;
Izz_lbin = Izz * 3410;

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


