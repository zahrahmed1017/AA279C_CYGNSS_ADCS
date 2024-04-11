%%% Inertia Tensor Calculations %%%%%%%%

clear; 

%% Define Reference Frames:
%{
1. Solar Panel Frame - origin at the center of the top frace of the solar
panel. This frame is what I'm basing all of my measurements off of.
2. Satellite Bus Frame - origin at the center of the satellite bus if it
were a rectangle.
3. Sensor frame - origin at center of sensor panel
4. Body Frame - origin at center of mass
** ALL frames are using the same orientation of axes, i.e. they are all
only offset by translation. Y axis is along the length of solar panel, Z
axis is positive nadir pointing and X axis is out the front of the
spacecraft completing the triad.
%}
%% Define Geometry, Volume and Masses:

% Defining (0, 0, 0) to be at the center of the top face of solar panel.
% Calling this the solar panel frame.
% length = along y axis
% width  = along x axis
% height = along z axis
sp_l = 1.6;
sp_w = 0.5;
sp_h = 0.01; 

% Satellite bus is a rectangular prism minus 2 triangular prisms
busRec_l     = 0.6;
busRec_w     = 0.5;
busRec_h     = 0.25;
busRec_side1 = 0.1;
bus_bottom   = 0.15;
busTri_l     = (busRec_l - bus_bottom)/2;
busTri_w     = 0.5;
busTri_h     = busRec_h - busRec_side1;

% Sensor
sensor_l = 0.15;
sensor_w = 0.3;
sensor_h = 0.01;

% Defining wrt to solar panel frame.
v_sp  = sp_l * sp_w * sp_h;
v_bus = (busRec_l * busRec_w * busRec_h) - (2 * 0.5 * busTri_l * busTri_w * busTri_h);
v_rec = busRec_l * busRec_w * busRec_h;
v_tri = 0.5 * busTri_l * busTri_w * busTri_h;
v_sen = sensor_l * sensor_w * sensor_h;
v_tot = v_sp + v_bus + v_sen;

m_tot = 27.96; %kg
m_sp  = m_tot * (v_sp / v_tot);
m_bus = m_tot * (v_bus / v_tot);
m_rec = m_tot * (v_rec / v_tot);
m_tri = m_tot * (v_tri / v_tot);
m_sen = m_tot * (v_sen / v_tot);

% Sanity check 
m_tot_check = m_bus + m_sp + m_sen;

%% Calculate Center of Mass Coordinates:

zCM = ((v_sp * sp_h/2) +...
       (v_rec * (sp_h + busRec_h/2)) + ...
       (v_sen * (sp_h + busRec_h + sensor_h/2)) - ...
       (2 * v_tri * (sp_h + busRec_side1 + (2/3 * busTri_h)))) / v_tot;
xCM = (v_sen * (busTri_w/2 - (sensor_w/2))) / v_tot;
yCM = 0;

CM = [xCM, yCM, zCM]; % This is defined wrt the solar panel frame

%% Inertia Tensor Calculations:

% Solar Panel Tensor:
Ixx_sp = (1/12) * m_sp * (sp_l^2 + sp_h^2);
Iyy_sp = (1/12) * m_sp * (sp_w^2 + sp_h^2);
Izz_sp = (1/12) * m_sp * (sp_l^2 + sp_w^2);
I_sp = [Ixx_sp, 0, 0; 0, Iyy_sp, 0; 0, 0, Izz_sp];
% Translation vector from sp frame to body frame:
T_sp_b = [xCM, yCM, zCM - sp_h/2];
% Solar Panel Inertia Tensor about CM
I_sp_b = ParallelAxisTheorem(I_sp, T_sp_b, m_sp);

% Satellite Bus: 
Ixx_rec = (1/12) * m_rec * (busRec_l^2 + busRec_h^2);
Iyy_rec = (1/12) * m_rec * (busRec_w^2 + busRec_h^2);
Izz_rec = (1/12) * m_rec * (busRec_l^2 + busRec_w^2);
I_rec   = [Ixx_rec, 0, 0; 0, Iyy_rec, 0; 0, 0, Izz_rec];
% % Subtract out the triangles:
Ixx_tri   = (m_tri/18) * (busTri_l^2 + busTri_h^2);
Iyy_tri   = (m_tri/36) * ((3*busTri_w^2) + 2*busTri_h^2);
Izz_tri   = (m_tri/36) * ((3*busTri_w^2) + 2*busTri_l^2);
I_tri     = [Ixx_tri, 0, 0; 0, Iyy_tri, 0; 0, 0, Izz_tri];
T_tri_rec = [0, busRec_l/2 - busTri_l/3, -busRec_h/2 + busTri_h/3];
T_tri_rec2 = [0, -busRec_l/2 + busTri_l/3, -busRec_h/2 + busTri_h/3];
I_tri_rec = ParallelAxisTheorem(I_tri, T_tri_rec, m_tri);
I_tri_rec2 = ParallelAxisTheorem(I_tri, T_tri_rec2, m_tri);
I_bus     = I_rec -  I_tri_rec - I_tri_rec2;

% Translation vector from bus frame to body frame:
bus_frame_origin = [0, 0, sp_h + busRec_h/2];
T_bus_b = bus_frame_origin - CM;
% Bus Inertia Tensor about CM:
I_bus_b = ParallelAxisTheorem(I_bus, T_bus_b, m_bus);

% Sensor:
Ixx_sen = (1/12) * m_sen * (sensor_l^2 + sensor_h^2);
Iyy_sen = (1/12) * m_sen * (sensor_w^2 + sensor_h^2);
Izz_sen = (1/12) * m_sen * (sensor_l^2 + sensor_w^2);
I_sen   = [Ixx_sen, 0, 0; 0, Iyy_sen, 0; 0, 0, Izz_sen];
% Translation vector from sensor frame to body frame:
sensor_frame_origin = [sp_w/2 - sensor_w/2, 0, sp_h + busRec_h + sensor_h/2];
T_sen_b = sensor_frame_origin - CM;
% Sensor Inertia Tensor about CM:
I_sen_b = ParallelAxisTheorem(I_sen, T_sen_b, m_sen);

I      = I_sp_b + I_sen_b + I_bus_b;
I_lbin = I .* 3410; % to compare with CAD





