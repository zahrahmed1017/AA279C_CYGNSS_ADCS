%% Outer Surface Discretization:
function cygnss = CYGNSS_Geometry()

    %{ 
    cygnss_geom_struct: 
    - dimensions (wrt their respective frames)
        - sp
            - l, w, h
        - bus
            - top length (tl), bottom length (bl), full height (fh), short
            height (sh), width (w)
        - sensor
            - l, w, h
    - CM (wrt sp frame)
    - sp
        - origin (locations of the origins for the following frames/CM wrt solar
               panel frame)
        - coordinates
            - top, bottom_left, bottom_right
        - surface_area
            - top, bottom_left, bottom_right
        - norm_vec 
            - top, bottom_left, bottom_right
        - barycenter (wrt CM)
            - top, bottom_left, bottom_right
    - bus
        - origin (locations of the origins for the following frames/CM wrt solar
               panel frame)
        - coordinates
            - bottom, front, back, left1, left2, right1, right2
        - surface_area
            - bottom, front, back, left1, left2, right1, right2
        - norm_vec
            - bottom, front, back, left1, left2, right1, right2
        - barycenter (wrt CM)
            - bottom, front, back, left1, left2, right1, right2
    - sensor
        - origin (locations of the origins for the following frames/CM wrt solar
               panel frame)
        - coordinates
            - bottom
        - surface_area
            - bottom
        - norm_vec
            - bottom
        - barycenter (wrt CM)
            - bottom
                 
    %}

    %%%%%%%%%% Initialize Struct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cygnss = struct();

    %%%%%%%%%% CM Vector - defined from solar panel frame %%%%%%%%%%%%%%%%%%%%%
    load('InertiaData.mat','CM')
    % cygnss.CM = [6.759294029290275e-04,0,0.101246714232069];
    cygnss.CM = CM;
    
    %%%%%%%%%% Dimensions:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % length = along y axis
    % width  = along x axis
    % height = along z axis
    sp_l = 1.6;
    sp_w = 0.5;
    sp_h = 0.01; 
    % Satellite bus is a rectangular prism minus 2 triangular prisms
    bus_tl = 0.6; % top length
    bus_bl = 0.15;
    bus_w  = 0.5;
    bus_fh = 0.25; % full height
    bus_sh = 0.1;  % short side height
    % Sensor
    sen_l = 0.15;
    sen_w = 0.3;
    sen_h = 0.01;

    % Add to struct:
    cygnss.dim.sp.l = sp_l;
    cygnss.dim.sp.w = sp_l;
    cygnss.dim.sp.h = sp_h;
    cygnss.dim.bus.tl = bus_tl;
    cygnss.dim.bus.bl = bus_bl;
    cygnss.dim.bus.w  = bus_w;
    cygnss.dim.bus.fh = bus_fh;
    cygnss.dim.bus.sh = bus_sh;
    cygnss.dim.sen.l  = sen_l;
    cygnss.dim.sen.w  = sen_w;
    cygnss.dim.sen.h  = sen_h;

    %%%%%%%%%% Surface Coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Each surface should have 4 rows for the 4 vertices [x, y, z]
    % Edges should be formed from vertices 1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1
    % Coordinates are all given in a frame attached to the top of the solar panel
    
    % Solar Panel: (Coordinates given in the solar panel frame)
    cygnss.sp.coord.top = [ sp_w/2,  sp_l/2 ,   0;...
                           -sp_w/2,  sp_l/2,    0;...
                           -sp_w/2, -sp_l/2,    0;...
                            sp_w/2, -sp_l/2,    0];
    cygnss.sp.coord.bl = [  sp_w/2,  sp_l/2,    sp_h;...
                           -sp_w/2,  sp_l/2,    sp_h;...
                           -sp_w/2,  bus_tl/2,  sp_h;...
                            sp_w/2,  bus_tl/2,  sp_h];
    cygnss.sp.coord.br = [  sp_w/2, -sp_l/2,    sp_h;...
                           -sp_w/2, -sp_l/2,    sp_h;...
                           -sp_w/2, -bus_tl/2,  sp_h;...
                            sp_w/2, -bus_tl/2,  sp_h];
    cygnss.sp.origin   = [0,0,0];
    
    % Satellite Bus: Coordinates given in satellite bus frame
    cygnss.bus.coord.bottom = [ bus_w/2,  bus_bl/2, bus_fh/2;...
                               -bus_w/2,  bus_bl/2, bus_fh/2;...
                               -bus_w/2, -bus_bl/2, bus_fh/2;...
                                bus_w/2, -bus_bl/2, bus_fh/2];
    cygnss.bus.coord.front  = [ bus_w/2,  bus_tl/2, -bus_fh/2;...
                                bus_w/2,  bus_tl/2, -(bus_fh/2 - bus_sh);...
                                bus_w/2,  bus_bl/2, bus_fh/2;...
                                bus_w/2, -bus_bl/2, bus_fh/2;...
                                bus_w/2, -bus_tl/2, -(bus_fh/2 - bus_sh);...
                                bus_w/2, -bus_tl/2, -bus_fh/2];
    cygnss.bus.coord.back =  [-bus_w/2,  bus_tl/2, -bus_fh/2;...
                               -bus_w/2,  bus_tl/2, -(bus_fh/2 - bus_sh);...
                               -bus_w/2,  bus_bl/2, bus_fh/2;...
                               -bus_w/2, -bus_bl/2, bus_fh/2;...
                               -bus_w/2, -bus_tl/2, -(bus_fh/2 - bus_sh);...
                               -bus_w/2, -bus_tl/2, -bus_fh/2];
    cygnss.bus.coord.left1 =  [ bus_w/2,  bus_tl/2, -bus_fh/2;...
                               -bus_w/2,  bus_tl/2, -bus_fh/2;...
                               -bus_w/2,  bus_tl/2, -(bus_fh/2 - bus_sh);...
                                bus_w/2,  bus_tl/2, -(bus_fh/2 - bus_sh)];
    cygnss.bus.coord.left2  = [ bus_w/2,  bus_tl/2, -(bus_fh/2 - bus_sh);...
                               -bus_w/2,  bus_tl/2, -(bus_fh/2 - bus_sh);...
                               -bus_w/2,  bus_bl/2, bus_fh/2;...
                                bus_w/2,  bus_bl/2, bus_fh/2];
    cygnss.bus.coord.right1 = [ bus_w/2, -bus_tl/2, -bus_fh/2;...
                               -bus_w/2, -bus_tl/2, -bus_fh/2;...
                               -bus_w/2, -bus_tl/2, -(bus_fh/2 - bus_sh);...
                                bus_w/2, -bus_tl/2, -(bus_fh/2 - bus_sh)];
    cygnss.bus.coord.right2 = [ bus_w/2, -bus_tl/2, -(bus_fh/2 - bus_sh);...
                               -bus_w/2, -bus_tl/2, -(bus_fh/2 - bus_sh);...
                               -bus_w/2, -bus_bl/2, bus_fh/2;...
                                bus_w/2, -bus_bl/2, bus_fh/2];
    
    cygnss.bus.origin = [0, 0, sp_h + bus_fh/2]; % wrt sp frame
    
    % Sensor: Coordinate given in sensor frame
    cygnss.sensor.coord.bottom = [ sen_w/2,  sen_l/2, sen_h/2;...
                                  -sen_w/2,  sen_l/2, sen_h/2;...
                                  -sen_w/2, -sen_l/2, sen_h/2;...
                                   sen_w/2, -sen_l/2, sen_h/2];
    cygnss.sensor.origin = [sp_w/2 - sen_w/2, 0, sp_h + bus_fh + sen_h/2]; % wrt sp frame

    %%%%%%%%%%%%% Face Direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cygnss.sp.faceDir.top       = 'out';
    cygnss.sp.faceDir.bl        = 'in';
    cygnss.sp.faceDir.br        = 'in';

    cygnss.bus.faceDir.bottom   = 'out';
    cygnss.bus.faceDir.front    = 'out';
    cygnss.bus.faceDir.back     = 'out';
    cygnss.bus.faceDir.left1    = 'out';
    cygnss.bus.faceDir.left2    = 'out';
    cygnss.bus.faceDir.right1   = 'out';
    cygnss.bus.faceDir.right2   = 'out';
    
    cygnss.sensor.faceDir.bottom  = 'out';


    %%%%%%%%%%%%% Surface Normal, Barycenter, and Area %%%%%%%%%%%%%%%%%%%%

    cygnss = SurfaceDiscretization(cygnss, 'sp');
    cygnss = SurfaceDiscretization(cygnss, 'bus');
    cygnss = SurfaceDiscretization(cygnss, 'sensor');

    save('GeometryMassProperties/cygnss.mat', 'cygnss');

end