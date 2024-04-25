%% CYGNSS Orbital Properties:

% Define Orbital Properties
a   = 6903;         % km
e   = 0.00162;
i   = deg2rad(35);  % rad     
w   = deg2rad(60);  % rad
O   = deg2rad(120); % rad
v   = deg2rad(0);   % rad        
muE = 398600;       % [km^3/s^2]

T      = 2*pi*sqrt(a^3/muE); % period in seconds
T_days = T/(24 * 60 * 60);

oe = [a; e; i; O; w; v];

% Convert to ECI position and velocity for orbit propagation 
state = OE2ECI(oe, muE);

% Propagate Orbit for multiple orbits
numPeriods  = 1;
% tspan       = 0 : 60 : T * numPeriods; % simulate once an hour?
tspan       = 0:10:2500;
initial     = state;
options     = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[tout,Xout] = ode113(@(t,State) PropagateOrbit(State, muE),...
                     tspan, initial, options);

figure()
plot3(Xout(:,1),Xout(:,2),Xout(:,3),'Color','r','LineWidth',3)
hold on;
grid on;
axis equal;
xlabel('X')
ylabel('Y')
zlabel('Z')
% title(['CYGNSS Orbit Propagation with No Perturbing Forces for ', num2str(numPeriods), ' orbits']);
fontsize(14,'points')

%plot Earth-sized sphere
rE = 6378;
[xE,yE,zE] = ellipsoid(0,0,0,rE,rE,rE,20);
surface(xE, yE , zE , 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha', 0.25);
view (3);

% Calculating the R-T-N unit vectors:
R_unitvec = zeros(length(tout),3);
T_unitvec = zeros(length(tout),3);
N_unitvec = zeros(length(tout),3);

for i = 1:length(tout)
    rx = Xout(i,1);
    ry = Xout(i,2);
    rz = Xout(i,3);
    vx = Xout(i,4);
    vy = Xout(i,5);
    vz = Xout(i,6);
    r_norm = sqrt(rx^2 + ry^2 + rz^2);
    v_norm = sqrt(vx^2 + vy^2 + vz^2);

    R_i = [rx/r_norm, ry/r_norm, rz/r_norm];
    T_i = [vx/v_norm, vy/v_norm, vz/v_norm];
    N_i = cross(R_i,T_i);

    R_unitvec(i,:) = R_i;
    T_unitvec(i,:) = T_i;
    N_unitvec(i,:) = N_i;
    
end

frames    = struct('cdata', [], 'colormap', []);
figHandle = figure;
origin    = [0,0,0];
offset    = 0.1;

for i = 1:length(tout)
    view(3)
    quiver3(origin(1), origin(2), origin(3), R_unitvec(i,1), R_unitvec(i,2), R_unitvec(i,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    hold on;
    quiver3(origin(1), origin(2), origin(3), T_unitvec(i,1), T_unitvec(i,2), T_unitvec(i,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), N_unitvec(i,1), N_unitvec(i,2), N_unitvec(i,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    text(R_unitvec(i,1), R_unitvec(i,2)+offset, R_unitvec(i,3), 'R', 'FontSize', 14, 'color', 'red');
    text(T_unitvec(i,1), T_unitvec(i,2)+offset, T_unitvec(i,3), 'T', 'FontSize', 14, 'color', 'red');
    text(N_unitvec(i,1), N_unitvec(i,2), N_unitvec(i,3)+offset, 'N', 'FontSize', 14, 'color', 'red');
    plot3(0,0,0,'.','MarkerSize', 20)
    % axis equal; 
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    frames(i) = getframe(figHandle);
    hold off;

end



movie(figHandle, frames, 1, 10);

videoFile = VideoWriter("Figures_and_Plots/RTN_animation.mp4",'MPEG-4');
videoFile.FrameRate = 10;
open(videoFile);
writeVideo(videoFile, frames);
close(videoFile);

save('Data/OrbitPropData.mat',"R_unitvec", "T_unitvec", "N_unitvec", "tout", "Xout");