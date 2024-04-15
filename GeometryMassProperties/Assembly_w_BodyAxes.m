%% Plot of Assembly with Body Frame

% For each component, for each face, plot the surface using "Patch".
load("GeometryMassProperties/cygnss.mat")
components = {'sp', 'bus', 'sensor'};

figure()
grid on;
hold on;
xlabel('X')
ylabel('Y')
zlabel('Z')
fontsize(16,'points')
% set(gca, 'YDir', 'reverse');
set(gca, 'ZDir', 'reverse');
set(gca, 'XDir', 'reverse');

for i = 1 : length(components)
    comp      = components{i};
    faces     = fieldnames(cygnss.(comp).coord);
    T_comp_CM = cygnss.(comp).origin - cygnss.CM;

    for j = 1 : length(faces)
        faceName = faces{j};
        vertices = cygnss.(comp).coord.(faceName);

        x = vertices(:,1) + T_comp_CM(1);
        y = vertices(:,2) + T_comp_CM(2);
        z = vertices(:,3) + T_comp_CM(3);

        fill3(x, y, z, 'blue', 'FaceAlpha', 0.1);
    end
    
end


% Body Frame Axes
origin = [0, 0, 0];
x_axis = [0.2, 0, 0];
y_axis = [0, 0.2, 0];
z_axis = [0, 0, 0.2];

% Plot arrows for the x, y, and z axes
quiver3(origin(1), origin(2), origin(3), x_axis(1), x_axis(2), x_axis(3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
quiver3(origin(1), origin(2), origin(3), y_axis(1), y_axis(2), y_axis(3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
quiver3(origin(1), origin(2), origin(3), z_axis(1), z_axis(2), z_axis(3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
text(0.3, 0, 0, 'X_b', 'FontSize', 14, 'color', 'red');
text(0, 0.3, 0, 'Y_b', 'FontSize', 14, 'color', 'red');
text(0, 0, 0.3, 'Z_b', 'FontSize', 14, 'color', 'red');
plot3(0,0,0,'.','MarkerSize', 20)
axis equal; 


load("GeometryMassProperties/InertiaData.mat","R_b_p");
scalefactor = 0.2;
x_axis_p = R_b_p(:,1) * scalefactor;
y_axis_p = R_b_p(:,2) * scalefactor;
z_axis_p = R_b_p(:,3) * scalefactor;
quiver3(origin(1), origin(2), origin(3), x_axis_p(1), x_axis_p(2), x_axis_p(3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
quiver3(origin(1), origin(2), origin(3), y_axis_p(1), y_axis_p(2), y_axis_p(3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
quiver3(origin(1), origin(2), origin(3), z_axis_p(1), z_axis_p(2), z_axis_p(3), 'b', 'LineWidth', 2, 'MarkerSize', 5);

text(x_axis_p(1)+0.1, x_axis_p(2), x_axis_p(3), 'X_p', 'FontSize', 14, 'color', 'blue');
text(y_axis_p(1), y_axis_p(2)+0.1, y_axis_p(3), 'Y_p', 'FontSize', 14, 'color', 'blue');
text(z_axis_p(1), z_axis_p(2), z_axis_p(3)+0.1, 'Z_p', 'FontSize', 14, 'color', 'blue');
hold off; 
view(3)
