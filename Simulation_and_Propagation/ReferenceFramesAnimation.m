%% Reference frames animation:
clear;
load("Data/OrbitPropData.mat")
load("Data/Body_PA_in_Inertial.mat")

b_x_i = b_x_i';
b_y_i = b_y_i';
b_z_i = b_z_i';
pa_x_i = pa_x_i';
pa_y_i = pa_y_i';
pa_z_i = pa_z_i';



%% Make Animation
frames    = struct('cdata', [], 'colormap', []);
figHandle = figure;
origin    = [0,0,0];
offset    = 0.1;

for i = 1:length(tout)
    view(3)

    %%%%%%% RTN FRAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quiver3(origin(1), origin(2), origin(3), R_unitvec(i,1), R_unitvec(i,2), R_unitvec(i,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    hold on;
    quiver3(origin(1), origin(2), origin(3), T_unitvec(i,1), T_unitvec(i,2), T_unitvec(i,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), N_unitvec(i,1), N_unitvec(i,2), N_unitvec(i,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    text(R_unitvec(i,1), R_unitvec(i,2)+offset, R_unitvec(i,3), 'R', 'FontSize', 14, 'color', 'red');
    text(T_unitvec(i,1), T_unitvec(i,2)+offset, T_unitvec(i,3), 'T', 'FontSize', 14, 'color', 'red');
    text(N_unitvec(i,1), N_unitvec(i,2), N_unitvec(i,3)+offset, 'N', 'FontSize', 14, 'color', 'red');

    %%%%%%% PRINCIPAL AXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quiver3(origin(1), origin(2), origin(3), pa_x_i(i,1), pa_x_i(i,2), pa_x_i(i,3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), pa_y_i(i,1), pa_y_i(i,2), pa_y_i(i,3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), pa_z_i(i,1), pa_z_i(i,2), pa_z_i(i,3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
    text(pa_x_i(i,1), pa_x_i(i,2)+offset, pa_x_i(i,3), 'X_p', 'FontSize', 14, 'color', 'blue');
    text(pa_y_i(i,1), pa_y_i(i,2)+offset, pa_y_i(i,3), 'Y_p', 'FontSize', 14, 'color', 'blue');
    text(pa_z_i(i,1), pa_z_i(i,2), pa_z_i(i,3)+offset, 'Z_p', 'FontSize', 14, 'color', 'blue');

    %%%%%% BODY AXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quiver3(origin(1), origin(2), origin(3), b_x_i(i,1), b_x_i(i,2), b_x_i(i,3), 'k', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), b_y_i(i,1), b_y_i(i,2), b_y_i(i,3), 'k', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), b_z_i(i,1), b_z_i(i,2), b_z_i(i,3), 'k', 'LineWidth', 2, 'MarkerSize', 5);
    text(b_x_i(i,1), b_x_i(i,2)+offset, b_x_i(i,3), 'X_p', 'FontSize', 14, 'color', 'black');
    text(b_y_i(i,1), b_y_i(i,2)+offset, b_y_i(i,3), 'Y_p', 'FontSize', 14, 'color', 'black');
    text(b_z_i(i,1), b_z_i(i,2), b_z_i(i,3)+offset, 'Z_p', 'FontSize', 14, 'color', 'black');

    plot3(0,0,0,'.','MarkerSize', 20)
    axis equal; 
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    title(['Reference Frames at t = ', num2str(tout(i)), ' s'])
    frames(i) = getframe(figHandle);
    hold off;

end



movie(figHandle, frames, 1, 10);

videoFile = VideoWriter("Figures_and_Plots/PS3/ReferenceFrames_animation.mp4",'MPEG-4');
videoFile.FrameRate = 10;
open(videoFile);
writeVideo(videoFile, frames);
close(videoFile);

%% Make Plot

%%%% Make a figure that shows at 2 different points in orbit:

figure()
% plot3(Xout(:,1),Xout(:,2),Xout(:,3),'Color','r','LineWidth',3)
% hold on;
% grid on;
% axis equal;
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% % title(['CYGNSS Orbit Propagation with No Perturbing Forces for ', num2str(numPeriods), ' orbits']);
% fontsize(14,'points')
% 
% %plot Earth-sized sphere
% rE = 5000; %6378;
% [xE,yE,zE] = ellipsoid(0,0,0,rE,rE,rE,20);
% surface(xE, yE , zE , 'FaceColor', 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% view (3);

tsample = [1, 100, 250];

for tsamp = 1:length(tsample)

    img = tsample(tsamp);
    view(3)

    %%%%%%% RTN FRAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quiver3(origin(1), origin(2), origin(3), R_unitvec(img,1), R_unitvec(img,2), R_unitvec(img,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    hold on;
    quiver3(origin(1), origin(2), origin(3), T_unitvec(img,1), T_unitvec(img,2), T_unitvec(img,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), N_unitvec(img,1), N_unitvec(img,2), N_unitvec(img,3), 'r', 'LineWidth', 2, 'MarkerSize', 5);
    text(R_unitvec(img,1), R_unitvec(img,2)+offset, R_unitvec(img,3), 'R', 'FontSize', 14, 'color', 'red');
    text(T_unitvec(img,1), T_unitvec(img,2)+offset, T_unitvec(img,3), 'T', 'FontSize', 14, 'color', 'red');
    text(N_unitvec(img,1), N_unitvec(img,2), N_unitvec(img,3)+offset, 'N', 'FontSize', 14, 'color', 'red');

    %%%%%%% PRINCIPAL AXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quiver3(origin(1), origin(2), origin(3), pa_x_i(img,1), pa_x_i(img,2), pa_x_i(img,3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), pa_y_i(img,1), pa_y_i(img,2), pa_y_i(img,3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), pa_z_i(img,1), pa_z_i(img,2), pa_z_i(img,3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
    text(pa_x_i(img,1), pa_x_i(img,2)+offset, pa_x_i(img,3), 'X_p', 'FontSize', 14, 'color', 'blue');
    text(pa_y_i(img,1), pa_y_i(img,2)+offset, pa_y_i(img,3), 'Y_p', 'FontSize', 14, 'color', 'blue');
    text(pa_z_i(img,1), pa_z_i(img,2), pa_z_i(img,3)+offset, 'Z_p', 'FontSize', 14, 'color', 'blue');

    %%%%%% BODY AXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quiver3(origin(1), origin(2), origin(3), b_x_i(img,1), b_x_i(img,2), b_x_i(img,3), 'k', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), b_y_i(img,1), b_y_i(img,2), b_y_i(img,3), 'k', 'LineWidth', 2, 'MarkerSize', 5);
    quiver3(origin(1), origin(2), origin(3), b_z_i(img,1), b_z_i(img,2), b_z_i(img,3), 'k', 'LineWidth', 2, 'MarkerSize', 5);
    text(b_x_i(img,1), b_x_i(img,2)+offset, b_x_i(img,3), 'X_p', 'FontSize', 14, 'color', 'black');
    text(b_y_i(img,1), b_y_i(img,2)+offset, b_y_i(img,3), 'Y_p', 'FontSize', 14, 'color', 'black');
    text(b_z_i(img,1), b_z_i(img,2), b_z_i(img,3)+offset, 'Z_p', 'FontSize', 14, 'color', 'black');

    plot3(0,0,0,'.','MarkerSize', 20)
    axis equal; 
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    title(['Reference Frames at t = ', num2str(tout(img)), ' s'])
    figname = ['ReferenceFrames_', num2str(tout(img)), 's.png'];
    saveas(gcf,['Figures_and_Plots/PS3', figname])
    hold off;
end