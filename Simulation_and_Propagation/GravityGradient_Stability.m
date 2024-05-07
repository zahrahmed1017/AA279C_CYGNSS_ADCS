%% PS5 - Q1 - Gravity Gradient Stability

close all; clear;

load("InertiaData.mat")
load("OrbitPropData_equil.mat")

%% Part a - Calculate k coefficients; plot stable/unstable regions

% create the stability plot - shade UNSTABLE regions
figure
hold on

alpha_patch = 0.2;

% K_R >= K_T 
p1_x = [-1, 1, -1, -1];
p1_y = [-1, 1, 1, -1];
p1 = patch(p1_x, p1_y, 'y');
p1.FaceVertexAlphaData = alpha_patch;
p1.FaceAlpha = 'flat';

% K_R * K_T <= 0
p2_x = [-1, -1, 0, 0, 1, 1, -1];
p2_y = [0, 1, 1, -1, -1, 0, 0];
p2 = patch(p2_x, p2_y, 'b');
p2.FaceVertexAlphaData = alpha_patch;
p2.FaceAlpha = 'flat';

% 1 + 3*K_T + K_R*K_T <= 4*sqrt(K_R*K_T)
p3_x = linspace(-1, -0.0505, 100);
p3_y = (-3* p3_x.^2 + 4*sqrt(3)*sqrt(p3_x.^2 - p3_x.^3) + 7*p3_x) ./ ...
    (p3_x.^2);
p3_x = [p3_x, -1, -1];
p3_y = [p3_y, -1, 0];
p3 = patch(p3_x, p3_y, 'r');
p3.FaceVertexAlphaData = alpha_patch;
p3.FaceAlpha = 'flat';

axis equal
ylim([-1, 1])
xlim([-1, 1])
xlabel("$K_T$", 'Interpreter', 'latex')
ylabel("$K_R$", 'Interpreter', 'latex')



% First K coefficients: z aligned w/ N

Ix = I_p(1,1);
Iy = I_p(2,2);
Iz = I_p(3,3);

K_N1 = (Iy - Ix) / Iz;
K_R1 = (Iz - Iy) / Ix;
K_T1 = (Iz - Ix) / Iy;

plot(K_T1, K_R1, 'kx')

% 
% K_N1 = (Iy - Ix) / Iz;
% K_R1 = (Iz - Iy) / Ix;
% K_T1 = (Iz - Ix) / Iy;


legend("$K_R \geq K_T$", "$K_R K_T \leq 0$", "$1 + 3 K_T + K_R K_T \leq 4\sqrt{K_R*K_T}$", ...
     "$Z_p$ aligned with N", 'Interpreter', 'latex', 'location', 'eastoutside')







