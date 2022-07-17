%% Displacements and Stresses due to surface load in 3D
clear; clf; close all;

% Parameters
N = 2.5e9; % Normal force applied at origin of half-space [N]
G = 35e9; % Shear Modulus [Pa]
nu = 0.25; % Poissons Ratio

% Plotting Bounds
step = 1;
x_limits = [-10, 10];
y_limits = [-10, 10];
z_limit = 10;

% Slices to Display
x_slice = 0;
y_slice = [-1 1];
z_slice = [];
colormap('jet')

%% Constructing Evaluation Grid
x = x_limits(1): step: x_limits(2); % [m]
y = y_limits(1): step: y_limits(2); % [m]
z = 0: step: z_limit; % depth >= 0, [m]

%% 3D Displacements
% timed
tic
[u, v, w] = Displacement_3D_NormalPointLoad(x, y, z, N, G, nu);
toc

%% 3D Stresses
% timed
tic
[tau_xx, ~, ~, ~, ~, ~] = Traction_3D_NormalPointLoad(x, y, z, N, nu);
toc

% Traction component to display
traction_component = tau_xx;

%% Interpolation of Traction Grid
interp_scale = 1/15;
[stress_interp, X_interp, Y_interp, Z_interp] = Interpolate_Traction_3D(x, y, z, traction_component, interp_scale);

%% Plotting Traction Component Feild and Displacement Feild
% Traction component Feild
figure(1)
slice(X_interp, Y_interp, Z_interp, stress_interp, x_slice, y_slice, z_slice)
shading flat
set(gca, 'ZDir', 'reverse')
cbar = colorbar();
title(cbar, "$N/m^2$", 'Interpreter','latex')
clim("auto")

xlabel('X axis','Interpreter','latex')
ylabel('Y axis','Interpreter','latex')
zlabel('Z axis','Interpreter','latex')
title(strcat("Stress Feild due to ", num2str(N, '%.1e'), " $N$ Point Force at Origin"),...
    'Interpreter','latex')

% Displacement Feild
figure(2)
[X, Y, Z] = meshgrid(x, y, z);
quiver3(X, Y, Z, u, v, -w, 2, 'Color', '#7E2F8E', 'LineWidth', 0.7, 'MaxHeadSize', 0.35)
set(gca, 'ZDir', 'reverse')

xlabel('X axis','Interpreter','latex')
ylabel('Y axis','Interpreter','latex')
zlabel('Z axis','Interpreter','latex')
title(strcat("Displacements Due to ", num2str(N, '%.1e'), " $N$ Point Force at Origin"), 'Interpreter', 'latex')
legend('Cannot Figure Out How to Scale', 'Location','southoutside')

%% What about a Tangenetial Force applied at the origin?
T = 2.5e9; % Tangential force applied in +x direction [N]

[tau_xx, tau_yy, tau_zz, tau_xy, tau_yz, tau_xz] = ...
    Traction_3D_TangentialPointLoad(x,y,z,T);
[u,v,w] = Displacement_3D_TangentialPointLoad(x,y,z,T);

traction_component = tau_xx;

[stress_interp, X_interp, Y_interp, Z_interp] = Interpolate_Traction_3D(x, y, z, traction_component, interp_scale);

%% Plotting Traction Component Feild and Displacement Feild
% Traction component Feild
figure(3)
slice(X_interp, Y_interp, Z_interp, stress_interp, x_slice, y_slice, z_slice)
shading flat
set(gca, 'ZDir', 'reverse')
cbar = colorbar();
title(cbar, "$N/m^2$", 'Interpreter','latex')
clim("auto")

xlabel('X axis','Interpreter','latex')
ylabel('Y axis','Interpreter','latex')
zlabel('Z axis','Interpreter','latex')
title(strcat("Stress Feild due to ", num2str(N, '%.1e'), " $N$ Tangential Force at Origin"),...
    'Interpreter','latex')

% Displacement Feild
figure(4)
quiver3(X, Y, Z, u, v, -w, 2, 'Color', '#000000', 'LineWidth', 0.7, 'MaxHeadSize', 0.35)
set(gca, 'ZDir', 'reverse')

xlabel('X axis','Interpreter','latex')
ylabel('Y axis','Interpreter','latex')
zlabel('Z axis','Interpreter','latex')
title(strcat("Displacements Due to ", num2str(N, '%.1e'), " $N$ Point Force at Origin"), 'Interpreter', 'latex')
legend('Cannot Figure Out How to Scale', 'Location','southoutside')
