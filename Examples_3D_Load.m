%% Example File Using Traction_3D_NormalLoad()

% User Inputs
glacial_height = 1000;
G = 35e9; % Shear Modulus [Pa]
nu = 0.25; % Poissons Ratio

% Half-Space Boundaries
x = linspace(100, 900, 25);
y = linspace(100, 1900, 55);
z = linspace(10, 1000, 25);

% Load Boundaries
xmin = 250;
xmax = 750;
ymin = 250;
ymax = 1750;

g = 9.8;
density_ice = 917;
sigma = glacial_height*density_ice*g; % Normal traction [N m^-2]

%% Employing Function
disp('start')
[tau_xx, tau_yy, tau_zz, tau_xy, tau_yz, tau_xz] =...
    Traction_3D_NormalLoad(x,y,z,sigma,xmin,xmax,ymin,ymax);

%% Normal Traction
n = [1,1,1]; % z axis
[Traction, Normal_Traction, Shear_Traction, Shear_Rake] = ...
    Traction_On_Plane(n,tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz);
disp('ok')

%% Interpolating Stress Feilds and Plotting
% Interpolation
[traction_interp, Y_interp, X_interp, Z_interp] =...
    Interpolate_Traction_3D(y,x,z,Normal_Traction,1/5);

%% Slicing Feild
figure()
a = n(1); b = n(2); c = n(3);
[xsurf, ysurf] = meshgrid(min(x):1:max(x),min(y):1:max(y));
zsurf3 = - a/c*(xsurf) - b/c*(ysurf) + 1200;
zsurf2 = - a/c*(xsurf) - b/c*(ysurf) + 1400;
zsurf1 = - a/c*(xsurf) - b/c*(ysurf) + 1000;

slice(X_interp, Y_interp, Z_interp, traction_interp, 500,...
    [200 400 1000], [], 'cubic')


shading flat
set(gca, 'ZDir', 'reverse')
cbar = colorbar();
colormap jet
title(cbar, "$N m^{-2}$", 'Interpreter','latex')
clim('auto')

axis equal
xlabel('X axis [m]','Interpreter','latex')
ylabel('Y axis [m]','Interpreter','latex')
zlabel('Z axis [m]','Interpreter','latex')
title(strcat("Stress Feild due to ", num2str(glacial_height), " m of Ice"),...
    'Interpreter','latex')
