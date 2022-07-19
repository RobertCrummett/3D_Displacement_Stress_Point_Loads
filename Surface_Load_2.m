%% repeating code
clear;

sigma = 2.5e9; % Normal traction applied at origin of half-space [N]
G = 35e9; % Shear Modulus [Pa]
nu = 0.25; % Poissons Ratio

x = linspace(-10, 10, 20);
y = linspace(-10, 10, 20);
z = linspace(0.1, 10, 10);

xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;
%%
disp('start')
[tau_xx, tau_yy, tau_zz, tau_xy, tau_yz, tau_xz] =...
    Traction_3D_NormalLoad(x,y,z,sigma,xmin,xmax,ymin,ymax);

%% 

[traction_interp, X_interp, Y_interp, Z_interp] =...
    Interpolate_Traction_3D(x,y,z,tau_yz,1/5);

figure()
slice(X_interp, Y_interp, Z_interp, traction_interp, 0, 0, [])
shading flat
set(gca, 'ZDir', 'reverse')
cbar = colorbar();
title(cbar, "$N/m^2$", 'Interpreter','latex')
clim("auto")

xlabel('X axis','Interpreter','latex')
ylabel('Y axis','Interpreter','latex')
zlabel('Z axis','Interpreter','latex')
title(strcat("Stress Feild due to ", num2str(sigma, '%.1e'), " $N m^{-2}$ Traction"),...
    'Interpreter','latex')
% figure(1)
% subplot(1,3,1)
% surface(X(:,:,1),Y(:,:,1),u(:,:,1))
% colorbar
% colormap jet
% xlabel('x'); ylabel('y'); zlabel('z')
% %set(gca, 'ZDir', 'reverse')


% subplot(1,3,2)
% surface(X(:,:,1),Y(:,:,1),v(:,:,1))
% colorbar
% colormap jet
% xlabel('x'); ylabel('y'); zlabel('z')
% %set(gca, 'ZDir', 'reverse')
% 
% subplot(1,3,3)
% surface(X(:,:,1),Y(:,:,1), w(:,:,1))
% colorbar
% colormap jet
% xlabel('x'); ylabel('y'); zlabel('z')
% set(gca, 'ZDir', 'reverse')
% 
% figure(2)
% subplot(1,3,1)
% quiver(X(:,:,1),Y(:,:,1),u(:,:,1), zeros(size(u(:,:,1))))
% xlabel('x'); ylabel('y');
% 
% subplot(1,3,2)
% quiver(X(:,:,1),Y(:,:,1),zeros(size(v(:,:,1))), v(:,:,1))
% xlabel('x'); ylabel('y');
% 
% subplot(1,3,3)
% quiver(X(:,:,1),Y(:,:,1),u(:,:,1),v(:,:,1))
% xlabel('x'); ylabel('y');