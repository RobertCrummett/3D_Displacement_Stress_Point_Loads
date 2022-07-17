%% Vertical displacement due to circle load
% uniform normal traction over a circular load
clear; clf; close all;
G = 35e9; % Shear modulus [Pa]
nu = 0.25; % Poissons constant [unitless]
sig = 2.5e9; % Uniform normal traction over load area [N m^-2]

a = 100; % Radius of load [m]
b = 0:1:a;

w = -2*a*(1- nu)*sig/(pi*G)*ellipticE((b/a).^2);

%% Plot Cross Section of vertical displacement vs distance from origin
% b ranges from 0 -> a, where a is the radius of the circular load.
figure(1);
plot(b, w, linewidth = 2, color = 'r')

title('Vertical Displacement vs Distance from Origin')
xlabel('b, Distance from Origin [m]')
ylabel('w, Vertical Displacement Under Load [m]')
ylim([(min(w) - 1) 0])

% Check figrure(1) with these
center_displacement = -a*(1 - nu)*sig/G;
edge_displacement = -2*a*(1 - nu)*sig/(G * pi);

%% Cylindrical Coords: vertical displacement vs distance from origin
theta = 0:1/32*pi:2*pi;
[B, Theta] = meshgrid(b, theta);

W =  -2*a*(1- nu)*sig/(pi*G)*ellipticE((B/a).^2);

figure(2);
surf(B.*cos(Theta), B.*sin(Theta), W)
colormap('autumn')
colorbar

xlabel('X axis [m]')
ylabel('Y axis [m]')
zlabel('Vertical Displacement [m]')
title('Deformation under circular load')

% The region plotted is directly under the load - there
% are solutions for the region to the sides of the load
% in Davis and Selvadurai 1996, p.120
% I do not currently have access legally or illegally to
% this publication, so the work ends here.