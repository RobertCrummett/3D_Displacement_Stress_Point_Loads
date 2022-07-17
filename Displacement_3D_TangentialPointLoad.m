function [u, v, w] = Displacement_3D_TangentialPointLoad(x, y, z, T, G, nu)
% Calculates displacements resulting from a point load at the origin
% 
% Jaeger, J.C., Cook, N.G.W. and Zimmerman, R.W. (2007) Fundamentals of Rock 
% Mechanics. 4th Edition, Chapman and Hall, London. pg. 410, 411
%
% Inputs:
% x - One dimensional vector of x values to evalute displacements over. [m]
% y - One dimensional vector of y values to evalute displacements over. [m]
% z - One dimensional vector of z values to evalute displacements over. [m]
% T - Point tangential force applied at the origin (0,0,0) in the +x direction. [N]
%
% Optional Inputs:
% G - Shear modulus. Default is 35e9 [Pa]
% nu - Poissons ratio. Default is 0.25 [unitless]
%
% Outputs:
% u - Resulting displacement in the x direction at every point in mesh. [m]
% v - Resulting displacement in the y direction at every point in mesh. [m]
% w - Resulting displacement in the z direction at every point in mesh. [m]

if ~exist('G', 'var')
    G = 35e9;
end

if ~exist('nu','var')
    nu = 0.25;
end

[X, Y, Z] = meshgrid(x, y, z); % 3D mesh grid [m] [m] [m] 

R = sqrt(X.^2 + Y.^2 + Z.^2); % distance from (0,0,0) [m]

u = T/(4*pi*G)*((X.^2 + R.^2)./(R.^3) + (1 - 2*nu)./(Z + R)...
    - ((1 - 2*nu)*X.^2)./(R.*(Z + R).^2)); % [m]
v = T/(4*pi*G)*((X.*Y)./(R.^3) - ((1 - 2*nu)*X.*Y)./(R.*(Z + R).^2)); % [m]
w = T/(4*pi*G)*((X.*Z)./(R.^3) + ((1 - 2*nu)*X)./(R.*(Z + R))); % [m]
end