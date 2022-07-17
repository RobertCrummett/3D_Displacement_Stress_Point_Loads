function [u, v, w] = Displacement_3D_NormalPointLoad(x, y, z, N, G, nu)
% Calculates displacements resulting from a point load at the origin
% 
% Jaeger, J.C., Cook, N.G.W. and Zimmerman, R.W. (2007) Fundamentals of Rock 
% Mechanics. 4th Edition, Chapman and Hall, London. pg. 408, 409
%
% Inputs:
% x - One dimensional vector of x values to evalute displacements over. [m]
% y - One dimensional vector of y values to evalute displacements over. [m]
% z - One dimensional vector of z values to evalute displacements over. [m]
% N - Point normal force applied at the origin (0,0,0). [N]
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

u = N/(4*pi*G)*(((1 - 2*nu)*X)./(R.*(Z + R)) - (X.*Z)./(R.^3)); % [m]
v = N/(4*pi*G)*(((1 - 2*nu)*Y)./(R.*(Z + R)) - (Y.*Z)./(R.^3)); % [m]
w = -N/(4*pi*G)*(((2*(1 - nu))./R) + (Z.^2)./(R.^3)); % [m]
end