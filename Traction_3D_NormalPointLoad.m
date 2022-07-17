function [tau_xx, tau_yy, tau_zz, tau_xy, tau_yz, tau_xz] = ...
    Traction_3D_NormalPointLoad(x, y, z, N, nu)
% Calculates tractions resulting from a point load at the origin
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
% nu - Poissons ratio. Default is 0.25 [unitless]
%
% Outputs:
% tau_xx - Traction on the x face, normal to the x face. [N m^-2]
% tau_yy - Traction on the y face, normal to the y face. [N m^-2]
% tau_zz - Traction on the z face, normal to the z face. [N m^-2]
% tau_xy - Traction on the x face, in the y direction. Change sign for the
%                 traction on the y face, in the x direction. [N m^-2]
% tau_yz - Traction on the y face, in the z direction. Change sign for the
%                 traction on the z face, in the y direction. [N m^-2]
% tau_xz - Traction on the x face, in the z direction. Change sign for the
%                 traction on the z face, in the x direction. [N m^-2]

if ~exist('nu','var')
    nu = 0.25;
end

[X, Y, Z] = meshgrid(x, y, z); % 3D mesh grid [m] [m] [m] 

R = sqrt(X.^2 + Y.^2 + Z.^2); % distance from (0,0,0) [m]

tau_xx = N/(2*pi)*(((3*X.^2.*Z)./(R.^5)) + ((1 - 2*nu)*(Y.^2 + Z.^2))./(R.^3.*(Z + R)) - ...
    ((1 - 2*nu)*Z)./(R.^3) - ((1 - 2*nu)*X.^2)./(R.^2.*(Z + R).^2)); % [N m^-2]
tau_yy = N/(2*pi)*(((3*Y.^2.*Z)./(R.^5)) + ((1 - 2*nu)*(X.^2 + Z.^2))./(R.^3.*(Z + R)) - ...
    ((1 - 2*nu)*Z)./(R.^3) - ((1 - 2*nu)*Y.^2)./(R.^2.*(Z + R).^2)); % [N m^-2]
tau_zz = N/(2*pi)*((3*Z.^3)./(R.^5)); % [N m^-2]
tau_xy = N/(2*pi)*((3*X.*Y.*Z)./(R.^5) - ((1 - 2*nu)*X.*Y.*(Z + 2*R))./...
    (R.^3.*(Z + R).^2)); % [N m^-2]
tau_yz = N/(2*pi)*((3*Y.*Z.^2)./(R.^5)); % [N m^-2]
tau_xz = N/(2*pi)*((3*X.*Z.^2)./(R.^5)); % [N m^-2]
end