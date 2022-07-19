function [u,v,w] = Displacement_3D_NormalLoad(x,y,z,traction,...
    xmin,xmax,ymin,ymax, G, nu)
% Calculates displacements resulting from a normal load
% 
% Jaeger, J.C., Cook, N.G.W. and Zimmerman, R.W. (2007) Fundamentals of Rock 
% Mechanics. 4th Edition, Chapman and Hall, London. pg. 408 - 411
%
% Inputs:
% x - One dimensional vector of x values to evalute displacements over. [m]
% y - One dimensional vector of y values to evalute displacements over. [m]
% z - One dimensional vector of z values to evalute displacements over. [m]
% traction - The normal stress exerted by surface on the ground. [N m^-2]
%            This can be uniform. (input scalar)
%            This can be a function of x and y location. (input function)
%            Positive is loading down. (+z direction, where z is into the
%            ground)
%
% Surface Input:
% xmin - x minimum boundary [m]
% xmax - x maximum boundary [m]
% ymin - y minimum boundary [m]
% ymax - y maximum boundary [m]
% In order to create surfaces more dynamic than squares, y can be 
%   parameterized in terms of x.
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
    G = 35e9; % [Pa]
end

if ~exist('nu','var')
    nu = 0.25;
end

[X,Y,Z] = meshgrid(x,y,z); % [m] [m] [m]
r = @(x,y,z,zeta,eta) sqrt((x - zeta).^2 + (y - eta).^2 + z.^2); % [m]

fun_u = @(x,y,z,zeta,eta) traction/(4*pi*G)*(((1 - 2*nu)*(x-zeta))./...
    (r(x,y,z,zeta,eta).*(z + r(x,y,z,zeta,eta))) - ...
    ((x-zeta).*z)./(r(x,y,z,zeta,eta).^3)); % [m]
fun_v = @(x,y,z,zeta,eta) traction/(4*pi*G)*(((1 - 2*nu)*(y-eta))./...
    (r(x,y,z,zeta,eta).*(z + r(x,y,z,zeta,eta))) - ...
    ((y-eta).*z)./(r(x,y,z,zeta,eta).^3)); % [m]
fun_w = @(x,y,z,zeta,eta) - traction/(4*pi*G)*(((2*(1 - nu))./...
    r(x,y,z,zeta,eta)) + (z.^2)./(r(x,y,z,zeta,eta).^3)); % [m]

% Preallocation
u = zeros(length(x),length(y),length(z));
v = zeros(length(x),length(y),length(z));
w = zeros(length(x),length(y),length(z));

tic
for i = 1:length(x)
    percenti = i/length(x)*100 - 5;
    for j = 1:length(y)
        percentdone = round(percenti + j/length(y)*5, 1, 'decimals');
        for k = 1:length(z)
            x_eval = X(i,j,k);
            y_eval = Y(i,j,k);
            z_eval = Z(i,j,k);

            fun_u_int = @(zeta, eta) fun_u(x_eval, y_eval, z_eval, zeta, eta);
            fun_v_int = @(zeta, eta) fun_v(x_eval, y_eval, z_eval, zeta, eta);
            fun_w_int = @(zeta, eta) fun_w(x_eval, y_eval, z_eval, zeta, eta);
            u(i,j,k) = integral2(fun_u_int, xmin, xmax, ymin, ymax);
            v(i,j,k) = integral2(fun_v_int, xmin, xmax, ymin, ymax);
            w(i,j,k) = integral2(fun_w_int, xmin, xmax, ymin, ymax);
        end
        disp(strcat(num2str(percentdone),"%"))
    end
end
disp('Done')
toc

% Sign flips
u = -u;
v = -v;
w = -w;

end