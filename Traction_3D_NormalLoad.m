function [tau_xx, tau_yy, tau_zz, tau_xy, tau_yz, tau_xz] =...
    Traction_3D_NormalLoad(x,y,z,traction,xmin,xmax,ymin,ymax, nu)
% Calculates tractions resulting from a normal load.
% Integrations are preformed with quad2d() function. Absolute error is set
% at 1e-10 and the 'Singular' keyword is set to true. 
% 
% Jaeger, J.C., Cook, N.G.W. and Zimmerman, R.W. (2007) Fundamentals of Rock 
% Mechanics. 4th Edition, Chapman and Hall, London. pg. 408 - 411
%
% Inputs:
% x - One dimensional vector of x values to evalute displacements over. [m]
% y - One dimensional vector of y values to evalute displacements over. [m]
% z - One dimensional vector of z values to evalute displacements over. [m]
%     Should not contain zeros in order for integrations to converge.
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

if nnz(z) ~= length(z)
    warning('If z contains 0 integration is not safe')
end

if ~exist('nu','var')
    nu = 0.25;
end

[X,Y,Z] = meshgrid(x,y,z); % [m] [m] [m]
r = @(x,y,z,zeta,eta) sqrt((x - zeta).^2 + (y - eta).^2 + z.^2); % [m]

fun_xx = @(x,y,z,zeta,eta) traction/(2*pi)*(((3*(x - zeta).^2.*z)./...
    (r(x,y,z,zeta,eta).^5)) + ((1 - 2*nu)*((y - eta).^2 + z.^2))./...
    (r(x,y,z,zeta,eta).^3.*(z + r(x,y,z,zeta,eta))) - ((1 - 2*nu)*z)./...
    (r(x,y,z,zeta,eta).^3) - ((1 - 2*nu)*(x - zeta).^2)./...
    (r(x,y,z,zeta,eta).^2.*(z + r(x,y,z,zeta,eta)).^2)); % [N m^-4]
fun_yy = @(x,y,z,zeta,eta) traction/(2*pi)*(((3*(y - eta).^2.*z)./...
    (r(x,y,z,zeta,eta).^5)) + ((1 - 2*nu)*((x - zeta).^2 + z.^2))./...
    (r(x,y,z,zeta,eta).^3.*(z + r(x,y,z,zeta,eta))) - ((1 - 2*nu)*z)./...
    (r(x,y,z,zeta,eta).^3) - ((1 - 2*nu)*(y - eta).^2)./...
    (r(x,y,z,zeta,eta).^2.*(z + r(x,y,z,zeta,eta)).^2)); % [N m^-4]
fun_zz = @(x,y,z,zeta,eta) traction/(2*pi)*((3*z.^3)./(r(x,y,z,zeta,eta).^5)); 
    % [N m^-4]
fun_xy = @(x,y,z,zeta,eta) traction/(2*pi)*((3*(x - zeta).*(y - eta).*z)./...
    (r(x,y,z,zeta,eta).^5) - ((1 - 2*nu)*(x - zeta).*(y - eta).*...
    (z + 2*r(x,y,z,zeta,eta)))./(r(x,y,z,zeta,eta).^3.*...
    (z + r(x,y,z,zeta,eta)).^2)); % [N m^-4]
fun_yz = @(x,y,z,zeta,eta) traction/(2*pi)*((3*(y - eta).*z.^2)./...
    (r(x,y,z,zeta,eta).^5)); % [N m^-4]
fun_xz = @(x,y,z,zeta,eta) traction/(2*pi)*((3*(x - zeta).*z.^2)./...
    (r(x,y,z,zeta,eta).^5)); % [N m^-4]

% Preallocation
tau_xx = zeros(length(x),length(y),length(z));
tau_yy = zeros(length(x),length(y),length(z));
tau_zz = zeros(length(x),length(y),length(z));
tau_xy = zeros(length(x),length(y),length(z));
tau_xz = zeros(length(x),length(y),length(z));
tau_yz = zeros(length(x),length(y),length(z));

tic
disp(0)
for i = 1:length(x)
    for j = 1:length(y)
        for k = 1:length(z)
            x_eval = X(i,j,k);
            y_eval = Y(i,j,k);
            z_eval = Z(i,j,k);

            fun_xx_int = @(zeta, eta) fun_xx(x_eval, y_eval, z_eval, zeta, eta);
            fun_yy_int = @(zeta, eta) fun_yy(x_eval, y_eval, z_eval, zeta, eta);
            fun_zz_int = @(zeta, eta) fun_zz(x_eval, y_eval, z_eval, zeta, eta);
            fun_xy_int = @(zeta, eta) fun_xy(x_eval, y_eval, z_eval, zeta, eta);
            fun_xz_int = @(zeta, eta) fun_xz(x_eval, y_eval, z_eval, zeta, eta);
            fun_yz_int = @(zeta, eta) fun_yz(x_eval, y_eval, z_eval, zeta, eta);

            tau_xx(i,j,k) = quad2d(fun_xx_int, xmin, xmax, ymin, ymax,...
                "AbsTol",1e-10,"RelTol",1e-10);
            tau_yy(i,j,k) = quad2d(fun_yy_int, xmin, xmax, ymin, ymax,...
                "AbsTol",1e-10,"RelTol",1e-10);
            tau_zz(i,j,k) = quad2d(fun_zz_int, xmin, xmax, ymin, ymax,...
                "AbsTol",1e-10,"RelTol",1e-10);
            tau_xy(i,j,k) = quad2d(fun_xy_int, xmin, xmax, ymin, ymax,...
                "AbsTol",1e-10,"RelTol",1e-10);
            tau_xz(i,j,k) = quad2d(fun_xz_int, xmin, xmax, ymin, ymax,...
                "AbsTol",1e-10,"RelTol",1e-10);
            tau_yz(i,j,k) = quad2d(fun_yz_int, xmin, xmax, ymin, ymax,...
                "AbsTol",1e-10,"RelTol",1e-10);
        end
        disp(strcat(num2str(round(i/length(x)*100 - 1/length(x)*100 +...
            j/(length(y)*length(x))*100, 2 ,'decimals')),"%"))
    end
end
disp('Done')
toc
disp('This is really freakin slow')
end