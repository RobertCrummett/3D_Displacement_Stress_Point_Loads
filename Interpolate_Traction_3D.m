function [traction_interp, X_interp, Y_interp, Z_interp] = Interpolate_Traction_3D(x, y, z, ...
    traction_component, interpolation_factor, interpolation_method)
% Interpolates output of Traction_3D_PointLoad() for smoother plots
% 
% Inputs:
% x - One dimensional vector of x values to evalute displacements over. [m]
% y - One dimensional vector of y values to evalute displacements over. [m]
% z - One dimensional vector of z values to evalute displacements over. [m]
% traction_component - Specific traction component to display. Output from
%                                          Traction_3D_PointLoad() function. [N m^-2]
%
% Optional Inputs:
% interpolation_factor - Factor to multiply the step size of the original
%                                          grid by. Should be < 1 to increase smoothness of output.
%                                          Scalar input for uniform
%                                          interpolation or (1,3) vector to
%                                          independantly scale the x, y, z
%                                          axis respectively                             
%                                          Default value is 1/10
% interpolation_method - Method to interpolate between grid points by. Input as a string. See
%                                             MATLAB documentation on interp3() function for options.
%                                             Default method is "spline"
%
% Outputs:
% traction_interp - Interpolated traction feild. [N m^-2]
% X_interp - Mesh of interpolated X. [m]
% Y_interp - Mesh of interpolated Y. [m]
% Z_interp - Mesh of interpolated Z. [m]

if ~exist("interpolation_factor", "var")
    interpolation_factor = 1/10;
end

if ~exist("interpolation_method", "var")
    interpolation_method = "spline";
    disp("Interpolation method set to 'spline' as default")
end

if length(interpolation_factor) == 1
    interpolation_factor = interpolation_factor*ones(1,3);
end

x_interp = x(1) : (x(end) - x(1) + 1)/length(x)*interpolation_factor(1) : x(end); % [m]
y_interp = y(1): (y(end) - y(1) + 1)/length(y)*interpolation_factor(2) : y(end); % [m]
z_interp = z(1): (z(end) - z(1) + 1)/length(z)*interpolation_factor(3) : z(end); % [m]

[X, Y, Z] = meshgrid(x, y, z); % 3D mesh grid [m] [m] [m]
[X_interp, Y_interp, Z_interp] = meshgrid(x_interp,y_interp,z_interp); % 3D mesh grid [m] [m] [m]

traction_interp = interp3(X, Y, Z, traction_component, X_interp, Y_interp, Z_interp, interpolation_method);
% [N m^-2]
end