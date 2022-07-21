function [traction, normal_traction, shear_traction, rake] = ...
    Traction_On_Plane(n,sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz)

% normal unit vector
n = n/norm(n);
stress = zeros(3,3,numel(sigma_xx));

stress(1,1,:) = reshape(sigma_xx,1,[]);
stress(2,2,:) = reshape(sigma_yy,1,[]);
stress(3,3,:) = reshape(sigma_zz,1,[]);
stress(1,2,:) = reshape(sigma_xy,1,[]);
stress(2,1,:) = reshape(sigma_xy,1,[]);
stress(1,3,:) = reshape(sigma_xz,1,[]);
stress(3,1,:) = reshape(sigma_xz,1,[]);
stress(2,3,:) = reshape(sigma_yz,1,[]);
stress(3,2,:) = reshape(sigma_yz,1,[]);

traction = pagemtimes(n,stress);
normal_traction = pagemtimes(traction,n');
nmat = repmat(n,1,1,numel(sigma_xx));
shear_traction = sqrt(sum(traction.*traction,[1,2]) - (normal_traction).^2);

Shear = traction - pagemtimes(normal_traction,nmat);
Shear = Shear./sqrt(Shear(1,1,:).^2 + Shear(1,2,:).^2 + Shear(1,3,:).^2);
z_axis = [0,0,1];

if n(1,3) >= 0
    reference = cross(z_axis,n);
else
    reference = -cross(z_axis,n);
end
reference = reference/norm(reference);

rake = zeros(1,1,numel(sigma_xx));
for i = 1:length(Shear)
    if Shear(1,3,i) >= 0
        rake(1,1,i) = acos(dot(reference,Shear(:,:,i)));
    else
        rake(1,1,i) = -acos(dot(reference,Shear(:,:,i)));
    end
end

normal_traction = reshape(normal_traction,size(sigma_xx));
shear_traction = reshape(shear_traction,size(sigma_xx));
rake = reshape(rake,size(sigma_xx));

end