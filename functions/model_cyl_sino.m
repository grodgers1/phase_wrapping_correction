function [m_sino] = model_cyl_sino(x0,r,npix,nproj)
%MODEL_CYL_SINO Generates a sinogram of projected thickness of a circle
%with specified radius and center
%   x0: center position at each angle (length(x0) should be nproj)
%   r: radius of cylinder
%   npix: number of detector pixels
%   nproj: number of projections

%   m_sino: will be a sinogram of the projected thickness of the cylinder

uvec = 0:(npix-1);
uvec = uvec - round(size(uvec,2)/2) + 0.5;
t = zeros(size(uvec)); % will contain projected thickness of uniform circle 
m_sino = zeros(npix,nproj);
for j = 1:nproj
    for i = 1:size(uvec,2)
        if abs(uvec(i)+x0(j)) > r % rays missing model circle
            t(i) = 0;
        else % rays going through model circle
            t(i) = 2 * sqrt(r^2 -(uvec(i)+x0(j))^2);
        end
    end
    m_sino(:,j) = t;
end

end

