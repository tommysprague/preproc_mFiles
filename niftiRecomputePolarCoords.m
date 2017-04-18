% niftiRecomputePolarAngle.m
%
% after surfsmooth, need to recompute polar angle - atan2(y0/x0)



function mynii = niftiRecomputePolarCoords(mynii,angle_vol,ecc_vol,x_vol, y_vol)

% TODO: introspect nii to find polar angle, etc briks
% ("default", as of april 2017, is for angle to be 2nd volume, and x/y to
% be end-2 and end-1)

if nargin < 2
    angle_vol = 2;
end

if nargin < 3
    ecc_vol = 3;
end

if nargin < 4
    x_vol = 6;
end

if nargin < 5
    y_vol = 7;
end

% make sure 5th dim is clean (assuming that's not useful for anything?)
mynii = niftiSqueeze(mynii);   % should this be separate?

mynii.data(:,:,:,angle_vol) = mod(atan2(mynii.data(:,:,:,y_vol),mynii.data(:,:,:,x_vol)),2*pi);
mynii.data(:,:,:,ecc_vol) = sqrt( mynii.data(:,:,:,x_vol).^2 + mynii.data(:,:,:,y_vol).^2 );


return