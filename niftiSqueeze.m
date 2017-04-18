% nifitSqueeze.m
%
% utility to fix SUMA surf2vol extraneous dimensions (also checks the
% fields of vistasoft nii structs, fixes them as necessary)
%
% optionally, fixes time dimension (need TR as argument for that); time
% units always sec
%
% Tommy Sprague, 4/4/2017



function nii = niftiSqueeze(nii, TR)

if nargin < 2
    TR = [];
end

singleton_dims = size(nii.data)==1;

% note - this shoudl be equal to below
% singleton_dims = find(nii.dim==1);

nii.data = squeeze(nii.data);
nii.dim = size(nii.data);
nii.ndim = length(nii.dim);

nii.pixdim = nii.pixdim(~singleton_dims);

if ~isempty(TR)
    nii.pixdim(4) = TR;
    nii.time_units = 'sec';
end

% TODO: make sure we have over 3 dims, etc.

return