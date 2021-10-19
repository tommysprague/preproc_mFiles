% make_grayplot.m
%
% inputs: filenames, not structures. uses vista niftiRead rather than
% matlab for now
%
% will use voxels ~= 0 for mask
%
% if out_img is empty, just plot. if out_img is defined, plot hidden and
% save
%
% TCS 10/19/2021

function make_grayplot(data_nii_file,mask_nii_file,out_img)

data_nii = niftiRead(data_nii_file);

% if no mask given, find std dev ~= 0 voxels to analyze
if nargin < 2 || isempty(mask_nii_file)
    mask_nii = data_nii;
    mask_nii.data = mask_nii.data(:,:,:,1);
    mask_nii.data = std(double(data_nii.data),[],4) ~= 0;
else
    mask_nii = niftiRead(mask_nii_file);
end

% if not single/double, cast to double
if ~ismember(class(data_nii.data),{'single','double'})
    data_nii.data = double(data_nii.data);
end

data_mat = niftiExtract(data_nii,mask_nii);


% normalize, detrend, etc
data_mat_norm = detrend(data_mat,2); % linear (quad?) detrend
data_mat_norm = data_mat_norm./std(data_mat_norm,[],1);
%data_mat_norm = data_mat;

if nargin < 3 || isempty(out_img)

    % figure out scale...
    clims = prctile(data_mat_norm(:),[5 95]);

    figure;
    imagesc((0:(data_nii.dim(4)))*data_nii.pixdim(4),1:size(data_mat_norm,2),data_mat_norm.');
    %set(gca,'CLim',clims,'TickDir','out');
    imax = gca;
    imax.CLim = clims; imax.TickDir = 'out'; imax.Box = 'off';
    colormap gray;

    xlabel('Time (s)'); ylabel('Voxel');
    th = title(data_nii_file,'Interpreter','none');

    

end

return