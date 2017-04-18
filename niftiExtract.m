% niftiExtract.m
%
% function for extracting timeseries/values from ROIs
%
% takes in a functional/task/RF nii struct, an ROI mask,
% and returns a matlab array (or cell array of matrices) containing voxel
% activity on each TR (rows) from each voxel (colums); or more geerally,
% from each brik of the nii.gz (rows). Later on, can re-format these as a
% struct of RF properties, etc. 
%
% (NOTE: only supports single ROIs, single nii's now; for speed reasons,
% loading these only once will be desirable)
%
% TCS, 4/18/2017
%



function data_mat = niftiExtract(nii_task,nii_roi)

% TODO: if nii_roi has 4 dimensions, loop over 4th dim and return an array
% of activation matrices?



% How to go from binary nii --> tSeries format (n_tpts x n_vox):
% 	1. Load ROI as nifti
% 		a. ENSURE ROI NIFTI AND DATA NIFTI ARE SAME ORI!!!!!
% 	2. Load dataset as nii
% 	3. Find ROI indices:
% 		a. find(myroi.data(:)==1)
% 	4. Reshape data in dataset to n_vox x n_timepoints:
% 		a. Dd = reshape(myd.data,prod(myd.dims(1:3)),myd.dim(4));
% Index into dd and transpose; done

% FIRST: check if dimesnsions (# pix) is same
if sum(nii_task.dim(1:3)==nii_roi.dim(1:3))~=3
    error('preproc:niftiExtract:mismatchDimensions','The number of voxels in ROI & data do not match');
end

% SECOND: check if size of voxels (first 3 dims, mm) is the same
if sum((nii_task.pixdim(1:3)-nii_roi.pixdim(1:3))<0.0001)~=3 % for floating point precision reasons...
    error('preproc:niftiExtract:mismatchVoxSize','The voxel size in ROI & data do not match');
end

% THIRD: check if each file is oriented the same (by comparing the qto_xyz
% matrix; could also comapre sign of this matrix)
% - if using AFNI processing pipeline, this should always be true, as the
%   nifti base used to fill ROIs is a RF data file. but good to be sure
if sum((nii_task.qto_xyz(:)-nii_roi.qto_xyz(:))<0.0001)~=16 % for floating point precision reasons...
    error('preproc:niftiExtract:incorrectOrientation','The ROI & task nii do not line up; check q/s transforms');
end


% FOURTH: if nii_task is 5dim, nifti-squeeze it
if length(nii_task.dim)>4
    fprintf('Input file %s has more than 4 dimensions; applying niftiSqueeze\n',nii_task.fname);
    nii_task = niftiSqueeze(nii_task);
end


data_mat = cell(size(nii_roi.data,4),1);

fprintf('%s\n',nii_roi.fname);
for ii = 1:size(nii_roi.data,4)
    
    %nvox = sum(reshape(nii_roi.data(:,:,:,ii),prod(nii_roi.dim(1:3)),1)>0);
    
    roiidx = find(nii_roi.data~=0);
    
    fprintf('ROI %i:\t%i vox\n',ii,length(roiidx));
    
    datatmp = reshape(nii_task.data,prod(nii_task.dim(1:3)),nii_task.dim(4));
    
    data_mat{ii} = datatmp(roiidx,:).';
    
    
end

% if only a single ROI found, return a matrix; otherwise return a cell
% array
if length(data_mat)==1
    data_mat_tmp = data_mat{1};
    clear data_mat;
    data_mat = data_mat_tmp;
end

return