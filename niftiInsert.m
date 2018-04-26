function newnii = niftiInsert(mask_nii,data_mat)
% inserts data into a nifti shaped like mask_nii from data_mat
%
% MASK_NII is a nii struct (like niftiRead) that's 1's and 0's. 
% DATA_MAT is a matrix that's n_vols x n_vox (number of non-zero entries in
%   mask_nii)
%
% This is essentially the opposite of niftiExtract.m
%
% Tommy Sprague, 11/5/2017 (draft)

newnii = mask_nii;

% first reshape mask to be 1 x n_vox
mask_reshape = reshape(mask_nii.data,1,prod(mask_nii.dim(1:3)));
data_all = zeros(size(data_mat,1),numel(mask_reshape));
data_all(:,mask_reshape==1) = data_mat;

newnii.data = reshape(data_all.',newnii.dim(1),newnii.dim(2),newnii.dim(3),size(data_mat,1));
newnii.dim(4) = size(newnii.data,4);



return