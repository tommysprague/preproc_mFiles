function fig_pointer = niftiPlotSlices(nii,coord,volnum)

if nargin < 2
    coord = round(size(nii.data)/2);
end

if nargin < 3
    volnum = 1;
end

fig_pointer = figure;

subplot(1,3,1);
imagesc(squeeze(nii.data(:,:,coord(3),volnum)));
axis equal tight;
ylabel('dim 1');
xlabel('dim 2');
title(sprintf('Slice %i',coord(3)));


subplot(1,3,2);
imagesc(squeeze(nii.data(:,coord(2),:,volnum)));
axis equal tight;
ylabel('dim 1');
xlabel('dim 3');
title(sprintf('Slice %i',coord(2)));


subplot(1,3,3);
imagesc(squeeze(nii.data(coord(1),:,:,volnum)));
axis equal tight;
ylabel('dim 2');
xlabel('dim 3');
title(sprintf('Slice %i',coord(1)));




return