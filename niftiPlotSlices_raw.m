% Plot 'raw' nii slices - don't do any adjustments, just plot the block of
% data. Useful for checking transforms. Now, niftiPlotSlices.m uses some
% transfomrations to make data a bit more easily oriented for making
% alignment comparison gifs, etc.
%
% TCS 7/30/2017


function fig_pointer = niftiPlotSlices_raw(nii,coord,volnum,transpose_mode)

if nargin < 2
    coord = round(size(nii.data)/2);
end

if nargin < 3
    volnum = 1;
end

if nargin < 4
    transpose_mode = 0; % 0: plot exactly the images, in matlab i,j coords; o
                        % 1: plot transposed/flip90'd versions
end

fig_pointer = figure;

subplot(1,3,1);
if transpose_mode == 1
    imagesc(squeeze(nii.data(:,:,coord(3),volnum)).');
else
    imagesc(squeeze(nii.data(:,:,coord(3),volnum)));
end
axis equal tight xy;
ylabel('dim 1');
xlabel('dim 2');
set(gca,'XTick',[],'YTick',[]);
title(sprintf('Slice %i',coord(3)));


subplot(1,3,2);
if transpose_mode==1
    imagesc(squeeze(nii.data(:,coord(2),:,volnum)).');
else
    imagesc(squeeze(nii.data(:,coord(2),:,volnum)));
end
axis equal tight xy;
ylabel('dim 1');
xlabel('dim 3');
set(gca,'XTick',[],'YTick',[]);
title(sprintf('Slice %i',coord(2)));


subplot(1,3,3);
if transpose_mode==1
    imagesc(squeeze(nii.data(coord(1),:,:,volnum)).');
else
    imagesc(squeeze(nii.data(coord(1),:,:,volnum)));
end
axis equal tight xy;
ylabel('dim 2');
xlabel('dim 3');
set(gca,'XTick',[],'YTick',[]);
title(sprintf('Slice %i',coord(1)));




return