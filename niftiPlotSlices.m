% Plot slices through defined coord(s) from a nifti file, loaded from
% niftiRead (requires vistasoft). Uses vistasoft orientation commands to
% try and orient the image sensically (up is up, anterior). Multiple coords
% can be provided, each row will correspond to one set of coords. Used by
% compare_scans.m, which creates gifs


function fig_pointer = niftiPlotSlices(nii,coords,volnum)

if nargin < 2
    coords = round(size(nii.data)/2);
end

if nargin < 3
    volnum = 1;
end

transpose_mode =1;

% if nargin < 4
%     transpose_mode = 0; % 0: plot exactly the images, in matlab i,j coords; o
%                         % 1: plot transposed/flip90'd versions
% end


% TODO: check coords, make sure they're inside vol dimensions...
% TODO: if coords is a single #, make this many even slices through each
% dimension


% if a filename given, load the nii
if ischar(nii)
    nii = niftiRead(nii);
end

curr_ori = niftiCurrentOrientation(nii);

desired_ori = 'LPS'; % NOTE: opposite reporting of afni...

tmpnii = niftiApplyXform(nii,niftiCreateXformBetweenStrings(curr_ori,desired_ori));  % reorient the nii so it's in the 'matched' ori as the model (PRS)

nii = tmpnii;


fig_pointer = figure;

for cc = 1:size(coords,1)
    subplot(size(coords,1),3,(cc-1)*3+1);
    if transpose_mode == 1
        imagesc(squeeze(nii.data(:,:,coords(cc,3),volnum)).');
    else
        imagesc(squeeze(nii.data(:,:,coords(cc,3),volnum)));
    end
    axis equal tight xy;
    ylabel('dim 1');
    xlabel('dim 2');
    set(gca,'XTick',[],'YTick',[]);
    title(sprintf('Slice %i',coords(cc,3)));
    
    
    subplot(size(coords,1),3,(cc-1)*3+2);
    if transpose_mode==1
        imagesc(squeeze(nii.data(:,coords(cc,2),:,volnum)).');
    else
        imagesc(squeeze(nii.data(:,coords(cc,2),:,volnum)));
    end
    axis equal tight xy;
    ylabel('dim 1');
    xlabel('dim 3');
    set(gca,'XTick',[],'YTick',[]);
    title(sprintf('Slice %i',coords(cc,2)));
    
    
    subplot(size(coords,1),3,(cc-1)*3+3);
    if transpose_mode==1
        imagesc(squeeze(nii.data(coords(cc,1),:,:,volnum)).');
    else
        imagesc(squeeze(nii.data(coords(cc,1),:,:,volnum)));
    end
    axis equal tight xy;
    ylabel('dim 2');
    xlabel('dim 3');
    set(gca,'XTick',[],'YTick',[]);
    title(sprintf('Slice %i',coords(cc,1)));
end

all_clim = get(get(fig_pointer,'Children'),'CLim');
all_clim = cell2mat(all_clim);
set(get(fig_pointer,'Children'),'CLim',[min(all_clim(:,1)) max(all_clim(:,2))]);
colormap gray;


return