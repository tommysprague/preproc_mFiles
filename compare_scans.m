% compare_scans.m
%
% alternates b/w two scans, two volumes, etc, makes a gif and saves it out
%
% for now, nii1, nii2 are nii structs loaded w/ niftiRead; should handle
% this more gracefully
%
% for now, doesn't change orientation of images, just takes them as they
% arrive
%
% TODO: allow for arbitrary # of nii's, so can make movies through entire
% session, etc. 
%
% (note that, on 7/30/2017, updated niftiPlotSlices.m - now does gray
% colormap, can plot multiple slices, orients to RAI, etc)
%
% Tommy Sprague



function fh = compare_scans(nii1,nii2,fout,vol1,vol2,coords,clims)

if nargin < 3 || isempty(fout)
    fout = 'myscans.gif';
end

if nargin < 4
    vol1 = 1;
end

if nargin < 5
    vol2 = 1;
end

if nargin < 6
%    coords = [];
    coords = round(size(nii1.data)/2);

end

if nargin < 7
    clims = [min(min(double(nii1.data(:))),min(double(nii2.data(:)))) max(max(double(nii1.data(:))),max(double(nii2.data(:))))];
end




%figure;

f1=niftiPlotSlices(nii1,coords,vol1,1);
set(get(gcf,'Children'),'CLim',clims);
colormap gray;
set(gcf,'Position',[571         939        1428         389]);
this_frame = getframe(f1);

im = frame2im(this_frame);
[imind,cm] = rgb2ind(im,256);

%close(f1);

% write first file
imwrite(imind,cm,fout,'gif', 'Loopcount',inf); 


f2=niftiPlotSlices(nii2,coords,vol2,1);
set(get(gcf,'Children'),'CLim',clims);
set(gcf,'Position',[571         939        1428         389]);
colormap gray;
this_frame = getframe(f2);
%close(f2);


im = frame2im(this_frame);
[imind,cm] = rgb2ind(im,256);


imwrite(imind,cm,fout,'gif','WriteMode','append'); 




return



% 
% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
% for n = 1:0.5:5
%     % Draw plot for y = x.^n
%     x = 0:0.01:1;
%     y = x.^n;
%     plot(x,y) 
%     drawnow 
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
%   end