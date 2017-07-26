% niftiCompareTimeSeries.m
%
% From 2 nii structs (or filenames, will load nifti's), and set of coords
% (1-indexed), plot both timeseries overlapping. 
%
% TODO: full grid like in afni, etc
%
% TODO: if either arg is a cell array of strings or nii structs, average &
% plot SE
%
% NOTE: detrending time series before plotting!!!!!! 
%
% TCS, 7/5/2017


function fig_handle = niftiCompareTimeSeries(nii1,nii2,coords)


% TODO: if not already opened
fig_handle = figure;

for cc = 1:size(coords,1)
    subplot(size(coords,1),1,cc);
    hold on;
    
    
    % TODO: always detrend? will make things easier....for now, yes
    plot(nii1.pixdim(4)*(0:size(nii1.data,4)-1), detrend( squeeze(nii1.data(coords(cc,1),coords(cc,2),coords(cc,3),:)) ) );
    
    plot(nii2.pixdim(4)*(0:size(nii2.data,4)-1), detrend( squeeze(nii2.data(coords(cc,1),coords(cc,2),coords(cc,3),:)) ) );
    
    xlabel('Time (s)');
    ylabel('BOLD response (detrended)');
    
    title(sprintf('coord: %i, %i, %i',coords(cc,1),coords(cc,2),coords(cc,3)));
    
    xlim([0 max((nii1.pixdim(4)*size(nii1.data,4)-1),(nii2.pixdim(4)*size(nii2.data,4)-1))])
    
    
end

hold off;


return