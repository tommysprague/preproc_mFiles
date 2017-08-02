%% do_RFs.m
%
% general function for computing RFs, given averaged EPIs.
%
% call from command line like:
% cd /my/func/dir
% matlab -nodesktop -nojvm -nosplash -r "do_RFs(arg1,arg2,arg3);quit();"
%
% NECESSARY ARGUMENTS:
% - SUBJ ID
% - session path (where the vista session will live)
% - EPIext: suffix used to identify EPIs (after bar_width_*)
%   of files as cell array)
% - stimExt (where to look for bar width files)
%
% TCS, 3/17/2017
% TCS, 7/30/2017 - improvements, like dynamically reading # TRs, cleaning up
% old bits

% optional ROIs arg?


function do_RFs_mb_tst(subjID,sessPath,EPIext,IPname,stimExt,myTR)

addpath ~/
startup;    % because sometimes running from command line is weird?


%% get all inputs in order

% if EPIext is an array, loop over those,?


if nargin < 3
    EPIext = {'bc','ss5','surf'};%.nii.gz';
end

% because we want to loop over EPI extensions
if ~iscell(EPIext)
    EPIext = {EPIext};
end

% if no IPname, use the surfanat_brainmask.nii.gz made by spatial_tcs.sh,
% which will be in same space/resolution as EPInames

% maybe use a 1.25 mm IP?
if nargin < 4 || isempty(IPname)
    IPname = '../../surfanat_brainmask_hires.nii.gz';
end

if nargin < 5
    stimExt = '';
end

if nargin < 6
    myTR = 1.2; % typically use this (s)
end

if ~isnumeric(myTR)
    myTR = str2num(myTR);
end

% outsourced this stuff to a separate command (setup_paths_RFs.m, which
% will live in this same directory and is run from matlab command)
%
% or maybe this happens here? sure, let's do that
setup_paths_RFs;


%% initialize vista session

cd(sessPath);


for ee = 1:length(EPIext)
    
    
    
    % Specify EPI files
    tmp = dir(sprintf('bar_width_*%s.nii.gz', EPIext{ee}));
    for ii = 1:numel(tmp)
        epi_file{ii} = tmp(ii).name;
        assert(exist(epi_file{ii}, 'file')>0)
        
    end
    
    
    % FOR NOW: assume all same # of TRs, so get nTRs from first nii
    % file
    niitmp = niftiRead(epi_file{ii});
    nTRs = niitmp.dim(4);
    fprintf('NII file %s has %i TRs\n',niitmp.fname,nTRs);
    clear niitmp;
    
    % Specify INPLANE file
    inplane_file = IPname;
    assert(exist(inplane_file, 'file')>0)
    
    % FOR NOW: assumign these are copied from FS direcotry, converted to NII,
    % and live here - going to see if we can gt by w/ anat_T1_brain.nii, but
    % the ribbon.mgz-sourced t1_class.nii.gz; if not, we'll go back to using
    % t1.nii.gz, which isn't being used right now
    % Specify 3DAnatomy file
    anat_file = '../../anat_T1_brain.nii';
    assert(exist(anat_file, 'file')>0)
    
    % Generate the expected generic params structure
    params = mrInitDefaultParams;
    
    % And insert the required parameters:
    params.inplane      = inplane_file;
    params.functionals  = epi_file;
    params.sessionDir   = sessPath;
    
    params.annotations = epi_file; % TCS 10/10/2016 - let's try this (getting errors)
    
    % Specify some optional parameters
    params.vAnatomy = anat_file;
    params.keepFrames = []; %We dropped frames in pre-processing
    params.subject = subjID;
    
    for ii = 1:length(epi_file)
        params.coParams{ii} = coParamsDefault;
        params.coParams{ii}.nCycles = 8;
    end
    
    
    % Run it:
    ok = mrInit(params);
    
    global dataTYPES;

    
    
    
    %% run RF model
    
    ip = initHiddenInplane;
    
    % initialize
    fnames = dir(sprintf('%s/Stimuli/*.mat',sessPath));
    if isempty(fnames), error('We need a file of images and image sequences within Stimuli directory'); end
    
    params(1).fliprotate=[0 0 0]; %% Set up pRF model
    params(1).stimType='StimFromScan';
    params(1).stimSize=29.35/2; % DEGREES VISUAL ANGLE (RADIUS)
    params(1).stimWidth=45; % ignored when you have 'stimFromScan'
    params(1).stimStart=0;  % ignored when you have 'stimFromScan'
    params(1).stimDir=0;    % ignored when you have 'stimFromScan'
    params(1).nCycles=6;    % ignored when you have 'stimFromScan'
    params(1).nStimOnOff=0; % ignored when you have 'stimFromScan'
    params(1).nUniqueRep=1; % ignored when you have 'stimFromScan'
    params(1).prescanDuration=0;
    params(1).nDCT=1;      % max frequency in in detrending using discrete cosine transform (1 = DC, plus 1/2 cycle, plus 1 cycle)
    params(1).hrfType='two gammas (SPM style)';
    params(1).hrfParams={[1.6800 3 2.0500]  [5.4000 5.2000 10.8000 7.3500 0.3500]};
    params(1).calcPC = 0; % tryign this out...
    
    
    % TODO: be smarter at this!!!!
    params(1).framePeriod=myTR; % TR in seconds TODO: make this automatic using nii info!!!!!
    params(1).nFrames=nTRs;        % number of volumes or 304
    
    params(1).imFile=sprintf('%s/Stimuli/bar_stimulus_masks_1200ms_images%s.mat',sessPath,stimExt);     % see makeStimFromScan
    params(1).jitterFile=sprintf('%s/Stimuli/none',sessPath);                 % ignore (for eye movements)
    params(1).paramsFile=sprintf('%s/Stimuli/bar_stimulus_masks_1200ms_params%s.mat',sessPath,stimExt); % see makeStimFromScan
    params(1).imFilter='none';                           % stim file is already a binary contrast mask
    
    %params(2) = params(1);
    %params(2).imFile = sprintf('%s/Stimuli/barWidth2_images%s.mat',sessPath,stimExt);
    %params(2).paramsFile = sprintf('%s/Stimuli/barWidth2_params%s.mat',sessPath,stimExt);
    
    %params(3) = params(1);
    %params(3).imFile = sprintf('%s/Stimuli/barWidth3_images%s.mat',sessPath,stimExt);
    %params(3).paramsFile = sprintf('%s/Stimuli/barWidth3_params%s.mat',sessPath,stimExt);
    
    %ip  = viewSet(ip,'rmParams',params);
    ip = viewSet(ip, 'Current DataTYPE', 'Original');
    
    dt = viewGet(ip,'dt struct');
    dt = dtSet(dt,'rmparams',params);
    %dataTYPES(1).retinotopyModelParams = params;
    dataTYPES(1).retinotopyModelParams=dt.retinotopyModelParams; % hack...trouble setting struct w/ new fields?
    % store it
    
    saveSession;
    mrvCleanWorkspace;
    clear ip dt;
    %
    
    % put the rm params into the view structure
    %ip = rmLoadParameters(ip);
    
    % check it
    %rmStimulusMatrix(viewGet(vw, 'rmparams'), [], [], 2, false);
    
    
    % Run it
    ip = initHiddenInplane;
    
    ip = viewSet(ip, 'Current DataTYPE', 'Original');
    
    ip = rmLoadParameters(ip);
    
    % Define scan/stim and analysis parameters
    params = rmDefineParameters(ip, 'model', {'onegaussiannonlinear_gpu'});
    %params = rmDefineParameters(ip, 'model', {'onegaussiannonlinear'});
    
    % make stimulus and add it to the parameters
    params = rmMakeStimulus(params);
    
    params.analysis.calcPC = 0; % need to echo it here????
    
    % store params in view struct
    ip  = viewSet(ip,'rmParams',params);
    
    
    RF_fn = sprintf('RF_%s',EPIext{ee});
    
    %ip = rmMain(ip, 'gray', 'coarse to fine', 'model', {'onegaussiannonlinear_gpu'},'matFileName',sprintf('wm_nonlinear_gray_bp_noDecimate_noC2F'),'coarseDecimate',0,'coarseToFine',0);
    ip = rmMain(ip, [], 'coarse to fine', 'model', {'onegaussiannonlinear_gpu'},'matFileName',RF_fn,'coarseDecimate',0,'coarseToFine',0,'calcPC',0);
    %ip = rmMain(ip, [], 'coarse to fine', 'model', {'onegaussiannonlinear'},'matFileName',RF_fn,'coarseDecimate',0,'coarseToFine',0,'calcPC',0);
    % store it
    saveSession;
    mrvCleanWorkspace;
    clear ip;
    
    
    
    
    %% save results into nii's
    basenii = niftiRead(epi_file{1});
    which_params = {'pol','ve','ecc','sigmamajor','exponent','x0','y0','b'};
    which_models = {'gFit','sFit','fFit'};
    for mm = 1:length(which_models)
        rf_model = load(fullfile(sessPath,'Inplane','Original',sprintf('%s-%s.mat', RF_fn,which_models{mm})));
        
        niftiWriteRF(rf_model.model{1},fullfile(sessPath,'Inplane','Original',sprintf('%s-%s.nii.gz', RF_fn,which_models{mm})), basenii, which_params);
        
        clear rf_model;
    end
    
    
end
