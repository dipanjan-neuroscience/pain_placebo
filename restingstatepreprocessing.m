% code for running bath analysis of resting state data in SPM12
% Things to be done before running the code: 
%1. Set the origin of the anatomical and functional images to anterior commissure ( see http://andysbrainblog.blogspot.com/2012/11/spm-setting-origin-and-normalization.html for rationale and the way to do it)
% by using the script acpc_coreg.m

clear all

root_dir = pwd;
SPM_PATH = fullfile(root_dir, 'spm12');
addpath(SPM_PATH)


%% Initialize SPM

spm('Defaults','fMRI');
spm_jobman('initcfg');

spm('CreateIntWin','on');
spm_figure('Create','Graphics','Graphics','on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition according to the specific data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

nd=5; % no of dummy scans

% for Slice time correction 
nslices = 25; 
tr = 2.5; % in seconds
ta = 2.4; % ta = tr - (tr/nslices)
so = [0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400]; % here we have ebtered slice timing in milliseconds
refslice = 0; % reference slice timing in milliseconds


% for smoothing
fwhm=[4 4 10]; % the thumb of rules says fwhm should be twice the voxel dimension

% Get a list of all files and folders in func folder.
files = dir('func');
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
for k = 1 : length(subFolders)
subNames{1,k}=subFolders(k).name;
end

clear k subFolders 

for sI = 1: length(subNames)

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% directories 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% define directories    
str_dir = fullfile(root_dir, subNames{sI},'t1');
func_dir = fullfile(root_dir, subNames{sI},'rsfmri');


% file select
f_or = spm_select('FPList',func_dir,'^vol.*\.nii$'); % original functional images
%f = spm_select('FPList',func_dir,'^vol.*\.nii$'); % functional images
s= spm_select('FPList',str_dir,'^defaced.*\.nii$'); % structural images 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create directory GLM & GLM2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
GLM_dir=fullfile(root_dir,'GLM');
mkdir(GLM_dir)
glm_dir=fullfile(GLM_dir, subNames{sI});
mkdir(glm_dir)

GLM2_dir=fullfile(root_dir,'GLM2');
mkdir(GLM2_dir)
glm2_dir=fullfile(GLM2_dir, subNames{sI});
mkdir(glm2_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONVERT FUNCTIONAL SCANS FROM 4D TO 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spm_file_split(f_or);



% clear matlabbatch
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(func_dir);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'original';
 
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(f_or); % remove the original functional file
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = cellstr(fullfile(func_dir,'original'));
 
spm_jobman('run',matlabbatch);


f = spm_select('FPList',func_dir,'^vol.*\.nii$'); % functional images

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DUMMY SCANS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(func_dir);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'dummy';

matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(f(1:nd,:)); % remove first 5 scans to allow for magnetization to be stable
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = cellstr(fullfile(func_dir,'dummy'));

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPATIAL PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Realignment
f = spm_select('FPList',func_dir,'^vol.*\.nii$'); % functional images

clear matlabbatch

matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(f)};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0]; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = ''; 
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1]; 
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4; 
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0]; 
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1; 
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r'; 
spm_jobman('run', matlabbatch)


%% creating friston24 regressor

rp = spm_select('FPList',func_dir,'^rp.*\.txt$'); % 6 head motion parameters
rp=readmatrix(rp);
[r,c]=size(rp);
rp1=vertcat(zeros(1,6),rp); % 6 head motion parameters from the previous time point
rp1(r+1,:)=[];
sq=rp.^2; % squared head motion parameters

% file select
rf = spm_select('FPList',func_dir,'^rvol.*\.nii$');% realigned images
meanf=spm_select('FPList',func_dir,'^meanvol.*\.nii$');% mean image

clear matlabbatch

%% Slice time correction
           
matlabbatch{1}.spm.temporal.st.scans =  {cellstr(rf)};
matlabbatch{1}.spm.temporal.st.nslices = nslices; 
matlabbatch{1}.spm.temporal.st.tr = tr; 
matlabbatch{1}.spm.temporal.st.ta = ta; 
matlabbatch{1}.spm.temporal.st.so = so; 
matlabbatch{1}.spm.temporal.st.refslice = refslice; 
matlabbatch{1}.spm.temporal.st.prefix = 'a'; 


%% Coregistration

matlabbatch{2}.spm.spatial.coreg.estimate.ref    = cellstr(meanf);
matlabbatch{2}.spm.spatial.coreg.estimate.source = cellstr(s);
matlabbatch{2}.spm.spatial.coreg.estimate.other = {''}; 
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi'; 
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2]; 
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]; 
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7]; 


%% Segmentation 
matlabbatch{3}.spm.spatial.preproc.channel.vols = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles')); 
matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.001; 
matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 60; 
matlabbatch{3}.spm.spatial.preproc.channel.write = [1 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm =cellstr(fullfile(SPM_PATH,'tpm\TPM.nii,1'));
matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 1; 
matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [1 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm =cellstr(fullfile(SPM_PATH,'tpm\TPM.nii,2')); 
matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 1; 
matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [1 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm = cellstr(fullfile(SPM_PATH,'tpm\TPM.nii,3')); 
matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2; 
matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [1 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm = cellstr(fullfile(SPM_PATH,'tpm\TPM.nii,4')); 
matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 3; 
matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [1 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm = cellstr(fullfile(SPM_PATH,'tpm\TPM.nii,5')); 
matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4; 
matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [1 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm = cellstr(fullfile(SPM_PATH,'tpm\TPM.nii,6')); 
matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 2; 
matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0]; 
matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0]; 
matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1; 
matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 1; 
matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2]; 
matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni'; 
matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0; 
matlabbatch{3}.spm.spatial.preproc.warp.samp = 3; 
matlabbatch{3}.spm.spatial.preproc.warp.write = [0 1]; 
matlabbatch{3}.spm.spatial.preproc.warp.vox = NaN; 
matlabbatch{3}.spm.spatial.preproc.warp.bb = [NaN NaN NaN 
NaN NaN NaN]; 




%% Normalisation
matlabbatch{4}.spm.spatial.normalise.write.subj.def = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'})); 
matlabbatch{4}.spm.spatial.normalise.write.subj.resample = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files')); 
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 
78 76 85]; 
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [2 2 2]; 
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4; 
matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix = 'w'; 





%% Smoothing
matlabbatch{5}.spm.spatial.smooth.data = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files')); 
matlabbatch{5}.spm.spatial.smooth.fwhm = fwhm; 
matlabbatch{5}.spm.spatial.smooth.dtype = 0; 
matlabbatch{5}.spm.spatial.smooth.im = 0; 
matlabbatch{5}.spm.spatial.smooth.prefix = 's'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM SPECIFICATION, ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Model specification
matlabbatch{6}.spm.stats.fmri_spec.dir = cellstr(fullfile(GLM_dir, subNames{sI}));
matlabbatch{6}.spm.stats.fmri_spec.timing.units = 'scans'; 
matlabbatch{6}.spm.stats.fmri_spec.timing.RT = tr; 
matlabbatch{6}.spm.stats.fmri_spec.timing.fmri_t = 16; 
matlabbatch{6}.spm.stats.fmri_spec.timing.fmri_t0 = 8; 
matlabbatch{6}.spm.stats.fmri_spec.sess.scans = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files')); 
matlabbatch{6}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}); 
matlabbatch{6}.spm.stats.fmri_spec.sess.multi = {''}; 
matlabbatch{6}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {}); 
matlabbatch{6}.spm.stats.fmri_spec.sess.multi_reg = {''}; 
matlabbatch{6}.spm.stats.fmri_spec.sess.hpf = 128; 
matlabbatch{6}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {}); 
matlabbatch{6}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; 
matlabbatch{6}.spm.stats.fmri_spec.volt = 1; 
matlabbatch{6}.spm.stats.fmri_spec.global = 'None'; 
matlabbatch{6}.spm.stats.fmri_spec.mthresh = -Inf; 
matlabbatch{6}.spm.stats.fmri_spec.mask = {''}; 
matlabbatch{6}.spm.stats.fmri_spec.cvi = 'AR(1)'; 


%% Estimation
matlabbatch{7}.spm.stats.fmri_est.spmmat = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat')); 
matlabbatch{7}.spm.stats.fmri_est.write_residuals = 0; 
matlabbatch{7}.spm.stats.fmri_est.method.Classical = 1; 
spm_jobman('run', matlabbatch)
%%

% file select
mask = spm_select('FPList',glm_dir,'^mask.*\.nii$'); % whole brain mask


clear matlabbatch


%% Creation of  white matter and CSF regressor
matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(GLM_dir, subNames{sI},'SPM.mat'));
matlabbatch{1}.spm.util.voi.adjust = NaN; 
matlabbatch{1}.spm.util.voi.session = 1; 
matlabbatch{1}.spm.util.voi.name = 'WM'; 
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [0 -24 -33]; 
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 6; 
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1; 
matlabbatch{1}.spm.util.voi.roi{2}.mask.image(1) = cellstr(mask);
matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5; 
matlabbatch{1}.spm.util.voi.expression = 'i1&i2'; 
spm_jobman('run',matlabbatch);

wm=Y;

clear matlabbatch


matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(GLM_dir, subNames{sI},'spm.mat'));
matlabbatch{1}.spm.util.voi.adjust = NaN; 
matlabbatch{1}.spm.util.voi.session = 1; 
matlabbatch{1}.spm.util.voi.name = 'CSF'; 
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [0 -40 -5]; 
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 6; 
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1; 
matlabbatch{1}.spm.util.voi.roi{2}.mask.image(1) = cellstr(mask); 
matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5; 
matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
spm_jobman('run', matlabbatch)


csf=Y;

clear matlabbatch


nuis_reg= horzcat(rp,rp1,sq,wm,csf); % friston regressor plus time series from wm and csf as nuissance regressor
save(fullfile(glm_dir,'nuis_reg.txt' ), 'nuis_reg','-ascii');

% file select
smooth= spm_select('FPList',func_dir,'^swarvol.*\.nii$'); % select smoothed files

%% Model specification
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(fullfile(GLM2_dir, subNames{sI}));
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans'; 
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.5; 
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16; 
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; 
matlabbatch{1}.spm.stats.fmri_spec.sess.scans=cellstr(smooth);
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}); 
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''}; 
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {}); 
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(glm_dir,'nuis_reg.txt' )};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; 
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {}); 
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; 
matlabbatch{1}.spm.stats.fmri_spec.volt = 1; 
matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; 
matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf; 
matlabbatch{1}.spm.stats.fmri_spec.mask = {''}; 
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; 


%% Estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat')); 
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0; 
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; 
%% 


spm_jobman('run', matlabbatch)


end
 
