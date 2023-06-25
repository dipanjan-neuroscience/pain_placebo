%initialize spm
root_dir = pwd;
SPM_PATH = fullfile(root_dir, 'spm12');
addpath(SPM_PATH)
spm('Defaults','fMRI');
spm_jobman('initcfg');


%% extract voi

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

clear  subFolders 

GLM2_dir=fullfile(root_dir,'GLM2');

%collect left mask file names
files = dir('Left');
files(ismember({files.name},{'.','..'}))=[]; %remove . and ..

for k = 1 : length(files)
lMaskNames{k}=files(k).name;
end

[pathstr, lMaskNames, ext] = fileparts(lMaskNames);
[pathstr, lVoiNames, ext] = fileparts(lMaskNames);

%maskNames={'lFp1.nii','lSSC.nii','lV1.nii','lA1.nii','lSMA.nii','lMC.nii','laIns.nii','lpIns.nii','lThal.nii','lBroca.nii','rFp1.nii','rSSC.nii','rV1.nii','rA1.nii','rSMA.nii','rMC.nii','raIns.nii','rpIns.nii'};
%voiNames={'lFp1','lSSC','lV1','lA1','lSMA','lMC','laIns','lpIns','lThal','lBroca','rFp1','rSSC','rV1','rA1','rSMA','rMC','raIns','rpIns'};


for sI = 1: length(subNames)

spm_dir = fullfile('GLM2_dir', subNames{sI},'SPM.mat');

for vI = 1: length(voiNames)

mask_dir= fullfile (root_dir, lMaskNames{vI});


clear matlabbatch;
matlabbatch{1}.spm.util.voi.spmmat = cellstr(spm_dir); % directory to spm.mat file
matlabbatch{1}.spm.util.voi.adjust = NaN;
matlabbatch{1}.spm.util.voi.session = 1; % Session index
matlabbatch{1}.spm.util.voi.name = lVoiNames{vI}; % name you want to give 
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(mask_dir); % directory to mask image
%matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.mtype = 0; % inclusion
matlabbatch{1}.spm.util.voi.expression = 'i1';
spm_jobman('run',matlabbatch);



end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% calculate correlation matrix among the time series
[lr,lp] = corrcoef(t{:,:});
%% Fisher's Z transformation of correlation values
% The Fisher Z-Transformation is a way to transform the sampling distribution of Pearson’s r 
% (i.e. the correlation coefficient) 
% so that it becomes normally distributed. The “z” in Fisher Z stands for a z-score.

lz=atanh(lr);

%%add to the multidimensional array
lzmatrix(:,:,sI)=lz;

cd ..
cd ..

end


%% Two smaple t test

[h,p,ci,stats] = ttest2(lzmatrix(:,:,m),lzmatrix(:,:,n),'Dim',3);
lzmatrix=lzmatrix.*h;

%% inverse transformation

lrmatrix=tanh(zmatrix);


