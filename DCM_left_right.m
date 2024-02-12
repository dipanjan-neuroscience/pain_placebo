%% define DCMs

clc
clear all
close all

cd( '/Users/sanjeev/Documents/Projects/Placebo_effect_pain_perception/Data')
clear all
root_dir = '/Users/sanjeev/Documents/Projects/Placebo_effect_pain_perception/Data';
SPM_PATH = fullfile(root_dir, 'spm12');
%addpath(SPM_PATH)
spm('Defaults','fMRI');
spm_jobman('initcfg');

%% 
% Get a list of all files and folders in func folder.
files = dir('sub*');
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
for k = 1 : length(subFolders)
subNames{1,k}=subFolders(k).name;
end

subNames(:,[2 3 26 73]) = []

%% 

GLM2_dir = ('/Users/sanjeev/Documents/Projects/Placebo_effect_pain_perception/Data/GLM2');
%voiNamesL={'VOI_lsth_1.mat','VOI_lSSC_1.mat','VOI_lSII_1.mat','VOI_lpIns_1.mat','VOI_laIns_1.mat','VOI_lACC_1.mat','VOI_lFP1_1.mat'};
%voiNamesR={'VOI_rsth_1.mat','VOI_rSSC_1.mat','VOI_rSII_1.mat','VOI_rpIns_1.mat','VOI_raIns_1.mat','VOI_rACC_1.mat','VOI_rFP1_1.mat'};

voiNamesL={'VOI_lpIns_1.mat','VOI_laIns_1.mat'};
voiNamesR={'VOI_rpIns_1.mat','VOI_raIns_1.mat'};


for sI = 1: length(subNames)

cd(fullfile(GLM2_dir, subNames{sI}));

model_name = 'L_two' ;% 'L_seven';

xY         = voiNamesL;

SPM        = 'SPM.mat';

n   = 2;    % number of regions

nu  = 1;    % number of inputs. For DCM for CSD we have one input: null

TR  = 2.5;    % volume repetition time (seconds)

TE  = 0.03; % echo time (seconds)

 

% Connectivity matrices

a  = ones(n,n);

b  = zeros(n,n,nu);

c  = zeros(n,nu);

d  = zeros(n,n,0);

 

% Specify DCM

s = struct();

s.name       = model_name;

s.u          = [];

s.delays     = repmat(TR/2, 1, n)';

s.TE         = TE;

s.nonlinear  = false;

s.two_state  = false;

s.stochastic = false;

s.centre     = false;

s.induced    = 1;       % indicates DCM for CSD

s.a          = a;

s.b          = b;

s.c          = c;

s.d          = d;

 
DCM = spm_dcm_specify(SPM,xY,s);

end

clear DCM

for sI = 1: length(subNames)

cd(fullfile(GLM2_dir, subNames{sI}));

model_name = 'R_two' ;% 'R_seven';

xY         = voiNamesR;

SPM        = 'SPM.mat';

n   = 2;    % number of regions

nu  = 1;    % number of inputs. For DCM for CSD we have one input: null

TR  = 2.5;    % volume repetition time (seconds)

TE  = 0.03; % echo time (seconds)

 

% Connectivity matrices

a  = ones(n,n);


b  = zeros(n,n,nu);

c  = zeros(n,nu);

d  = zeros(n,n,0);

% Specify DCM

s = struct();

s.name       = model_name;

s.u          = [];

s.delays     = repmat(TR/2, 1, n)';

s.TE         = TE;

s.nonlinear  = false;

s.two_state  = false;

s.stochastic = false;

s.centre     = false;

s.induced    = 1;       % indicates DCM for CSD

s.a          = a;

s.b          = b;

s.c          = c;

s.d          = d;

 
DCM = spm_dcm_specify(SPM,xY,s);



end

clear DCM

%%estimate DCMs

for h=1: length(subNames) 
   
%GCM_L_seven(h,1) = {fullfile(GLM2_dir, subNames{h},'DCM_L_seven.mat')}; 
 
GCM_L_two(h,1) = {fullfile(GLM2_dir, subNames{h},'DCM_L_two.mat')}; 
 
end
  
for h=1: length(subNames) 
   
GCM_R_two(h,1) = {fullfile(GLM2_dir, subNames{h},'DCM_R_two.mat')}; 
 
end

use_parfor = true ;

cd('/Users/sanjeev/Documents/Projects/Placebo_effect_pain_perception/Data/GLM2/GCM');
GCM_L_two = spm_dcm_fit(GCM_L_two);

save('GCM_L_two_pIns_aIns.mat','GCM_L_two');


GCM_R_two = spm_dcm_fit(GCM_R_two);
save('GCM_R_two_pIns_aIns.mat','GCM_R_two');