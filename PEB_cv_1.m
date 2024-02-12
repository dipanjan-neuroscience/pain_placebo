%% This code runs PEB analysis for the connectivity Estimates obtained from DCM analysis %%

data_info = '/Users/sanjeev/Documents/Projects/Placebo_effect_pain_perception/extra '


data_sheet = xlsread(fullfile(data_info,'participants_modified_t.xlsx'));%
% Column 1: 0- Control; 1: study 1 patients; 2:Study 2 patients 
% colimne 2 : Gender:: 1 - Female ; 2 - Male 
% column 3 : age 
% column 4 : Drug 0 for controls, 1 for placebo, 2 for dulexotin 
% column 5 : -1 for controls, 0 for non-responsders, 1 for responsers

% From this will selectall 72 participants (healthy controls and Patients in study 1 & study 2)

% First column with all ones 
covariate = ones(length(data_sheet),1); 

for ii = 1:length(covariate)
    
    % second column for healthy control vs  patients
    
    if data_sheet(ii,1) == 0
        covariate(ii,2) = -1
    else 
        covariate(ii,2) = 1
    end 
    
    % third column for gender 
    
    if data_sheet(ii,2) == 1
        covariate(ii,3) = -1
    else 
        covariate(ii,3) = +1
    end 
    
    % Fifth column for age 
    
       
    
end 
covariate(:,4) = data_sheet(:,3); 
% normalize the age 
covariate(:,4)=covariate(:,4)-mean(covariate(:,4));  %mean centering the age column 


%% Run the analysis for right ROIs
gcm_dir = '/Users/sanjeev/Documents/Projects/Placebo_effect_pain_perception/Data/GLM2/GCM'

load(fullfile(gcm_dir,'GCM_R_three.mat'));

spm_dcm_fmri_check(GCM_R_three)

M   = struct();
M.Q = 'all'; 
% Specify design matrix for N subjects. It should start with a constant column
M.X= covariate;
M.maxit  = 256; % maximum number of iterations
% Choose field
field = {'A'};
% Estimate model
PEB_R_three    = spm_dcm_peb(GCM_R_three,M,field);
%save('PEB_R_three.mat','PEB_R_three'); 
BMA_R_three=spm_dcm_peb_bmc(PEB_R_three);
%save('BMA_R_three.mat','BMA_R_three');
spm_dcm_peb_review(BMA_R_three,GCM_R_three);
%% leave one out cross validation 
M   = struct();
M.Q = 'all'; 
% Specify design matrix for N subjects. It should start with a constant column
M.X = covariate;
[qE,qC,Q]=spm_dcm_loo(GCM_R_three,M,{'A(3,1)'});


%
%% Run the analysis for Left ROIs
gcm_dir = '/Users/sanjeev/Documents/Projects/Placebo_effect_pain_perception/Data/GLM2/GCM'

load(fullfile(gcm_dir,'GCM_L_three.mat'));

spm_dcm_fmri_check (GCM_L_three)

M   = struct();
M.Q = 'all'; 
% Specify design matrix for N subjects. It should start with a constant column
M.X= covariate;
M.maxit  = 256; % maximum number of iterations
% Choose field
field = {'A'};
% Estimate model
PEB_L_three    = spm_dcm_peb(GCM_L_three,M,field);
%save('PEB_R_three.mat','PEB_R_three'); 
BMA_L_three=spm_dcm_peb_bmc(PEB_L_three);
%save('BMA_R_three.mat','BMA_R_three');
spm_dcm_peb_review(BMA_L_three,GCM_L_three);
%% leave one out cross validation 
M   = struct();
M.Q = 'all'; 
% Specify design matrix for N subjects. It should start with a constant column
M.X = covariate;

[qE,qC,Q]=spm_dcm_loo(GCM_L_three,M,{'A(3,1)'}); %(to,from)

