% Get a list of all files and folders in ts folder.
files = dir('ts');
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
for k = 1 : length(subFolders)
subNames{1,k}=subFolders(k).name;
end

clear  subFolders 



for sI=1: length(subNames)


cd(fullfile('ts\', subNames{sI}))    
%%convert csv file to matlab table
t = readtable("out_parc_timeseries.tsv", "FileType","text",'Delimiter', '\t');
% calculate correlation matrix among the time series
[r,p] = corrcoef(t{:,:});
%% Fisher's Z transformation of correlation values
% The Fisher Z-Transformation is a way to transform the sampling distribution of Pearson’s r 
% (i.e. the correlation coefficient) 
% so that it becomes normally distributed. The “z” in Fisher Z stands for a z-score.
z=atanh(r);



%%add to the multidimensional array
zmatrix(:,:,sI)=z;

cd ..
cd ..

end


%% Two smaple t test

[h,p,ci,stats] = ttest2(zmatrix(:,:,m),zmatrix(:,:,n),'Dim',3);
zmatrix=zmatrix.*h;

%% inverse transformation

rmatrix=tanh(zmatrix);

%% plotting result on the brain
