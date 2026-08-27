%% LINEAR MODEL IMAGING vs. EEG
% Bahne H Bahners, BWH 2024
clear variables
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
rmappath=dryeeg_result_path(repoRoot, 'primary_crossvalidations_leftchannels');
savepath=dryeeg_result_path(repoRoot, 'imaging_imagingvseeg_linearmodel');
if~exist(savepath,'file')
mkdir(savepath)
end

eeg=load([rmappath,'Results_leftchannels_','LOO-CV','.mat'])%,'RegressorHat','Regressor','R','p');
eeg10=load([rmappath,'Results_leftchannels_','10-fold-CV','.mat'])%,'RegressorHat','Regressor','R','p');
eeg5=load([rmappath,'Results_leftchannels_','5-fold-CV','.mat'])%,'RegressorHat','Regressor','R','p');

dist=load(dryeeg_data_path(repoRoot, 'Rmaps/Figures/STNDistancetoSweetspot.mat'))%,'Distance','Regressor','subs')
ovlp=load(dryeeg_data_path(repoRoot, 'Rmaps/Figures/STNoverlap.mat'))%,'OverlapMotorSTN','Regressor','subs')

EReg=eeg.RegressorHat;
EReg([11,40],:)=[];
subid=[dist.subs;dist.subs];
statstb=table(EReg,dist.Distance,ovlp.OverlapMotorSTN,dist.Regressor,categorical(subid),'VariableNames',{'EP','Proximity','Overlap','UPDRS','SubID'});
statstbz=table(zscore(EReg),zscore(dist.Distance),zscore(ovlp.OverlapMotorSTN),zscore(dist.Regressor),categorical(subid),'VariableNames',{'EP','Proximity','Overlap','UPDRS','SubID'});

lmem=fitlme(statstbz,'UPDRS~EP+Proximity+Overlap+(1|SubID)');
diary([savepath,'linearmixedmodel_LOOCV.txt'])
lmem
disp(['Rsquared: ',num2str(lmem.Rsquared.Ordinary)]);
diary off


lm=fitlm(statstbz,'UPDRS~EP+Proximity+Overlap');
diary([savepath,'linearmodel_LOOCV.txt'])
lm
disp(['Rsquared: ',num2str(lm.Rsquared.Ordinary)]);
diary off

EReg=eeg10.RegressorHat;
EReg([11,40],:)=[];
subid=[dist.subs;dist.subs];
statstb=table(EReg,dist.Distance,ovlp.OverlapMotorSTN,dist.Regressor,categorical(subid),'VariableNames',{'EP','Proximity','Overlap','UPDRS','SubID'});
statstbz=table(zscore(EReg),zscore(dist.Distance),zscore(ovlp.OverlapMotorSTN),zscore(dist.Regressor),categorical(subid),'VariableNames',{'EP','Proximity','Overlap','UPDRS','SubID'});

lmem=fitlme(statstbz,'UPDRS~EP+Proximity+Overlap+(1|SubID)');
diary([savepath,'linearmixedmodel_10fold.txt'])
lmem
disp(['Rsquared: ',num2str(lmem.Rsquared.Ordinary)]);
diary off


lm=fitlm(statstbz,'UPDRS~EP+Proximity+Overlap');
diary([savepath,'linearmodel_10fold.txt'])
lm
disp(['Rsquared: ',num2str(lm.Rsquared.Ordinary)]);
diary off


EReg=eeg5.RegressorHat;
EReg([11,40],:)=[];
subid=[dist.subs;dist.subs];
statstb=table(EReg,dist.Distance,ovlp.OverlapMotorSTN,dist.Regressor,categorical(subid),'VariableNames',{'EP','Proximity','Overlap','UPDRS','SubID'});
statstbz=table(zscore(EReg),zscore(dist.Distance),zscore(ovlp.OverlapMotorSTN),zscore(dist.Regressor),categorical(subid),'VariableNames',{'EP','Proximity','Overlap','UPDRS','SubID'});

lmem=fitlme(statstbz,'UPDRS~EP+Proximity+Overlap+(1|SubID)');
diary([savepath,'linearmixedmodel_5fold.txt'])
lmem
disp(['Rsquared: ',num2str(lmem.Rsquared.Ordinary)]);
diary off


lm=fitlm(statstbz,'UPDRS~EP+Proximity+Overlap');
diary([savepath,'linearmodel_5fold.txt'])
lm
disp(['Rsquared: ',num2str(lm.Rsquared.Ordinary)]);
diary off




