% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
%% Crossprediction analysis symptoms
foldcell={'LOOCV', '10-fold-CV','7-fold-CV','5-fold-CV','Split-half-CV'};
% Brad=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp_Brady.xlsx'));
% Rigid=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp_Rigidity.xlsx'));
% Trem=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp_Tremor.xlsx'));

load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat'));

rmappathbrady=dryeeg_result_path(repoRoot, 'symptoms_crossvalidations_Bradykinesia'); % results path
rmappathrigid=dryeeg_result_path(repoRoot, 'symptoms_crossvalidations_Rigidity'); % results path
rmappathtrem=dryeeg_result_path(repoRoot, 'symptoms_crossvalidations_Tremor');

Trem=load ([rmappathtrem,'R_results_',foldcell{1},'_leftchannels.mat']);

Brad=load ([rmappathbrady,'R_results_',foldcell{1},'_leftchannels.mat']);

Rigid=load ([rmappathrigid,'R_results_',foldcell{1},'_leftchannels.mat']);



h1=ea_corrplot(Trem.RegressorHat,Brad.Regressor,'no',{'Similarity to Tremor Map','Similarity to R-Matrix','Bradykinesia Improvement'},[],[],cmap(10,:))
h2=ea_corrplot(Trem.RegressorHat,Rigid.Regressor,'no',{'Similarity to Tremor Map','Similarity to R-Matrix','Rigidity Improvement'},[],[],cmap(10,:))
h3=ea_corrplot(Trem.RegressorHat,Trem.Regressor,'no',{'Similarity to Tremor Map','Similarity to R-Matrix','Tremor Improvement'},[],[],cmap(10,:))


h4=ea_corrplot(Brad.RegressorHat,Brad.Regressor,'no',{'Similarity to Bradykinesia Map','Similarity to R-Matrix','Bradykinesia Improvement'},[],[],cmap(30,:))
h5=ea_corrplot(Brad.RegressorHat,Rigid.Regressor,'no',{'Similarity to Bradykinesia Map','Similarity to R-Matrix','Rigidity Improvement'},[],[],cmap(30,:))
h6=ea_corrplot(Brad.RegressorHat,Trem.Regressor,'no',{'Similarity to Bradykinesia Map','Similarity to R-Matrix','Tremor Improvement'},[],[],cmap(30,:))

h7=ea_corrplot(Rigid.RegressorHat,Brad.Regressor,'no',{'Similarity to Rigidity Map','Similarity to R-Matrix','Bradykinesia Improvement'},[],[],cmap(50,:))
h8=ea_corrplot(Rigid.RegressorHat,Rigid.Regressor,'no',{'Similarity to Rigidity Map','Similarity to R-Matrix','Rigidity Improvement'},[],[],cmap(50,:))
h9=ea_corrplot(Rigid.RegressorHat,Trem.Regressor,'no',{'Similarity to Rigidity Map','Similarity to R-Matrix','Tremor Improvement'},[],[],cmap(50,:))


savepath=dryeeg_result_path(repoRoot, 'symptoms_crosssymptomprediction_leftchannels');
if~exist(savepath,'file')
mkdir(savepath)
end
saveas(h1,[savepath,'TremorMapBrady.png']);
saveas(h2,[savepath,'TremorMapRigid.png'])
saveas(h3,[savepath,'TremorMapTremo.png'])
saveas(h4,[savepath,'BradyMapBrady.png'])
saveas(h5,[savepath,'BradyMapRigid.png'])
saveas(h6,[savepath,'BradyMapTremo.png'])
saveas(h7,[savepath,'RigidMapBrady.png'])
saveas(h8,[savepath,'RigidMapRigid.png'])
saveas(h9,[savepath,'RigidMapTremo.png'])






