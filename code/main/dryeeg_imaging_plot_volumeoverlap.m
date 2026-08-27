% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
%% Plot Volume OVERLAP
% Author: Bahne Bahners, BWH
tab=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp.xlsx'));

rmappath=dryeeg_result_path(repoRoot, 'imaging_plot_volumeoverlap');
if ~exist(rmappath,"file")
    mkdir(rmappath)
end
% Historical personal path removed for public repository.
tab2=readtable(dryeeg_data_path(repoRoot, 'Native_Volume_Overlap.csv'));
tab([11,27],:)=[];
%tab2([find(contains(tab2.ID,'27')),find(contains(tab2.ID,'001'))],:)=[];
for pt=1:height(tab)
    ord(pt)=find(contains(tab2.ID,tab(pt,:).SubID));
end
tab2=tab2(ord,:);
Regressor=[tab.DeltaPercL;tab.DeltaPercR];
OverlapMotorSTN=[tab2.STN_RH_Volume_Overlap_MM;tab2.STN_LH_Volume_Overlap_MM];
OverlapMotorSweetSpot=[tab2.TOR_PSM_Sweetspot_RH_Volume_Overlap_MM;tab2.TOR_PSM_Sweetspot_LH_Volume_Overlap_MM];

for i=1:height(tab)*2
    Groupcell{i}='1';
end

group1.idx=Groupcell;
group1.tag='';

h=ea_corrplot(Regressor*100,OverlapMotorSTN,'noperm',{'VTA Overlap with STN','%UPDRS Improvement','Volume Overlap (mm)'},...
    group1,[],cmap(10,:));
h.Children(2).FontSize=14;
h.Children(2).FontWeight='bold';
h.Children(2).TickDir="out";
h.Position=[100 100 550 580];

saveas(h,[rmappath,'/STNoverlap.png'])
subs=tab.SubID;

save([rmappath,'/STNoverlap.mat'],'OverlapMotorSTN','Regressor','subs')


h=ea_corrplot(Regressor*100,OverlapMotorSweetSpot,'noperm',{'VTA Overlap with STN Sweetspot','%UPDRS Improvement','Volume Overlap (mm)'},...
    group1,[],cmap(10,:));
h.Children(2).FontSize=14;
h.Children(2).FontWeight='bold';
h.Children(2).TickDir="out";
h.Position=[100 100 550 580];

saveas(h,[rmappath,'/SweetSpotOverlap.png'])




