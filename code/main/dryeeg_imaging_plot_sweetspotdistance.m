%% Distance to PD Sweetspot Dry EEG Cohort
% Authors: Bahne H Bahners, Andreas Horn

clear variables
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
 % Machine-specific addpath removed; see config/local_paths.m
% Machine-specific addpath removed; see config/local_paths.m
load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat')); 
groupfile=dryeeg_data_path(repoRoot, 'dataset-DryEEGLeadDBS240715_analysis-20240724175347.mat');

rmappath=dryeeg_result_path(repoRoot, 'imaging_plot_sweetspotdistance');
if ~exist(rmappath,"file")
    mkdir(rmappath)
end
tab=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp.xlsx'));

load (groupfile);

mnisweet={[12.58 -13.41 -5.87],[-12.58 -13.41 -5.87]};

allpts=1:length(M.patient.list);
for pt=allpts
    for side=1:2
        activectx{side}(pt,:)=mean(M.elstruct(pt).coords_mm{side}(...
            logical(M.S(pt).activecontacts{side}),:),1);
    end
end

for pt=allpts
    for side=1:2

        Dist{side}(pt,:)=norm(mnisweet{side}-activectx{side}(pt,:))*-1;

    end
end


%% Exclude Subject 27
tab([11,27],:)=[];

allpts=1:size(tab,1);

Distance=[Dist{1};Dist{2}];
Regressor=[tab.DeltaPercL;tab.DeltaPercR];
for i=1:height(tab)*2

    Groupcell{i}='1';

end

group1.idx=Groupcell;
group1.tag='';
[h,R,p,g]=ea_corrplot(Regressor*100,Distance,'no',...
    {['Distance to STN Sweetspot'],'% UPDRS Improvement','Distance to Sweetspot (mm)'},...
    group1,[],cmap(10,:),[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Map'}

h.Children(2).FontWeight='bold';
h.Children(2).TickDir="out";
h.Children(2).FontSize=14;

h.Position=[100 100 550 580];
subs=tab.SubID;
saveas(h,[rmappath,'/STNDistancetoSweetspot.png'])
save([rmappath,'/STNDistancetoSweetspot.mat'],'Distance','Regressor','subs')

%h.Position=[100 100 1000 1000];
%exportgraphics(h,[rmappath,'R_results_DistanceSweetspot.jpg'],'Resolution',300)




