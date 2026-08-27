%% Dry EEG DBS Network Mapping Analysis
%  Based on LeadDBS LOOCV code
%  Author: Bahne H. Bahners

clear variables
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
%close all
% Machine-specific addpath removed; see config/local_paths.m
% Machine-specific addpath removed; see config/local_paths.m
cohmappath=[dryeeg_data_path(repoRoot, 'EPmaps_Discovery'), filesep]; % connectivity map path
rmappath=dryeeg_result_path(repoRoot, 'primary_circular_leftchannels');

% Machine-specific addpath removed; see config/local_paths.m
load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat')); 
load(dryeeg_data_path(repoRoot, 'timevec.mat')); 
load(dryeeg_data_path(repoRoot, 'channels.mat')); 
tab=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp.xlsx'));
corrtype='Spearman'; % use which type of correlation metric - for fMRI could also do Pearson
mirror=false;
% Make sure resultpath exists
if ~exist([rmappath])
    mkdir([rmappath]);
end

%% Analysis time window:
ini=0.01; % beginning of time window in seconds
fin=0.2; % end of time window in seconds
win=[min(find(time>=ini)):min(find(time>=fin))];
Tmap=repmat(time(win),[32,1]); %time map
time=time(win); % time based on window defined

%% Load Clinical Data
% left and right scores & improvements are defined in Table

tab([27],:)=[]; % Exclude P014, P019, P022 (Line noise artefacts) & P027 (left out for external validation)
%tab(isnan(tab.DeltaPerc),:)=[]; % Exclude Patients with missing subitems
allpts=1:size(tab,1); % all patients vector (cave: only used for first loop), then hemispheres (allconds)

%% Prepare Channel Order & Channel Flip
channels(:,33:35)=[];
for ii=1:size(channels,2)
    channlabels{ii}=channels(ii).Name; % create cell array with EEG channel labels
end
chanord= bb_chanord(channlabels);% now get a new channelorder for the plots (Frontal to Occipital)
chanlist_new=char(channlabels{chanord}); % now reorder the channellist
% new channellabels are assigned to essentially flip the EEG data
chanside={'right','right','right','right','right','right','right','right','right','right','left','left','z','right','z','right','left','left','z','right','right','z','left','left','left','left','left','left','left','left','left','left'};
chanflip={'P7',	'T7', 'CP5',	'FC5',	'F7',	'F3',	'C3',	'P3',	'AF3',	'Fp1',	'Fp2',	'AF4',	'Fz',	'FC1',	'Cz',	'CP1',	'PO4',	'O2',   'Oz',	'O1',	'PO3',	'Pz',	'CP2',	'FC2',	'P4',	'C4',	'F4',	'F8',	'FC6',	'CP6',	'T8',	'P8'};
[~, flip] = ismember(chanflip,channlabels); % Find the corresponding indices for channel reordering
ridx=contains(chanside,'right');
lidx=contains(chanside,'left');

%% 0. create a matrix of EEG Evoked Potential Files
for pt=allpts % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
    list=dir([cohmappath]);
    files=char(list.name);
    clear fn
    fn=files(contains(string(files),tab.SubID{pt}),:);
    left=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'le'),:),' ','')]);
    right=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'ri'),:),' ','')]);
    L=left.m(:,win);
    L(ridx,:)=NaN;
    R=right.m(:,win);
    R(lidx,:)=NaN;
    patConnectivityFiles{pt,1}=R(flip,:);
    patConnectivityFiles{size(allpts,2)+pt,1}=L;
    patConnectivitySubject{pt,1}=tab.SubID{pt};
    patConnectivitySubject{size(allpts,2)+pt,1}=tab.SubID{pt};
    %% CAVE: improvement for other body half !
    Regressor{pt,1}=tab.DeltaPercL(pt); % percentage improvement right hemiscores
    Regressor{size(allpts,2)+pt,1}=tab.DeltaPercR(pt); % percentage improvement left hemiscores
    %Regressor=randn(length(allpts),1); % random regressor
    Groupcell{pt}='1';
end

Regressor=cell2mat(Regressor);
channelvec=1:size(patConnectivityFiles{1},1);
allconds=1:size(Regressor,1);
[d,e,subidx]=unique(patConnectivitySubject);

%% 1. generate R-Maps
% usage: Rmap=bb_Rmap_ep(fingerprints,regressor,output,corrtype)
Rmap=bb_Rmap_ep(patConnectivityFiles,Regressor,...
    [rmappath,'R_',sprintf('%02.0f',0),'.mat'],corrtype);
Chmap=repmat(channelvec,[size(patConnectivityFiles{1},2),1])';
h2=figure;
imagesc(time,channelvec,Rmap(chanord,:))
colormap(cmap)
clim([-0.4,0.4]);
xlabel('Time (s)','FontSize',16)
yticklabels(chanlist_new)
yticks(channelvec)
cb = colorbar();
ylabel(cb,'R','FontSize',16,'Rotation',0)
cb.Label.Position(1) = 3;
h2.Color='white';
ax=gca;
ax.TickDir="out";
saveas(h2,[rmappath,'R_',sprintf('%02.0f',0),'.fig']);
saveas(h2,[rmappath,'R_',sprintf('%02.0f',0),'.png']);

patRmap=load([rmappath,'R_',sprintf('%02.0f',0),'.mat']);
thisRmap=patRmap.Rmap(:);

for pt=allconds % iterate through hemispheres, leaving out one each time. In first iteration ("0"), we will leave nothing out and generate an R-map over all hemispheres (usually denoted "R0map" in our nomenclature).
    %% 2. now compare that R-map with the left-out patient's connectivity
    patConn=patConnectivityFiles{pt};
    thisPatConn=patConn(:);
    RegressorHat(pt)=corr(thisPatConn,thisRmap,'rows','pairwise','type',corrtype); % estimate of how similar this patient's connectivity is to the "optimal" connectivity profile denoted by the R-map (that is based on all patients except this particular one).
end

group1.idx=[Groupcell,Groupcell];
group1.tag='';
%% 3. show correlation between similarities to "optimal" connectivity and empirical improvement
%usage: [h,R,p,g] = ea_corrplot(X,Y,permutation,labels,group1,group2,colors,markers,plottype,h)
[h,R,p,g]=ea_corrplot(Regressor(allconds)*100,RegressorHat(allconds)','no',{'Similarity to R-Matrix','% UPDRS Improvement','Similarity to R-Matrix'},group1,[],cmap(210,:),[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Map'}

h.Children(2).FontWeight='bold';
h.Children(2).TickDir="out";
h.Children(2).FontSize=14;
h.Position(4)=h.Position(3)+h.Position(3)*0.025;

%h.Position=[100 100 1000 1000];
exportgraphics(h,[rmappath,'R_results_leftchannels.jpg'],'Resolution',300)
%saveas(h,[rmappath,'R_results_leftchannels.png']);
%saveas(h,[rmappath,'R_results_leftchannels.fig']);

save([rmappath,filesep,'stats_results_leftchannels.mat'],'R','p','RegressorHat','Regressor');

%% 4. Create Average Channels x Time Map

Amap=mean(cat(3,patConnectivityFiles{:}),3);
col=cmap;
h4=figure;
imagesc(time,channelvec,Amap(chanord,:)*10e6)
colormap(col)
clim([-15,15]);
xlabel('Time (s)','FontSize',16)
yticklabels(chanlist_new)
yticks(channelvec)
cb = colorbar();
ylabel(cb,'Amplitude (ÂµV)','FontSize',16,'Rotation',90)
cb.Label.Position(1) = 3;
h4.Color='white';
ax=gca;
ax.TickDir="out";

saveas(h4,[rmappath,'Amap_',sprintf('%02.0f',pt),'.fig']);
saveas(h4,[rmappath,'Amap_',sprintf('%02.0f',pt),'.png']);

save([rmappath,'Amap.mat'],"Amap")
save([rmappath,'Rmap.mat'],"Rmap")
save([rmappath,'Ch.mat'],"channlabels")
save([rmappath,'time.mat'],"time")
save([rmappath,'Chanord.mat'],"chanord")
save([rmappath,'channew.mat'],"chanlist_new")

% %% 0. create a matrix of EEG Evoked Potential Files
% for pt=allpts % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
%     list=dir([cohmappath]);
%     files=char(list.name);
%     clear fn
%     fn=files(contains(string(files),tab.SubID{pt}),:);
%     left=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'le'),:),' ','')]);
%     right=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'ri'),:),' ','')]);
%     L=left.m(:,win);
%     %L(ridx,:)=NaN;
%     R=right.m(:,win);
%    % R(lidx,:)=NaN;
%     patConnectivityFiles{pt,1}=R(flip,:);
%     patConnectivityFiles{size(allpts,2)+pt,1}=L;
%     patConnectivitySubject{pt,1}=tab.SubID{pt};
%     patConnectivitySubject{size(allpts,2)+pt,1}=tab.SubID{pt};
%     %% CAVE: improvement for other body half !
%     Regressor{pt,1}=tab.DeltaPercL(pt); % percentage improvement right hemiscores
%     Regressor{size(allpts,2)+pt,1}=tab.DeltaPercR(pt); % percentage improvement left hemiscores
%     %Regressor=randn(length(allpts),1); % random regressor
%     Groupcell{pt}='1';
% end


