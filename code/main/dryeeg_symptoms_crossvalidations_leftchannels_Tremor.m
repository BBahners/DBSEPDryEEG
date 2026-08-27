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
rmappath=dryeeg_result_path(repoRoot, 'symptoms_crossvalidations_Tremor');

% Make sure resultpath exists
if ~exist([rmappath])
    mkdir([rmappath]);
end
% Machine-specific addpath removed; see config/local_paths.m
load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat')); 
load(dryeeg_data_path(repoRoot, 'timevec.mat')); 
load(dryeeg_data_path(repoRoot, 'channels.mat')); 
corrtype='Spearman'; % use which type of correlation metric - for fMRI could also do Pearson

%% Load Regressor Table
tab=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp_Tremor.xlsx'));
tab([27],:)=[]; % P027 (left out for external validation)
allpts=1:size(tab,1); % all patients vector (cave: only used for first loop), then hemispheres (allconds)

%% Analysis time window:
ini=0.01; % beginning of time window in seconds
fin=0.2 ; % end of time window in seconds
win=[min(find(time>=ini)):min(find(time>=fin))];
Tmap=repmat(time(win),[32,1]); %time map
time=time(win); % time based on window defined

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
    Ri=right.m(:,win);
    Ri(lidx,:)=NaN;
    patConnectivityFiles{pt,1}=Ri(flip,:);
    patConnectivityFiles{size(allpts,2)+pt,1}=L;
    Regressor{pt,1}=tab.DeltaPercL(pt); % percentage improvement right hemiscores
    Regressor{size(allpts,2)+pt,1}=tab.DeltaPercR(pt); % percentage improvement left hemiscores
    Groupcell{pt}='1';

end

Regressor=cell2mat(Regressor);
channelvec=1:size(patConnectivityFiles{1},1);
allconds=1:size(Regressor,1);
%RegressorHat=NaN(size(allconds,2),1);
for pt=[0,allconds] % iterate through hemispheres, leaving out one each time. In first iteration ("0"), we will leave nothing out and generate an R-map over all hemispheres (usually denoted "R0map" in our nomenclature).

    otherpts=allconds;

    if pt>0 % in all but the first iteration, delete one hemisphere from "otherpts".
        otherpts(otherpts==pt)=[]; % delete hemisphere
    end

    %% 1. generate R-Maps
    % usage: Rmap=bb_Rmap_ep(fingerprints,regressor,output,corrtype)
    Rmap=bb_Rmap_ep(patConnectivityFiles(otherpts),Regressor(otherpts),...
        [rmappath,'R_',sprintf('%02.0f',pt),'.mat'],corrtype);

    %% 2. now compare that R-map with the left-out patient's connectivity

    if pt>0
        patConn=patConnectivityFiles{pt};
        thisPatConn=patConn(:);
        patRmap=load([rmappath,'R_',sprintf('%02.0f',pt),'.mat']);
        thisRmap=patRmap.Rmap(:);
        RegressorHat(pt)=corr(thisPatConn,thisRmap,'rows','pairwise','type',corrtype); % estimate of how similar this patient's connectivity is to the "optimal" connectivity profile denoted by the R-map (that is based on all patients except this particular one).
    end
end

group1.idx=[Groupcell,Groupcell];
group1.tag='';
%% 3. show correlation between similarities to "optimal" connectivity and empirical improvement
%usage: [h,R,p,g] = ea_corrplot(X,Y,permutation,labels,group1,group2,colors,markers,plottype,h)
[h,R,p,g]=ea_corrplot(Regressor(allconds)*100,RegressorHat(allconds)','no',{'Leave-One-Out Cross-Validation','% UPDRS Improvement','Similarity to R-Matrix'},group1,[],cmap(210,:),[],'linear');
h.Children(2).FontWeight='bold';
h.Children(2).TickDir="out";
h.Children(2).FontSize=14;

h.Children(3).Children.String{1}= [h.Children(3).Children.String{1}, ' (N =',sprintf('%d',size(Regressor,1)),')'];
h.Children(3).Children.String{2}=h.Children(3).Children.String{2}(12:end-1);
h.Children(3).Children.String(3)=[] ;

exportgraphics(h,[rmappath,'R_results_LOOCV_leftchannels.jpg'],'Resolution',600)
save([rmappath,'R_results_LOOCV_leftchannels.mat'],'RegressorHat','Regressor','R','p');
% set(gcf, 'Color', 'none');
% set(gca, 'Color', 'none');
% set(gca,'Box','off','Color','none')
% export_fig([savepath,'R_results_LOOCV_leftchannels_transparent.png'], '-png', '-transparent');


%saveas(h,[rmappath,'R_results_LOOCV_leftchannels.png']);
%saveas(h,[rmappath,'R_results_LOOCV_leftchannels.fig']);

save([rmappath,'Rmap.mat'],"Rmap")
save([rmappath,'Ch.mat'],"channlabels")
save([rmappath,'time.mat'],"time")
save([rmappath,'Chanord.mat'],"chanord")
save([rmappath,'channew.mat'],"chanlist_new")


%% 0. create a matrix of EEG Evoked Potential Files
for pt=allpts % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
    list=dir([cohmappath]);
    files=char(list.name);
    clear fn
    fn=files(contains(string(files),tab.SubID{pt}),:);
    left=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'le'),:),' ','')]);
    right=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'ri'),:),' ','')]);
    L=left.m(:,win);
    %L(ridx,:)=NaN;
    Ri=right.m(:,win);
    %Ri(lidx,:)=NaN;
    patConnectivityFiles{pt,1}=Ri(flip,:);
    patConnectivityFiles{size(allpts,2)+pt,1}=L;
    

end
Rmap=bb_Rmap_ep(patConnectivityFiles(allconds),Regressor(allconds),...
        [rmappath,'R_',sprintf('%02.0f',pt),'.mat'],corrtype);

% Time windows of interest
% Time window 1
ini=0.02;
fin=0.04;
win1=[min(find(time>=ini)):min(find(time>=fin))];
% Time window 2
ini=0.04;
fin=0.08;
win2=[min(find(time>=ini)):min(find(time>=fin))];
% Time window 3
ini=0.12;
fin=0.15;
win3=[min(find(time>=ini)):min(find(time>=fin))];

timwin={win1,win2,win3};
timwinlab={'20 to 40 ms', '50 to 80 ms','120 to 150 ms'};
intp=1000;
%% Topoplots

    for iwin=1:size(timwin,2)
%subplot(2,3,iwin)
h5=figure;
[h,ch_x,ch_y]=plot_topography(channlabels,(mean(Rmap(1:32,timwin{iwin})')),false(1),'10-20',false(1),false(1),intp);%,'false')%,1,1,0)
% hold on
% scatter(ch_x(lidx,:), ch_y(lidx,:), 60,'white', 'LineWidth',1);
ax=gcf;
%ax.Colormap=cmap;
cmp=ea_colorgradient(256,[1 1 1],[1 1 1],[0.3176 0.5216 0.2078]);
ax.Colormap=cmp;
title(timwinlab{iwin});
cb = colorbar();
        ylabel(cb,'Amplitude (ÂµV)','FontSize',10,'Rotation',90)
        clim([-0.4,0.4])
        %cb.Label.Position(1) = 3;
        h5.Color='white';
        ax=gca;
        ax.TickDir="out";

        saveas(h5,[rmappath,'topo_',timwinlab{iwin},'_green','.png']);
    end


