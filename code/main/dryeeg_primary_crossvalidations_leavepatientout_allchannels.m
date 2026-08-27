%% Dry EEG DBS Network Mapping Analysis
%  Based on LeadDBS Kfold code
%  Author: Bahne H. Bahners, Lukas L. Goede
clear variables
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
close all

%% Paths and Metadata
% Machine-specific addpath removed; see config/local_paths.m
% Machine-specific addpath removed; see config/local_paths.m
% Machine-specific addpath removed; see config/local_paths.m
cohmappath=[dryeeg_data_path(repoRoot, 'EPmaps_Discovery'), filesep]; % connectivity map path
rmappath=dryeeg_result_path(repoRoot, 'primary_crossvalidations_allchannels'); % results path
if ~exist([rmappath])
    mkdir([rmappath]);
end
% Machine-specific addpath removed; see config/local_paths.m
load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat')); 
load(dryeeg_data_path(repoRoot, 'timevec.mat')); 
load(dryeeg_data_path(repoRoot, 'channels.mat')); 
corrtype='Spearman'; % use which type of correlation metric - for fMRI could also do Pearson

% Load Improvement Data
tab=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp.xlsx'));
tab([27],:)=[]; % P027 used for prospective validation
allpts=1:size(tab,1); % all patients vector (cave: only used for first loop), then hemispheres (allconds)

%% Analysis time window:
ini=0.01; % beginning of time window in seconds
fin=0.2 ;% end of time window in seconds
win=[min(find(time>=ini)):min(find(time>=fin))];
Tmap=repmat(time(win),[32,1]); %time map
time=time(win); % time based on window defined

%% Prepare Channel Order & Channel Flip
channels(:,33:35)=[];
for ii=1:size(channels,2)
    channlabels{ii}=channels(ii).Name; % create cell array with EEG channel labels
end
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
    %L(ridx,:)=NaN;
    R=right.m(:,win);
    %R(lidx,:)=NaN;
    % Connectivity Files:
    patConnectivityFiles{pt,1}=R(flip,:);
    patConnectivityFiles{size(allpts,2)+pt,1}=L;
    % Regressor
    Regressor(pt,1)=tab.DeltaPercL(pt);
    Regressor(size(allpts,2)+pt,1)=tab.DeltaPercR(pt); % percentage improvement left hemiscores
    Groupcell{pt}='1';
end
folds=[29];%[58,10, 7, 5, 2];
foldcell={'LOPO-CV'};%{'LOO-CV', '10-fold-CV','7-fold-CV','5-fold-CV','Split-half-CV'};

  
    RegressorHat=nan(length(Regressor),1);
    
 for pt=allpts% iterate through hemispheres, leaving out one each time. In first iteration ("0"), we will leave nothing out and generate an R-map over all hemispheres (usually denoted "R0map" in our nomenclature).

          allconds=1:max(allpts)*2;
         
          idx_test=allconds(allconds==pt|allconds==pt+max(allpts))
          idx_training=allconds(allconds~=pt&allconds~=pt+max(allpts));
        
          %% 1. generate R-Map for this fold
        Rmap=bb_Rmap_ep(patConnectivityFiles(idx_training),Regressor(idx_training),[rmappath,'Rmap_train',num2str(pt),'.mat'],corrtype);
        %% 2. now compare that R-map with the patient connectivity
        testsubs=idx_test;
        for pts=testsubs
            patConn=patConnectivityFiles{pts};
            %RegressorHattmp(pt)=corr(patConn(:),Rmap(:),'rows','pairwise','type',corrtype); % estimate of how similar this patient's connectivity is to the "optimal" connectivity profile denoted by the R-map (that is based on all patients except this particular one).
            RegressorHat(pts)=corr(patConn(:),Rmap(:),'rows','pairwise','type',corrtype); % estimate of how similar this patient's connectivity is to the "optimal" connectivity profile denoted by the R-map (that is based on all patients except this particular one).

        end
    end
    %RegressorHatout(:,iter)=RegressorHattmp;
    %end
    %MRegressorHat=mean(RegressorHatout,2);

   for pt=allpts% iterate through hemispheres, leaving out one each time. In first iteration ("0"), we will leave nothing out and generate an R-map over all hemispheres (usually denoted "R0map" in our nomenclature).

          allconds=1:max(allpts)*2;
         
          idx_test=allconds(allconds==pt|allconds==pt+max(allpts)); 
          idx_training=allconds(allconds~=pt&allconds~=pt+max(allpts));
          testsubs=idx_test;

        for pts=testsubs
            mdl{pts}=fitlm(RegressorHat(idx_training),Regressor(idx_training)*100);
            Pred(pts)=predict(mdl{pts},RegressorHat(pts));
            MAEpts(pts) = mean(abs(Pred(pts)-Regressor(pts)));
            RMSpts(pts)=mdl{pts}.RMSE;
        end
    end

    group1.idx=[Groupcell,Groupcell];
    group1.tag='';

    [h,R,p,g]=ea_corrplot(Regressor*100,Pred','no',{foldcell{1},'% UPDRS Improvement (Empirical)','% UPDRS Improvement (Estimated)'},group1,[],cmap(210,:),[],'linear');
    RMS=mean(RMSpts);
    MAE=mean(MAEpts);
    h.Children(2).FontWeight="bold";
    h.Children(2).TickDir= "out";
    h.Children(2).FontSize=14;
    h.Children(3).Children.String{3}=[];%['R^{2} = ',sprintf('%.2f',windraw.Rsquared.Ordinary)];
    h.Children(3).Children.Interpreter='tex';
    if  p.spearman<0.001
        h.Children(3).Children.String{2}= ['R = ', sprintf('%.2f',R.spearman),', p = ', sprintf('%.2e',p.spearman)];
        h.Children(3).Children.String{3}= ['R^{2} = ','; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];

    else
        h.Children(3).Children.String{2}= ['R = ', sprintf('%.2f',R.spearman),', p = ', sprintf('%.3f',p.spearman)];
        h.Children(3).Children.String{3}= ['R^{2} = ','; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];

    end
    h.Position(4)=h.Position(3)+h.Position(3)*0.025;

    saveas(h,[rmappath,'Results_leftchannels_pred_',foldcell{1},'.png']);
    save([rmappath,'Results_leftchannels_pred_',foldcell{1},'.mat'],'RegressorHat','Regressor','MAE','RMS','R','p');


    group1.idx=[Groupcell,Groupcell];
    group1.tag='';
    %% 3. show correlation between similarities to "optimal" connectivity and empirical improvement
    %usage: [h,R,p,g] = ea_corrplot(X,Y,permutation,labels,group1,group2,colors,markers,plottype,h)
    %[h,R,p,g]=ea_corrplot(Regressor*100,MRegressorHat','no',{'10-fold Cross-Validation (Average across 1000 iterations) Left Channels','% UPDRS Improvement','Similarity to R-Matrix'},group1,[],cmap(25,:),[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Map'}
    [h,R,p,g]=ea_corrplot(Regressor*100,RegressorHat','no',{foldcell{1},'% UPDRS Improvement','Similarity to R-Matrix'},group1,[],cmap(25,:),[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Map'}
    h.Position(4)=h.Position(3)+h.Position(3)*0.025;
 h.Children(2).FontWeight="bold";
    h.Children(2).TickDir= "out";
    h.Children(2).FontSize=14;
    saveas(h,[rmappath,'Results_leftchannels_',foldcell{1},'.png']);
    save([rmappath,'Results_leftchannels_',foldcell{1},'.mat'],'RegressorHat','Regressor','R','p');



