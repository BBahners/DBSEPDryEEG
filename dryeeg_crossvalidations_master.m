%% Dry EEG DBS Network Mapping Analysis Master Script Primary Analysis
%  Based on LeadDBS Kfold code
%  Author: Bahne H. Bahners, Lukas L. Goede
clear variables
close all

%% Paths and Metadata
addpath ./spm12/
addpath(genpath('./leaddbs'));
addpath /functions/ 

load /data/timevec.mat
load /data/channels.mat
corrtype='Spearman'; % use which type of correlation metric - for fMRI could also do Pearson

% Load Improvement Data
SubID={...
    'P001'
    'P002'
    'P003'
    'P004'
    'P005'
    'P006'
    'P007'
    'P008'
    'P009'
    'P010'
    'P011'
    'P012'
    'P013'
    'P014'
    'P015'
    'P016'
    'P017'
    'P018'
    'P019'
    'P020'
    'P021'
    'P022'
    'P023'
    'P024'
    'P025'
    'P026'
    'P027'
    'P028'
    'P029'};
DeltaPercR=randn(29,1); % right hemibody improvement (random placeholder)
DeltaPercL=randn(29,1); % left hemibody improvement (random placeholder)
tab=table(SubID,DeltaPercR,DeltaPercL);
allpts=1:size(tab,1); % all patients vector (cave: only used for first loop), then hemispheres (allconds)

%% Analysis time window:
ini=0.01; % beginning of time window in seconds
fin=0.2 ;% end of time window in seconds
win=[min(find(time>=ini)):min(find(time>=fin))];
Tmap=repmat(time(win),[32,1]); %time map
time=time(win); % time based on window defined

%% Prepare Channel Order & Channel Flip
channels(:,33:35)=[]; % Accelerometer channels deleted
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
    left.m=randn(32,251); % random placeholder
    right.m=randn(32,251); % random placeholder
    L=left.m(:,win);
    L(ridx,:)=NaN;
    R=right.m(:,win);
    R(lidx,:)=NaN;
    % Connectivity Files:
    patConnectivityFiles{pt,1}=R(flip,:);
    patConnectivityFiles{size(allpts,2)+pt,1}=L;
    % Regressor
    Regressor(pt,1)=tab.DeltaPercL(pt);
    Regressor(size(allpts,2)+pt,1)=tab.DeltaPercR(pt); % percentage improvement right hemiscores for left sided stimulation ('L')
    Groupcell{pt}='1';
end

%% Cross-Validation

folds=[58,10,5];
foldcell={'LOO-CV', '10-fold-CV','5-fold-CV'};

    for ifolds=1:length(folds)

    rngseed = 'default';
    rng(rngseed);

    c=cvpartition(length(Regressor),'Kfold',folds(ifolds)); % can define different size of folds
    RegressorHat=nan(length(Regressor),1);

        for i=1:c.NumTestSets
        idx_test=test(c,i); % defining test indices of improvement score
        idx_training=training(c,i); % define the training set (will be used to train ur model as in leave one out u were traning on eg. 58 improvement our of 59

        %% 1. generate R-Map for this fold
        Rmap=bb_Rmap_ep(patConnectivityFiles(idx_training),Regressor(idx_training),[],corrtype,'dontsave');
       
        if i==1 && ifolds==1 % in first iteration plot exemplary R-Map 
            %% 1.1 Plot R-Map
            %Preparing the plot /channel order etc.
            chanord= bb_chanord(channlabels);% now get a new channelorder for the plots (Frontal to Occipital)
            chanlist_new=char(channlabels{chanord}); % now reorder the channellist
            halfidx=contains(chanside,'left')|contains(chanside,'z') ; % only plot left sided channels or central ones
            %Rmap(~halfidx,:)=NaN;
            Rmapplot=Rmap(chanord,:);
            chanlist_new(sum(isnan(Rmapplot),2)>1,:)=[];
            channelvechalf=1:length(chanlist_new);
            Rmapplot(sum(isnan(Rmapplot),2)>1,:)=[];
            Chmap=repmat(channelvechalf,[size(patConnectivityFiles{1},2),1])';
            % Plot image here
            figure;
            imagesc(time,channelvechalf,Rmapplot)
            clim([min(min(Rmapplot)),max(max(Rmapplot))]);
            xlabel('Time (s)','FontSize',16)
            yticklabels(chanlist_new)
            yticks(channelvechalf)
            cb = colorbar();
            ylabel(cb,'R','FontSize',10,'Rotation',0)

        end

            %% 2. now compare that R-map with the patient connectivity
            testsubs=find(idx_test)';
            for pt=testsubs
                patConn=patConnectivityFiles{pt};
                RegressorHat(pt)=corr(patConn(:),Rmap(:),'rows','pairwise','type',corrtype); % estimate of how similar this patient's connectivity is to the "optimal" connectivity profile denoted by the R-map (that is based on all patients except this particular one).
            end

        end


        %% 3. show correlation between similarities to "optimal" connectivity and empirical improvement
        % usage: [h,R,p,g] = ea_corrplot(X,Y,permutation,labels,group1,group2,colors,markers,plottype,h)
        group1.idx=[Groupcell,Groupcell];
        group1.tag='';
        [h,R,p,g]=ea_corrplot(Regressor*100,RegressorHat','no',{foldcell{ifolds},'% UPDRS Improvement','Similarity to R-Matrix'},group1,[],[],[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Map'}
        h.Position(4)=h.Position(3)+h.Position(3)*0.025;
        h.Children(2).FontWeight="bold";
        h.Children(2).TickDir= "out";
        h.Children(2).FontSize=14;

    end
