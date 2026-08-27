%% LOAD CONTACT TEST DATA AND PLOT FIGURES
% Author: Bahne H. Bahners, MD
clear variables
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
close all
% Machine-specific addpath removed; see config/local_paths.m
% Machine-specific addpath removed; see config/local_paths.m
rmappath=dryeeg_result_path(repoRoot, 'primary_circular_leftchannels'); % results path
savepath=dryeeg_result_path(repoRoot, 'validation_alltrials_lme');

% Machine-specific addpath removed; see config/local_paths.m
if ~exist(savepath,'file')
    mkdir(savepath);
end
load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat')); 
load(dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/channels.mat')); 
corrtype='Spearman'; % use which type of correlation metric - for fMRI could also do Pearson
load([rmappath,'R_00.mat']); % load Rmap
rmaptime=load([rmappath,'time.mat']); % load time definition for Rmap
fptime=load(dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/timevec.mat'));
tbs=readtable(dryeeg_data_path(repoRoot, 'metadata/DryEEGOutsample_nonans_new.xlsx'));
files={...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S027_left_contact8_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S027_left_contact9_data_stim.mat')	,...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S027_left_contact10_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S027_left_contact11_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S027_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S027_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S027_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S027_right_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S031_left_contact8_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S031_left_contact9_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S031_left_contact10_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S031_left_contact11_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S032_left_contact8_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S032_left_contact9_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S032_left_contact10_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S032_left_contact11_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S032_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S032_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S032_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S032_right_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S033_left_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S033_left_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S033_left_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S033_left_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S033_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S033_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S033_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S033_right_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S034_left_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S034_left_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S034_left_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S034_left_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S034_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S034_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S034_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S034_right_contact3_data_stim.mat'),...
     dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S035_left_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S035_left_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S035_left_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S035_left_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S035_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S035_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S035_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S035_right_contact3_data_stim.mat'),...
     dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S036_left_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S036_left_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S036_left_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S036_left_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S036_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S036_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S036_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S036_right_contact3_data_stim.mat'),...
     dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S037_left_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S037_left_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S037_left_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S037_left_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S037_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S037_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S037_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S037_right_contact3_data_stim.mat'),...
     dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S038_left_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S038_left_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S038_left_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S038_left_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S038_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S038_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S038_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S038_right_contact3_data_stim.mat'),...
     dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S039_left_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S039_left_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S039_left_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S039_left_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S039_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S039_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S039_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S039_right_contact3_data_stim.mat'),...
     dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S040_left_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S040_left_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S040_left_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S040_left_contact3_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S040_right_contact0_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S040_right_contact1_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S040_right_contact2_data_stim.mat'),...
    dryeeg_data_path(repoRoot, 'EPmapsDUSnoH/S040_right_contact3_data_stim.mat'),...
    };



effectthreshold=tbs.effectthreshold;
sidethreshold=tbs.sidethreshold;
effectlog=tbs.effectlog;



%% Prepare Channel Order & Channel Flip
%channels(:,33:35)=[];
for ii=1:size(channels,2)
    channlabels{ii}=channels(ii).Name; % create cell array with EEG channel labels
end
chanord= bb_chanord(channlabels);% now get a new channelorder for the plots (Frontal to Occipital)
chanlist_new=char(channlabels{chanord}); % now reorder the channellist

%% Channelflip
% new channellabels are assigned to essentially flip the EEG data (see
% original channellabels below in commented line and below the
% corresponding channels for the other hemisphere)
%'P8'	'T8'  'CP6'	    'FC6'	'F8'	'F4'	'C4'	'P4'	'AF4'	'Fp2'	'Fp1'	'AF3'	'Fz'	'FC2'	'Cz'	'CP2'	'PO3'	'O1'	'Oz'	'O2'	'PO4'	'Pz'	'CP1'	'FC1'	'P3'	'C3'	'F3'	'F7'	'FC5'	'CP5'	'T7'	'P7'
chanside={'right','right','right','right','right','right','right','right','right','right','left','left','z','right','z','right','left','left','z','right','right','z','left','left','left','left','left','left','left','left','left','left'};
chanflip={'P7',	'T7', 'CP5',	'FC5',	'F7',	'F3',	'C3',	'P3',	'AF3',	'Fp1',	'Fp2',	'AF4',	'Fz',	'FC1',	'Cz',	'CP1',	'PO4',	'O2',   'Oz',	'O1',	'PO3',	'Pz',	'CP2',	'FC2',	'P4',	'C4',	'F4',	'F8',	'FC6',	'CP6',	'T8',	'P8'};
[~, flip] = ismember(chanflip,channlabels); % Find the corresponding indices for channel reordering
ridx=contains(chanside,'right');
lidx=contains(chanside,'left');
%          chanflip={'P7',	'T7', 'CP5',	'FC5',	'F7',	'F3',	'C3',	'P3',	'AF3',	'Fp1',	'Fp2',	'AF4',	'Fz',	'FC1',	'Cz',	'CP1',	'PO4',	'O2',   'Oz',	'O1',	'PO3',	'Pz',	'CP2',	'FC2',	'P4',	'C4',	'F4',	'F8',	'FC6',	'CP6',	'T8',	'P8'};
% [~, flip] = ismember(chanflip,channlabels); % Find the corresponding indices for channel reordering

%% Plot R-map
figure;
Chmap=repmat(1:32,[32,1])';
imagesc(rmaptime.time,1:32,Rmap(chanord,:))
colormap(cmap)
clim([-0.4,0.4]);

%% 0. create a matrix of connectivity files
for pt=1:size(files,2) % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
    tmp=load(files{pt});
    if contains(files{pt},'left')
        L=tmp.m(1:32,:);
        L(ridx,:)=NaN;
        patConnectivityFiles{pt,1}=L;
        patConnectivityFileName{pt,1}=files{pt};
    elseif contains(files{pt},'right')
        Ri=tmp.m(1:32,:);
        Ri(lidx,:)=NaN;
        patConnectivityFiles{pt,1}=Ri(flip,:);
        patConnectivityFileName{pt,1}=files{pt};
    end
end
% for pt=1:size(files,2) % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
% 
%         patConnectivityFiles{pt,1}(find(~isnan(patConnectivityFiles{pt,1}(:,1))),:)=zscore(patConnectivityFiles{pt,1}(find(~isnan(patConnectivityFiles{pt,1}(:,1))),:));
% 
% end
%% Align Time Vectors
% Identify subset of Patient Channel X Amplitude map (fptime) in the time domain
% that fits the Rmap time definition (rmaptime)
TidxPatConn=ismember(round(fptime.time,5),round(rmaptime.time,5)); %time idices of fingerprint map matching rmap
TidxRmap=ismember(round(rmaptime.time,5),round(fptime.time,5)); % time indices of Rmap matching fingerprint map

%% Time Vectors
timevectorscomp=[fptime.time(TidxPatConn);rmaptime.time(TidxRmap)]';
Rmap=Rmap(:,TidxRmap(1:end)); % select subset of Rmap matching fingerprint map

allconds=1:size(files,2);
for pt=allconds
    patConn=patConnectivityFiles{pt};
    subplot(size(allconds,2)/2,size(allconds,2)/2,pt)
    plot(fptime.time(TidxPatConn),patConn(:,TidxPatConn)')
    title(patConnectivityFileName{pt})
    patConn=patConn(:,TidxPatConn);
    thisPatConn=patConn(:);
    subplot(size(allconds,2)/2,size(allconds,2)/2,pt+1)
    plot(rmaptime.time(TidxRmap),Rmap')
    thisRmap=Rmap(:);
    RegressorHat{pt,1}=corr(thisPatConn,thisRmap,'rows','pairwise','type',corrtype); % estimate of how similar this patient's connectivity is to the "optimal" connectivity profile denoted by the R-map (that is based on all patients except this particular one).
    [a,b]=fileparts(files{pt});
    RegressorHat{pt,2}=regexprep(b,'_',' ');
end

% figure;
% for pt=allconds
%     patConn=patConnectivityFiles{pt};
%     subplot(size(allconds,2)/2,size(allconds,2)/2,pt)
%     imagesc(fptime.time(TidxPatConn),1:32,patConn(chanord,TidxPatConn)*10e6')
%     clim([-10 10])
%     colormap(cmap)
%     title(patConnectivityFileName{pt}(end-15:end))
%     patConn=patConn(:,TidxPatConn);
%     thisPatConn=patConn(:);
%     subplot(size(allconds,2)/2,size(allconds,2)/2,pt+1)
%     imagesc(rmaptime.time(TidxRmap),1:32,Rmap(chanord,:))
% end

tab=table(string(char(RegressorHat{:,2})),vertcat(RegressorHat{:,1}),effectthreshold,sidethreshold,(((sidethreshold-effectthreshold)./(effectthreshold))*100),'VariableNames',{'Contact','SpatialCorrelation','ClinicalThreshold','SideEffectThreshold','TherapeuticWindow'});
%tab=table(string(char(RegressorHat{:,2})),vertcat(RegressorHat{:,1}),effectthreshold',sidethreshold',(sidethreshold-effectthreshold)','VariableNames',{'Contact','SpatialCorrelation','ClinicalThreshold','SideEffectThreshold','TherapeuticWindow'})
writetable(tab,[savepath,'R_results_SpatialCorrs.xlsx']);

%Groupcell4={'1','1','1','1','2','2','2','2','3','3','3','3','4','4','4','4','5','5','5','5',};
for hh=1:height(tbs)
Groupcell4{hh}=char([tbs.SubID{hh},tbs.side{hh}]);
end
%Groupcell4=tbs.SubID;

group4.idx=Groupcell4';
group4.tag='';

% Groupcell5=deltaperc;
% group5.idx=Groupcell5';
% group5.tag='';

tab.Hemi=Groupcell4';
%tab.Hemi=Groupcell4;
%tab.UPDRS=Groupcell5';
tab.ClinicalThresholdz=tab.ClinicalThreshold;
tab.SideEffectThresholdz=tab.SideEffectThreshold;
tab.TherapeuticWindowz=tab.TherapeuticWindow;
tab.SpatialCorrelationz=tab.SpatialCorrelation;

tab.ClinicalThresholdz(~isnan(tab.ClinicalThreshold),:)=zscore(tab.ClinicalThreshold(~isnan(tab.ClinicalThreshold),:));
tab.SideEffectThresholdz(~isnan(tab.SideEffectThreshold),:)=zscore(tab.SideEffectThreshold(~isnan(tab.SideEffectThreshold),:));
tab.TherapeuticWindowz(~isnan(tab.TherapeuticWindow),:)=zscore(tab.TherapeuticWindow(~isnan(tab.TherapeuticWindow),:));
tab.SpatialCorrelationz(~isnan(tab.SpatialCorrelation),:)=zscore(tab.SpatialCorrelation(~isnan(tab.SpatialCorrelation),:));

clinraw= fitlme(tab,'ClinicalThreshold~SpatialCorrelation+(1|Hemi)')
diary([savepath,filesep,'lmeClinicalThreshold_raw.txt'])
clinraw
diary off
sidraw=fitlme(tab,'SideEffectThreshold~SpatialCorrelation+(1|Hemi)')
diary([savepath,filesep,'lmeSideEffect_raw.txt'])
sidraw
diary off
windraw= fitlme(tab,'TherapeuticWindow~SpatialCorrelation+(1|Hemi)')
diary([savepath,filesep,'lmeWindow_raw.txt'])
windraw
diary off


%+ (-1 + X1 | g1)
clinrawranslope= fitlme(tab,'ClinicalThreshold~SpatialCorrelation+(1|Hemi)+(-1+SpatialCorrelation|Hemi) ');
diary([savepath,filesep,'lmeClinicalThreshold_raw_ranslope.txt'])
clinrawranslope
diary off
sidrawranslope=fitlme(tab,'SideEffectThreshold~SpatialCorrelation+(1|Hemi)+(-1+SpatialCorrelation|Hemi)');
diary([savepath,filesep,'lmeSideEffect_raw_ranslope.txt'])
sidrawranslope
diary off
windrawranslope= fitlme(tab,'TherapeuticWindow~SpatialCorrelation+(1|Hemi)+(-1+SpatialCorrelation|Hemi)');
diary([savepath,filesep,'lmeWindow_raw_ranslope.txt'])
windrawranslope
diary off


clin=fitlme(tab,'ClinicalThresholdz~SpatialCorrelationz+(1|Hemi)');
diary([savepath,filesep,'lmeClinicalThreshold.txt'])
clin
diary off
sid=fitlme(tab,'SideEffectThresholdz~SpatialCorrelationz+(1|Hemi)');
diary([savepath,filesep,'lmeSideEffect.txt'])
sid
diary off
wind=fitlme(tab,'TherapeuticWindowz~SpatialCorrelationz+(1|Hemi)');
diary([savepath,filesep,'lmeWindow.txt'])
wind
diary off

% Random slopes
clinranslope=fitlme(tab,'ClinicalThresholdz~SpatialCorrelationz+(1|Hemi)+(-1+SpatialCorrelation|Hemi)');
diary([savepath,filesep,'lmeClinicalThreshold_ranslope.txt'])
clinranslope
diary off
sidranslope=fitlme(tab,'SideEffectThresholdz~SpatialCorrelationz+(1|Hemi)+(-1+SpatialCorrelation|Hemi)');
diary([savepath,filesep,'lmeSideEffect_ranslope.txt'])
sidranslope
diary off
windranslope=fitlme(tab,'TherapeuticWindowz~SpatialCorrelationz+(1|Hemi)+(-1+SpatialCorrelation|Hemi)');
diary([savepath,filesep,'lmeWindow_ranslope.txt'])
windranslope
diary off
hemicols=mandrill(length(unique(Groupcell4)));


[h,R,p,g]=ea_corrplot(tab.TherapeuticWindow,tab.SpatialCorrelation,'no',{'Out of Sample Validation','Therapeutic Window (%)','Similarity to R-Matrix'},group4,[],hemicols,[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Matrix'}
h.Children(2).FontWeight="bold";
h.Children(2).TickDir= "out";
h.Children(2).FontSize=12;
h.Children(3).Children.String{3}=h.Children(3).Children.String{2};
h.Children(3).Children.Interpreter='tex';
if  wind.Coefficients.pValue(2)<0.001
    h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', sprintf('%.2e',wind.Coefficients.pValue(2))];
else
    h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', sprintf('%.2f',wind.Coefficients.pValue(2))];
end
axes(h.Children(2))

% el=lsline
% len=1:length(hemicols);
% el(1).Color=hemicols(5,:);
% el(2).Color=hemicols(4,:);
% el(3).Color=hemicols(3,:);
% el(4).Color=hemicols(2,:);
% el(5).Color=hemicols(1,:);
% for pp=len
% el(pp).LineStyle=':';
% el(pp).LineWidth=1;
% el(pp).XData=[min(tab.TherapeuticWindow),max(tab.TherapeuticWindow)];
% end
h.Position(4)=h.Position(3)+h.Position(3)*0.025;
saveas(h,[savepath,'LMEresults_TherapeuticWindow','.png']);
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')

export_fig([savepath,'LMEresults_TherapeuticWindow_transparent.png'], '-png', '-transparent');

[h,R,p,g]=ea_corrplot(tab.ClinicalThreshold,tab.SpatialCorrelation,'no',{'Out of Sample Validation','Clinical Effect Threshold (mA)','Similarity to R-Matrix'},group4,[],hemicols,[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Matrix'}
h.Children(2).FontWeight="bold";
h.Children(2).TickDir= "out";
h.Children(2).FontSize=12;
h.Children(3).Children.String{3}=h.Children(3).Children.String{2};
h.Children(3).Children.Interpreter='tex';
if  clin.Coefficients.pValue(2)<0.001
    h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', sprintf('%.2e',clin.Coefficients.pValue(2))];
else
    h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', sprintf('%.2f',clin.Coefficients.pValue(2))];
end
hold on

% axes(h.Children(2))
% el=lsline
% ylim([-0.5 0.3])
% len=1:length(hemicols);
% el(1).Color=hemicols(5,:);
% el(2).Color=hemicols(4,:);
% el(3).Color=hemicols(3,:);
% el(4).Color=hemicols(2,:);
% el(5).Color=hemicols(1,:);
% for pp=len
% el(pp).LineStyle=':';
% el(pp).LineWidth=1;
% el(pp).XData=[min(tab.ClinicalThreshold),max(tab.ClinicalThreshold)];
% end
h.Position(4)=h.Position(3)+h.Position(3)*0.025;
axes(h.Children(2))
set (gca,'xdir','reverse')
saveas(h,[savepath,'LMEresults_ClinicalThreshold','.png']);
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')

export_fig([savepath,'LMEresults_ClinicalThreshold_transparent.png'], '-png', '-transparent');

% %% Prediction plots:
% 
[h,R,p,g]=ea_corrplot(response(clinraw),fitted(clinraw),'no',...
    {'Out of Sample Validation','Clinical Effect Threshold (mA)','Fitted Effect Threshold (mA)'},...
    group4,[],hemicols,[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Matrix'}
h.Children(2).FontWeight="bold";
h.Children(2).TickDir= "out";
h.Children(2).FontSize=12;
% axes(h.Children(2));
% % text(max(response(clinraw))-max(response(clinraw))*0.6+min(response(clinraw)),...
% %     (max(fitted(clinraw)))-max(fitted(clinraw))*0.8+min(fitted(clinraw)),{['R^{2} = ',sprintf('%.2f',clinraw.Rsquared.Ordinary)],...
% % ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', ...
% %     sprintf('%.2e',clin.Coefficients.pValue(2))]},...
% %     'FontWeight','bold','FontSize',14);

% axes(h.Children(2))
% hold on
% 
% el=lsline
% len=1:length(hemicols);
% % el(1).Color=hemicols(5,:);
% % el(2).Color=hemicols(4,:);
% % el(3).Color=hemicols(3,:);
% % el(4).Color=hemicols(2,:);
% % el(5).Color=hemicols(1,:);
% for pp=len
% el(pp).LineStyle=':';
% el(pp).LineWidth=1;
% el(pp).XData=[min(response(clinraw)),max(response(clinraw))];
% el(pp).Color=hemicols(pp,:);
% end
% ylim([0 6]);
RMS=sqrt(mean(clinraw.residuals.^2));
MAE= mean(abs(clinraw.predict-clinraw.response));
% 
h.Children(2).FontWeight="bold";
h.Children(2).TickDir= "out";
h.Children(2).FontSize=12;
h.Children(3).Children.String{3}=[];%['R^{2} = ',sprintf('%.2f',windraw.Rsquared.Ordinary)];
h.Children(3).Children.Interpreter='tex';
% if  clin.Coefficients.pValue(2)<0.001
    h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', sprintf('%.2e',clin.Coefficients.pValue(2))];
h.Children(3).Children.String{3}= ['R^{2} = ',sprintf('%.2f',clinraw.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];

% else
%     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', sprintf('%.2f',clin.Coefficients.pValue(2)),'; R^{2} = ',sprintf('%.2f',clinraw.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];
% h.Children(3).Children.String{3}= ['R^{2} = ',sprintf('%.2f',clinraw.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];
% 
% end
h.Position(4)=h.Position(3)+h.Position(3)*0.025;
saveas(h,[savepath,'LMEresults_ClinicalThreshold_pred','.png']);
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')
export_fig([savepath,'LMEresults_ClinicalThreshold_pred_transparent.png'], '-png', '-transparent');


[h,R,p,g]=ea_corrplot(response(windraw),fitted(windraw),'no',{'Out of Sample Validation','Therapeutic Window (%)','Fitted Therapeutic Window (%)'},group4,[],hemicols,[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Matrix'}

h.Children(2).FontWeight="bold";
h.Children(2).TickDir= "out";
h.Children(2).FontSize=12;
axes(h.Children(2));
% text(max(response(windraw))*0.5,max(fitted(windraw))*0.25,{['R^{2} = ',sprintf('%.2f',windraw.Rsquared.Ordinary)],...
% ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', ...
%     sprintf('%.2e',wind.Coefficients.pValue(2))]},...
%     'FontWeight','bold','FontSize',14);

% axes(h.Children(2))
% 
% el=lsline
% len=1:length(hemicols);
% % el(1).Color=hemicols(5,:);
% % el(2).Color=hemicols(4,:);
% % el(3).Color=hemicols(3,:);
% % el(4).Color=hemicols(2,:);
% % el(5).Color=hemicols(1,:);
% for pp=len
% el(pp).LineStyle=':';
% el(pp).LineWidth=1;
% el(pp).XData=[min(response(windraw)),max(response(windraw))];
% el(pp).Color=hemicols(pp,:);
% end
RMS=sqrt(mean(windraw.residuals.^2));
MAE= mean(abs(windraw.predict-windraw.response));

h.Children(2).FontWeight="bold";
h.Children(2).TickDir= "out";
h.Children(2).FontSize=12;
%h.Children(3).Children.String{3}=[];%['R^{2} = ',sprintf('%.2f',windraw.Rsquared.Ordinary)];
h.Children(3).Children.Interpreter='tex';
% if  clin.Coefficients.pValue(2)<0.001
    h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', sprintf('%.2e',wind.Coefficients.pValue(2))];
h.Children(3).Children.String{3}= ['R^{2} = ',sprintf('%.2f',wind.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];

% else
%     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', sprintf('%.2f',wind.Coefficients.pValue(2))];
% h.Children(3).Children.String{3}= ['R^{2} = ',sprintf('%.2f',wind.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];
% 
% end
h.Position(4)=h.Position(3)+h.Position(3)*0.025;
saveas(h,[savepath,'LMEresults_TherapeuticWindow_pred','.png']);
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')
export_fig([savepath,'LMEresults_TherapeuticWindow_pred_transparent.png'], '-png', '-transparent');


