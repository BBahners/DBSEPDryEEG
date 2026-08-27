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
rmappath=[dryeeg_data_path(repoRoot, 'Rmaps/LOOCV_leftchannels'), filesep]; % results path
%savepath=dryeeg_data_path(repoRoot, 'Rmaps/Figures/Outsample_omnidirectional_trialvalidation_2/');
savepath=dryeeg_result_path(repoRoot, 'validation_recordingtime_hitratios');
if ~exist(savepath,'file')
    mkdir(savepath);
end
% Machine-specific addpath removed; see config/local_paths.m
if ~exist(savepath,'file')
    mkdir(savepath);
end
load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat')); 
load(dryeeg_data_path(repoRoot, 'EPmaps_Contacts_trials/channels.mat')); 
corrtype='Spearman'; % use which type of correlation metric - for fMRI could also do Pearson
load([rmappath,'R_00.mat']); % load Rmap
rmaptime=load([rmappath,'time.mat']); % load time definition for Rmap
fptime=load(dryeeg_data_path(repoRoot, 'EPmaps_Contacts_alltrials/timevec.mat'));

Subjects={...
'S027',...
'S031',...
'S032',...
 'S033',...
    'S034',...
    'S035',...
    'S036',...
    'S037',...
    'S038',...
    'S039',...
    'S040',...
};
load(dryeeg_data_path(repoRoot, 'TrialValidationAnalysis/AverageMaps.mat'))
Avg=Avglarge;
Avg=horzcat(Avg{:});
load(dryeeg_data_path(repoRoot, 'TrialValidationAnalysis/AverageMapsNames.mat'))
load(dryeeg_data_path(repoRoot, 'TrialValidationAnalysis/AverageMapsTrials.mat'))
AvgName=horzcat(AvglargeName{:});

triallength={1:60, 1:100, 1:200, 1:300, 1:400, 1:500, 1:600, 1:700, 1:800, 1:900};


tbs=readtable(dryeeg_data_path(repoRoot, 'metadata/DryEEGOutsample_nonans_new.xlsx'));
files={...
    '027_stn_le_8mA_contact8',...
    '027_stn_le_8mA_contact9'	,...
    '027_stn_le_8mA_contact10',...
    '027_stn_le_8mA_contact11',...
    '027_stn_ri_8mA_contact0',...
    '027_stn_ri_8mA_contact1',...
    '027_stn_ri_8mA_contact2',...
    '027_stn_ri_8mA_contact3',...
    '031_stn_le_6mA_contact8',...
    '031_stn_le_6mA_contact9',...
    '031_stn_le_6mA_contact10',...
    '031_stn_le_6mA_contact11',...
    '032_stn_le_8mA_contact8',...
    '032_stn_le_8mA_contact9',...
    '032_stn_le_8mA_contact10',...
    '032_stn_le_8mA_contact11',...
    '032_stn_ri_8mA_contact0',...
    '032_stn_ri_8mA_contact1',...
    '032_stn_ri_8mA_contact2',...
    '032_stn_ri_8mA_contact3',...
    '033_stn_le_contact0',...
    '033_stn_le_contact1',...
    '033_stn_le_contact2',...
    '033_stn_le_contact3',...
    '033_stn_ri_contact0',...
    '033_stn_ri_contact1',...
    '033_stn_ri_contact2',...
    '033_stn_ri_contact3',...
    '034_stn_le_contact0',...
    '034_stn_le_contact1',...
    '034_stn_le_contact2',...
    '034_stn_le_contact3',...
    '034_stn_ri_contact0',...
    '034_stn_ri_contact1',...
    '034_stn_ri_contact2',...
    '034_stn_ri_contact3',...
     '035_stn_le_contact0',...
    '035_stn_le_contact1',...
    '035_stn_le_contact2',...
    '035_stn_le_contact3',...
    '035_stn_ri_contact0',...
    '035_stn_ri_contact1',...
    '035_stn_ri_contact2',...
    '035_stn_ri_contact3',...
     '036_stn_le_contact0',...
    '036_stn_le_contact1',...
    '036_stn_le_contact2',...
    '036_stn_le_contact3',...
    '036_stn_ri_contact0',...
    '036_stn_ri_contact1',...
    '036_stn_ri_contact2',...
    '036_stn_ri_contact3',...
     '037_stn_le_contact0',...
    '037_stn_le_contact1',...
    '037_stn_le_contact2',...
    '037_stn_le_contact3',...
    '037_stn_ri_contact0',...
    '037_stn_ri_contact1',...
    '037_stn_ri_contact2',...
    '037_stn_ri_contact3',...
     '038_stn_le_contact0',...
    '038_stn_le_contact1',...
    '038_stn_le_contact2',...
    '038_stn_le_contact3',...
    '038_stn_ri_contact0',...
    '038_stn_ri_contact1',...
    '038_stn_ri_contact2',...
    '038_stn_ri_contact3',...
     '039_stn_le_contact0',...
    '039_stn_le_contact1',...
    '039_stn_le_contact2',...
    '039_stn_le_contact3',...
    '039_stn_ri_contact0',...
    '039_stn_ri_contact1',...
    '039_stn_ri_contact2',...
    '039_stn_ri_contact3',...
     '040_stn_le_contact0',...
    '040_stn_le_contact1',...
    '040_stn_le_contact2',...
    '040_stn_le_contact3',...
    '040_stn_ri_contact0',...
    '040_stn_ri_contact1',...
    '040_stn_ri_contact2',...
    '040_stn_ri_contact3',...
    };



effectthreshold=tbs.effectthreshold;
sidethreshold=tbs.sidethreshold;
effectlog=tbs.effectlog;




for ifil=1:length(files)
idx(ifil,:)=find(contains(string(char(AvgName{1,:})),files{ifil}));
end

Avg=Avg(:,idx);
channels(:,33:35)=[];
for itrials=1:length(triallength) 





%% Prepare Channel Order & Channel Flip

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
% %% Plot R-map
% figure;
% Chmap=repmat(1:32,[32,1])';
% imagesc(rmaptime.time,1:32,Rmap(chanord,:))
% colormap(cmap)
% clim([-0.4,0.4]);

%% 0. create a matrix of connectivity files
for pt=1:size(files,2) % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
    tmp=Avg{itrials,pt};
    if contains(files{pt},'le')
        L=tmp(1:32,:);
        L(ridx,:)=NaN;
        patConnectivityFiles{pt,1}=L;
        patConnectivityFileName{pt,1}=files{pt};
    elseif contains(files{pt},'ri')
        Ri=tmp(1:32,:);
        Ri(lidx,:)=NaN;
        patConnectivityFiles{pt,1}=Ri(flip,:);
        patConnectivityFileName{pt,1}=files{pt};
    end
end

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

for hh=1:height(tbs)
Groupcell4{hh}=char([tbs.SubID{hh},tbs.side{hh}]);
end
%Groupcell4=tbs.SubID;

group4.idx=Groupcell4';
group4.tag='Hemisphere';

% Groupcell5=deltaperc;
% group5.idx=Groupcell5';
% group5.tag='';

tab.Hemi=Groupcell4';
tab.ClinicalThresholdz=zscore(tab.ClinicalThreshold);
tab.SideEffectThresholdz=zscore(tab.SideEffectThreshold);
tab.TherapeuticWindowz=zscore(tab.TherapeuticWindow);
tab.SpatialCorrelationz=zscore(tab.SpatialCorrelation);
tab.effectlog=effectlog;

clinraw{itrials}= fitlme(tab,'ClinicalThreshold~SpatialCorrelation+(1|Hemi)')
diary([savepath,filesep,'lmeClinicalThreshold_raw.txt'])
clinraw{itrials}
diary off
windraw{itrials}= fitlme(tab,'TherapeuticWindow~SpatialCorrelation+(1|Hemi)')
diary([savepath,filesep,'lmeWindow_raw.txt'])
windraw{itrials}
diary off

clin{itrials}=fitlme(tab,'ClinicalThresholdz~SpatialCorrelationz+(1|Hemi)');
diary([savepath,filesep,'lmeClinicalThreshold.txt'])
clin{itrials}
diary off

wind{itrials}=fitlme(tab,'TherapeuticWindowz~SpatialCorrelationz+(1|Hemi)');
diary([savepath,filesep,'lmeWindow.txt'])
wind{itrials}
diary off

hemicols=mandrill(length(unique(Groupcell4)));


% [h,Rwin{itrials},pwin{itrials},g]=ea_corrplot(tab.TherapeuticWindow,tab.SpatialCorrelation,'no',{'Out of Sample Validation','Therapeutic Window (%)','Similarity to R-Matrix'},group4,[],hemicols,[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Matrix'}
% h.Children(2).FontWeight="bold";
% h.Children(2).TickDir= "out";
% h.Children(2).FontSize=14;
% % h.Children(3).Children.String{3}=h.Children(3).Children.String{2};
% % h.Children(3).Children.Interpreter='tex';
% % if  wind.Coefficients.pValue(2)<0.001
% %     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', sprintf('%.2e',wind.Coefficients.pValue(2))];
% % else
% %     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', sprintf('%.2f',wind.Coefficients.pValue(2))];
% % end
%  axes(h.Children(2))
% 
% % el=lsline
% % len=1:length(hemicols);
% % el(1).Color=hemicols(5,:);
% % el(2).Color=hemicols(4,:);
% % el(3).Color=hemicols(3,:);
% % el(4).Color=hemicols(2,:);
% % el(5).Color=hemicols(1,:);
% % for pp=len
% % el(pp).LineStyle=':';
% % el(pp).LineWidth=1;
% % el(pp).XData=[min(tab.TherapeuticWindow),max(tab.TherapeuticWindow)];
% % end
% h.Position(4)=h.Position(3)+h.Position(3)*0.025;
% saveas(h,[savepath,'LMEresults_TherapeuticWindow',num2str(max(triallength{itrials})),'.png']);
% set(gcf, 'Color', 'none');
% set(gca, 'Color', 'none');
% set(gca,'Box','off','Color','none')
% 
% export_fig([savepath,'LMEresults_TherapeuticWindow_transparent.png'], '-png', '-transparent');

% [h,Rclin{itrials},pclin{itrials},g]=ea_corrplot(tab.ClinicalThreshold,tab.SpatialCorrelation,'no',{'Out of Sample Validation','Clinical Effect Threshold (mA)','Similarity to R-Matrix'},group4,[],hemicols,[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Matrix'}
% h.Children(2).FontWeight="bold";
% h.Children(2).TickDir= "out";
% h.Children(2).FontSize=14;
% 
% 
% % h.Children(3).Children.String{3}=h.Children(3).Children.String{2};
% % h.Children(3).Children.Interpreter='tex';
% % if  clin.Coefficients.pValue(2)<0.001
% %     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', sprintf('%.2e',clin.Coefficients.pValue(2))];
% % else
% %     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', sprintf('%.2f',clin.Coefficients.pValue(2))];
% % end
% % axes(h.Children(2))
% % el=lsline
% % len=1:length(hemicols);
% % el(1).Color=hemicols(5,:);
% % el(2).Color=hemicols(4,:);
% % el(3).Color=hemicols(3,:);
% % el(4).Color=hemicols(2,:);
% % el(5).Color=hemicols(1,:);
% % for pp=len
% % el(pp).LineStyle=':';
% % el(pp).LineWidth=1;
% % el(pp).XData=[min(tab.ClinicalThreshold),max(tab.ClinicalThreshold)];
% % end
% h.Position(4)=h.Position(3)+h.Position(3)*0.025;
% axes(h.Children(2))
% set (gca,'xdir','reverse')
% saveas(h,[savepath,'LMEresults_ClinicalThreshold',num2str(max(triallength{itrials})),'.png']);
% set(gcf, 'Color', 'none');
% set(gca, 'Color', 'none');
% set(gca,'Box','off','Color','none')

%export_fig([savepath,'LMEresults_ClinicalThreshold_transparent.png'], '-png', '-transparent');



mdl{itrials} = fitglme(tab, 'effectlog ~ SpatialCorrelationz+ (1|Hemi)', ...
'Distribution', 'Binomial', ...
'Link', 'logit');

% Suppose we have the following data for 2 subjects (4 predictions each)
Y_true = effectlog% True labels for all predictions
%Y_pred = tab.SpatialCorrelation;  % Predicted probabilities for all predictions
Y_pred = fitted(mdl{itrials});  % Predicted probabilities for all predictions





%Y_pred = fitted(clinraw{itrials});
% Generate the ROC curve
[X,Y,T,AUC{itrials}] = perfcurve(Y_true,Y_pred,1);




% % Plot the ROC curve
% figure;
% plot(X,Y);
% xlabel('False positive rate');
% ylabel('True positive rate');
% title(['ROC Curve (AUC = ' num2str(AUC{itrials}) ')']);
% grid on;
% 
% youden_index = Y + (1 - X) - 1
% Youden(itrials)=max(youden_index);


% Cummulative Hit Ratio
tbs.Comb=strcat(string(tbs.SubID),string(tbs.side));
electrode_ids=categorical(tbs.Comb);
active_contact_flags=tbs.effectlog;

predicted_windows=zscore(tab.SpatialCorrelation);
actual_windows=zscore(tab.TherapeuticWindow);

% % 1. Pearson Correlation
% pearson_r = corr(predicted_windows, actual_windows, 'Type', 'Pearson');
% fprintf('Pearson correlation: %.4f\n', pearson_r);

% 2. Ranking contacts by electrode and computing cumulative hit ratio
unique_electrodes = unique(electrode_ids);
num_electrodes = numel(unique_electrodes);

max_rank = 4; % assume max 8 contacts per electrode
cumulative_hits = zeros(max_rank, 1);

for i = 1:num_electrodes
    eid = unique_electrodes(i);
    idx = electrode_ids == eid;
    
    predicted = predicted_windows(idx);
    is_active = active_contact_flags(idx);
    
    [~, sorted_idx] = sort(predicted, 'descend');
    active_pos = find(is_active(sorted_idx));
    
    if isempty(active_pos)
        continue;
    end

    for r = active_pos:max_rank
        cumulative_hits(r) = cumulative_hits(r) + 1;
    end
end

cumulative_hit_ratio{itrials} = cumulative_hits / num_electrodes;

end


% 3. Empirical null distribution (random ranking)
rng('default');  % For reproducibility
num_shuffles = 10000;
null_dist = zeros(max_rank, num_shuffles);

for s = 1:num_shuffles
    hits = zeros(max_rank, 1);
    
    for i = 1:num_electrodes
        eid = unique_electrodes(i);
        idx = electrode_ids == eid;
        
        is_active = active_contact_flags(idx);
        num_contacts = sum(idx);
        
        rand_order = randperm(num_contacts);
        active_pos = find(is_active(rand_order));
        
        if isempty(active_pos)
            continue;
        end
        
        for r = active_pos:max_rank
            hits(r) = hits(r) + 1;
        end
    end
    
    null_dist(:, s) = hits / num_electrodes;
end

% Compute 95th percentile of the null distribution
null_95 = prctile(null_dist, 95, 2);
cmap2=mandrill(length(triallength));

% Compute null mean and Â±1 standard deviation
% null_mean = mean(null_dist, 2);
% null_std = std(null_dist, 0, 2);  % std across shuffles
% null_lower = null_mean - null_std;
% null_upper = null_mean + null_std;
null_mean = mean(null_dist, 2);
null_lower = prctile(null_dist, 2.5, 2);
null_upper = prctile(null_dist, 97.5, 2);

 x = 1:max_rank;
% Plot
figure;
hold on;

% Confidence interval patch (no border)
cp=ciplot(null_lower*100, null_upper*100, x', [0.8 0.8 0.8]);  % Light gray fill
cp.EdgeColor='none';

% Null mean line and dots
plot(x, null_mean*100, 'k--', 'LineWidth', 1.5);
plot(x(1:4), null_mean(1:4)*100, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);  % Dots

for itrials=1:length(triallength) 
% Check significance
significant{itrials} = cumulative_hit_ratio{itrials} > null_95;
% Observed line and dots
plot(x, cumulative_hit_ratio{itrials}*100, 'Color',cmap2(itrials,:), 'LineWidth', 2);
plot(x(1:4), cumulative_hit_ratio{itrials}(1:4)*100, 'o', 'MarkerFaceColor', cmap2(itrials,:), 'MarkerEdgeColor', cmap2(itrials,:), 'MarkerSize', 6);  % Dots
% v=x(1:4);
% plot(v(significant{itrials}),cumulative_hit_ratio{itrials}(significant{itrials}),'*')
end

xlabel('Contact Rank');
ylabel('Cumulative Hit Ratio (%)');
legend('Null 95% CI', 'Null Mean', '', '10s','','16s','','33s','','50s','','66s','','83s','','100s','','116s','','133s','','150s','','Location','Southeast');
%title('Cumulative Hit Ratio vs Null Distribution');
xlim([0.5 4.5])
xticks([1 2 3 4]);
%grid on;

ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.Legend.Box='off';
ax.FontWeight="bold";
 
 ax.FontSize=14;
f=gcf;
f.Color='white';
f.Position(4)=f.Position(3)+f.Position(3)*0.025;
 
saveas(f,[savepath,'CumulativeHitRatio_Trials.png']);

set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')

export_fig([savepath,'CumulativeHitRatio_Trials_transp.png'], '-png', '-transparent','-r600');
% for kk=1:length(Rclin)
% Rclinical(kk)=Rclin{kk}.spearman;
% pclinical(kk)=pclin{kk}.spearman;
% end
% 
% for kk=1:length(Rclin)
% Rwindow(kk)=Rwin{kk}.spearman;
% pwindow(kk)=pwin{kk}.spearman;
% end
for kk=1:length(clin)
lmeclin(kk)=clin{kk}.Coefficients(2, 2).Estimate;
lmeclinp(kk)=clin{kk}.Coefficients(2, 6).pValue;
lmeclinrsq(kk)=clin{kk}.Rsquared.Ordinary;
lmeclinbic(kk)=clin{kk}.ModelCriterion.BIC;
lmeclinaic(kk)=clin{kk}.ModelCriterion.AIC;

end

for kk=1:length(clin)
lmewin(kk)=wind{kk}.Coefficients(2, 2).Estimate;
lmewinp(kk)=wind{kk}.Coefficients(2, 6).pValue;
lmewinrsq(kk)=wind{kk}.Rsquared.Ordinary;
lmeclinbic(kk)=wind{kk}.ModelCriterion.BIC;
lmeclinaic(kk)=wind{kk}.ModelCriterion.AIC;
end

for kk=1:length(clin)
lmeeff(kk)=mdl{kk}.Coefficients(2, 2).Estimate;
lmeeffp(kk)=mdl{kk}.Coefficients(2, 6).pValue;
lmeeffrsq(kk)=mdl{kk}.Rsquared.Ordinary;
end

for iii=1:length(triallength)
trls(iii)=max(triallength{iii});
end
AUCcon=vertcat(AUC{:});
% figure;
% plot(trls/6,Rclinical*-1)
% hold on
% plot(pclinical)

% figure;
% plot(Rwindow)
% hold on
% plot(pwindow)

%trls=(trls*16)/60
h=figure;
plot(trls/6,lmeclin*-1,'Color',cmap(10,:),'LineWidth',2)
scatter(trls/6,lmeclin*-1,'Color',cmap(10,:),'LineWidth',2)
hold on
plot(trls/6,lmeclinp,'Color',cmap(220,:),'LineWidth',2)
plot(trls/6,lmewin,'Color',cmap(50,:),'LineWidth',2)
plot(trls/6,lmewinp,'Color',cmap(210,:),'LineWidth',2)
yline(0.05,'LineStyle','--')
ylabel('Î²_{std}');
xlabel('Recording Duration per DBS Contact (s)');
xlim([trls(2)/6,trls(end-1)/6])

ax=gca;
ax.TickDir='out';
ax.Box='off';
h.Color='white';


h=figure;

yline(0.05,'LineStyle','--','Color',[.5 .5 .5])

hold on

sigindx=min(find((lmewinp*10<0.05&lmeclinp*10<0.05)));
area(trls(sigindx:end)/6,repmat(10,[1,length(trls(sigindx:end))]),0,"FaceColor",[.8 .8 .8],'FaceAlpha',0.4)
text(67,0.07,'p = 0.05')


scatter(trls/6,lmewin,'MarkerFaceColor',cmap(210,:),'MarkerEdgeColor',cmap(210,:),'LineWidth',2)
plot(trls/6,lmewin,'Color',cmap(210,:),'LineWidth',2)

scatter(trls/6,lmeclin*-1,'MarkerFaceColor',cmap(10,:),'MarkerEdgeColor',cmap(10,:),'LineWidth',2)
plot(trls/6,lmeclin*-1,'Color',cmap(10,:),'LineWidth',2)

scatter(trls/6,lmewinp*10,'MarkerFaceColor',cmap(150,:),'MarkerEdgeColor',cmap(150,:),'LineWidth',2)
plot(trls/6,lmewinp*10,'Color',cmap(150,:),'LineWidth',2)

scatter(trls/6,lmeclinp*10,'MarkerFaceColor',cmap(70,:),'MarkerEdgeColor',cmap(70,:),'LineWidth',2)
plot(trls/6,lmeclinp*10,'Color',cmap(70,:),'LineWidth',2)

% scatter(trls/6,AUCcon,'MarkerFaceColor',cmap(10,:),'MarkerEdgeColor',cmap(10,:),'LineWidth',2)
% plot(trls/6,AUCcon,'Color',cmap(10,:),'LineWidth',2)

% scatter(trls/6,lmeclinbic,'MarkerFaceColor',cmap(210,:),'MarkerEdgeColor',cmap(10,:),'LineWidth',2)
% plot(trls/6,lmeclinbic,'Color',cmap(210,:),'LineWidth',2)
% 
% scatter(trls/6,lmeclinaic,'MarkerFaceColor',cmap(10,:),'MarkerEdgeColor',cmap(10,:),'LineWidth',2)
% plot(trls/6,lmeclinaic,'Color',cmap(10,:),'LineWidth',2)
% 
% scatter(trls/6,lmeeff,'MarkerFaceColor',cmap(10,:),'MarkerEdgeColor',cmap(10,:),'LineWidth',2)
% plot(trls/6,lmeeff,'Color',cmap(10,:),'LineWidth',2)
%scatter(trls/6,Youden)
%scatter(trls/6,lmewinp,'MarkerEdgeColor',cmap(150,:),'LineWidth',2)

ylabel('Î²_{std} / p_{bonf}');

%ylabel('Î²_{std}');
xlabel('Recording Duration / DBS Contact (s)');
xlim([trls(2)/6,trls(end-1)/6])
ylim([0 0.6])
legend('',...
    '',...
    'Î²_{std} Therapeutic Window',...
    '',...
    'Î²_{std} Clinical Threshold',...
    '',...
    'p_{bonf} Therapeutic Window',...
    '',...
    'p_{bonf} Clinical Threshold',...
     '')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.Legend.Box='off';
ax.FontSize=14;
ax.FontWeight='bold';

%ax.Legend.Position(1)=ax.Legend.Position(1)+0.1
% ax.Legend.Location='southeast'
% ax.Legend.Position(2)=ax.Legend.Position(2)+0.1
h.Color='white';

h.Position(4)=h.Position(3)+h.Position(3)*0.025;
saveas(h,[savepath,'TrialComparison','.png']);
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')

export_fig([savepath,'TrialComparison_transp.png'], '-png', '-transparent','-r600');
% %% Prediction plots:
% 
% [h,R,p,g]=ea_corrplot(response(clinraw),fitted(clinraw),'no',...
%     {'Out of Sample Validation','Empirical Effect Threshold (mA)','Fitted Effect Threshold (mA)'},...
%     group4,[],hemicols,[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Matrix'}
% h.Children(2).FontWeight="bold";
% h.Children(2).TickDir= "out";
% h.Children(2).FontSize=14;
% % axes(h.Children(2));
% % % text(max(response(clinraw))-max(response(clinraw))*0.6+min(response(clinraw)),...
% % %     (max(fitted(clinraw)))-max(fitted(clinraw))*0.8+min(fitted(clinraw)),{['R^{2} = ',sprintf('%.2f',clinraw.Rsquared.Ordinary)],...
% % % ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', ...
% % %     sprintf('%.2e',clin.Coefficients.pValue(2))]},...
% % %     'FontWeight','bold','FontSize',14);
% 
% axes(h.Children(2))
% hold on
% 
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
% el(pp).XData=[min(response(clinraw)),max(response(clinraw))];
% end
% RMS=sqrt(mean(clinraw.residuals.^2));
% MAE= mean(abs(clinraw.predict-clinraw.response));
% 
% h.Children(2).FontWeight="bold";
% h.Children(2).TickDir= "out";
% h.Children(2).FontSize=14;
% h.Children(3).Children.String{3}=[];%['R^{2} = ',sprintf('%.2f',windraw.Rsquared.Ordinary)];
% h.Children(3).Children.Interpreter='tex';
% if  clin.Coefficients.pValue(2)<0.001
%     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', sprintf('%.2e',clin.Coefficients.pValue(2))];
% h.Children(3).Children.String{3}= ['R^{2} = ',sprintf('%.2f',clinraw.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];
% 
% else
%     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',clin.Coefficients.Estimate(2)),', p = ', sprintf('%.2f',clin.Coefficients.pValue(2)),'; R^{2} = ',sprintf('%.2f',clinraw.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];
% h.Children(3).Children.String{3}= ['R^{2} = ',sprintf('%.2f',clinraw.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];
% 
% end
% h.Position(4)=h.Position(3)+h.Position(3)*0.025;
% saveas(h,[savepath,'LMEresults_ClinicalThreshold_pred','.png']);
% set(gcf, 'Color', 'none');
% set(gca, 'Color', 'none');
% set(gca,'Box','off','Color','none')
% export_fig([savepath,'LMEresults_ClinicalThreshold_pred_transparent.png'], '-png', '-transparent');
% 
% 
% [h,R,p,g]=ea_corrplot(response(windraw),fitted(windraw),'no',{'Out of Sample Validation','Empirical Therapeutic Window (%)','Fitted Therapeutic Window (%)'},group4,[],hemicols,[],'linear');%{'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Matrix'}
% 
% h.Children(2).FontWeight="bold";
% h.Children(2).TickDir= "out";
% h.Children(2).FontSize=14;
% % axes(h.Children(2));
% % text(max(response(windraw))*0.5,max(fitted(windraw))*0.25,{['R^{2} = ',sprintf('%.2f',windraw.Rsquared.Ordinary)],...
% % ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', ...
% %     sprintf('%.2e',wind.Coefficients.pValue(2))]},...
% %     'FontWeight','bold','FontSize',14);
% 
% axes(h.Children(2))
% 
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
% el(pp).XData=[min(response(windraw)),max(response(windraw))];
% end
% RMS=sqrt(mean(windraw.residuals.^2));
% MAE= mean(abs(windraw.predict-windraw.response));
% 
% h.Children(2).FontWeight="bold";
% h.Children(2).TickDir= "out";
% h.Children(2).FontSize=14;
% h.Children(3).Children.String{3}=[];%['R^{2} = ',sprintf('%.2f',windraw.Rsquared.Ordinary)];
% h.Children(3).Children.Interpreter='tex';
% if  clin.Coefficients.pValue(2)<0.001
%     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', sprintf('%.2e',wind.Coefficients.pValue(2))];
% h.Children(3).Children.String{3}= ['R^{2} = ',sprintf('%.2f',wind.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];
% 
% else
%     h.Children(3).Children.String{2}= ['Î²_{std} = ', sprintf('%.2f',wind.Coefficients.Estimate(2)),', p = ', sprintf('%.2f',wind.Coefficients.pValue(2))];
% h.Children(3).Children.String{3}= ['R^{2} = ',sprintf('%.2f',wind.Rsquared.Ordinary),'; RMS = ',sprintf('%.2f',RMS),'; MAE = ',sprintf('%.2f',MAE)];
% 
% end
% h.Position(4)=h.Position(3)+h.Position(3)*0.025;
% saveas(h,[savepath,'LMEresults_TherapeuticWindow_pred','.png']);
% set(gcf, 'Color', 'none');
% set(gca, 'Color', 'none');
% set(gca,'Box','off','Color','none')
% export_fig([savepath,'LMEresults_TherapeuticWindow_pred_transparent.png'], '-png', '-transparent');
% 

