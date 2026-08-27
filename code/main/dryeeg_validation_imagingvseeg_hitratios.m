%% Distance / Spatial Correlation vs Cumulative Hit Ratio
clear variables; close all;
% Cumulative-hit analysis updated for EEG, LFP, and imaging predictors.
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
% Machine-specific addpath removed; see config/local_paths.m

% Exact theoretical mean for random selection among four contacts. Set to
% "matchedActiveContacts" to preserve multiple acceptable contacts in the
% Monte Carlo null instead.
null_model = "singleContactTheoretical";

%% Sweetspot definitions
sweetspots = {
    {[12.5 -12.72 -5.38],[-12.5 -12.72 -5.38]}, ...   % Dembek et al. 2019
    {[12.58 -13.41 -5.87],[-12.58 -13.41 -5.87]}, ... % Caire et al. symmetric left
    {[10.83 -13.31 -7],[-11 -14 -7]} ...            % Akram et al.
    
};

sweetspot_labels = {'Dembek et al.','Caire et al.','Akram et al.'};

leadfold=[dryeeg_data_path(repoRoot, 'DryEEGLeadsrefined/derivatives/leaddbs'), filesep];
savepath=dryeeg_result_path(repoRoot, 'validation_imagingvseeg_hitratios');
if ~exist(savepath,'file')
    mkdir(savepath);
end

list = dir(leadfold); list(1:2)=[];

tab = readtable(dryeeg_data_path(repoRoot, ...
    'Rmaps/Figures/Outsample_omnidirectional_V2_noH/R_results_SpatialCorrs_bcn_sid_new_LFP_completed_new_values_with_mono_effectlog.xlsx'));
tab(85:end,:) = [];
tab2 = readtable(dryeeg_data_path(repoRoot, 'metadata/DryEEGOutsample_nonans_new_2.xlsx'));

tab.effectlog = tab2.effectlog2;
tab.effectclin = tab2.effectclin;

%% Load coordinates
for ii = 1:length(list)
    fils = dir([leadfold,list(ii).name,filesep,'reconstruction']);
    fidx = find(contains(string(char(fils.name)),'reconstruction'));
    tmp = load([fils(fidx).folder,filesep,fils(fidx).name]);

    for side = 1:2
        Coords{ii,side}(1,:) = tmp.reco.mni.coords_mm{side}(1,:);
        Coords{ii,side}(2,:) = mean(tmp.reco.mni.coords_mm{side}(2:4,:));
        Coords{ii,side}(3,:) = mean(tmp.reco.mni.coords_mm{side}(5:7,:));
        Coords{ii,side}(4,:) = tmp.reco.mni.coords_mm{side}(8,:);
    end
end

%% Build data tables for each sweet spot
results = struct();
for ss = 1:length(sweetspots)
    mnisweet = sweetspots{ss};

    % Compute distance to sweet spot
    for side = 1:2
        for pt = 1:length(list)
            for cont = 1:4
                Dist{pt,side}(cont,:) = norm(mnisweet{side}-Coords{pt,side}(cont,:))*-1;
            end
        end
    end

    % Collect behavioral data
    for ii = 1:length(list)
        % right
        sbidx = find(contains(string(char(tab.SubID)),list(ii).name(6:end)) & ...
                     contains(string(char(tab.side)),'right'));
        if isempty(sbidx)
            Clin{ii,1} = nan(4,1);
             Clin2{ii,1} = nan(4,1);
            Sim{ii,1}  = nan(4,1);
            Sub{ii,1}  = repmat({NaN},4,1);
            Sid{ii,1}  = repmat({NaN},4,1);
            Eff{ii,1}  = nan(4,1);
            LBeta{ii,1} = nan(4,1);
            HBeta{ii,1} = nan(4,1);
            Beta{ii,1} = nan(4,1);
        else
            Clin{ii,1} = tab(sbidx,:).TherapeuticWindow;
            Clin2{ii,1} = tab(sbidx,:).ClinicalThreshold;
            Sim{ii,1}  = tab(sbidx,:).SpatialCorrelation;
            Sub{ii,1}  = tab(sbidx,:).SubID;
            Sid{ii,1}  = tab(sbidx,:).side;
            Eff{ii,1}  = tab(sbidx,:).effectlog;
            LBeta{ii,1} = tab(sbidx,:).LowBeta_euclidean;
            HBeta{ii,1} = tab(sbidx,:).HighBeta_euclidean;
            Beta{ii,1} = tab(sbidx,:).Beta_euclidean;
        end

        % left
        sbidx = find(contains(string(char(tab.SubID)),list(ii).name(6:end)) & ...
                     contains(string(char(tab.side)),'left'));
        if isempty(sbidx)
            Clin{ii,2} = nan(4,1);
             Clin2{ii,2} = nan(4,1);
            Sim{ii,2}  = nan(4,1);
            Sub{ii,2}  = repmat({NaN},4,1);
            Sid{ii,2}  = repmat({NaN},4,1);
            Eff{ii,2}  = nan(4,1);
            LBeta{ii,2} = nan(4,1);
            HBeta{ii,2} = nan(4,1);
            Beta{ii,2} = nan(4,1);
        else
            Clin{ii,2} = tab(sbidx,:).TherapeuticWindow;
           Clin2{ii,2} = tab(sbidx,:).ClinicalThreshold;
            Sim{ii,2}  = tab(sbidx,:).SpatialCorrelation;
            Sub{ii,2}  = tab(sbidx,:).SubID;
            Sid{ii,2}  = tab(sbidx,:).side;
            Eff{ii,2}  = tab(sbidx,:).effectlog;
            LBeta{ii,2} = tab(sbidx,:).LowBeta_euclidean;
            HBeta{ii,2} = tab(sbidx,:).HighBeta_euclidean;
            Beta{ii,2} = tab(sbidx,:).Beta_euclidean;
        end
    end

    % Combine
    Distc = [vertcat(Dist{:,2}); vertcat(Dist{:,1})];
    Clinc = [vertcat(Clin{:,2}); vertcat(Clin{:,1})];
    Clinc2 = [vertcat(Clin2{:,2}); vertcat(Clin2{:,1})];
    Simc  = [vertcat(Sim{:,2});  vertcat(Sim{:,1})];
    Subc  = [vertcat(Sub{:,2});  vertcat(Sub{:,1})];
    Sidc  = [vertcat(Sid{:,2});  vertcat(Sid{:,1})];
    Effc  = [vertcat(Eff{:,2});  vertcat(Eff{:,1})];
    if ss == 1
        Betac = [vertcat(LBeta{:,2}); vertcat(LBeta{:,1})];
    elseif ss == 2
        Betac = [vertcat(HBeta{:,2}); vertcat(HBeta{:,1})];
    else
        Betac = [vertcat(Beta{:,2}); vertcat(Beta{:,1})];
    end

    tbs = table(Distc,Clinc,Clinc2,Simc,Subc,Sidc,Effc,Betac, ...
        'VariableNames', {'Distance','TherapeuticWindow', ...
        'ClinicalThreshold','SpatialCorrelation','SubID','side', ...
        'effectlog','Beta'});
    tbs(isnan(tbs.effectlog),:) = [];

    % Store table
    results(ss).tbs = tbs;
end

%% Compare five predictors
colors = mandrill(8);
colors = colors(end:-1:1,:);
colors(3:5,:) = colors(6:8,:);

plot_tables = {results(1).tbs,results(1).tbs,results(3).tbs, ...
    results(2).tbs,results(1).tbs};
plot_variables = {'SpatialCorrelation','Beta','Distance','Distance','Distance'};
plot_labels = {'EEG','LFP','Akram','Caire','Dembek'};

num_curves = numel(plot_labels);
cumhit_curves = nan(4,num_curves);
num_electrodes_curves = nan(num_curves,1);

for p = 1:num_curves
    tbs_pred = plot_tables{p};
    ids_pred = categorical(strcat(string(tbs_pred.SubID),string(tbs_pred.side)));
    predictor_vals = tbs_pred{:,plot_variables{p}};
    [cumhit_curves(:,p),null_mean,null_lower,null_upper, ...
        num_electrodes_curves(p)] = compute_hit_ratio( ...
        predictor_vals,ids_pred,tbs_pred.effectlog,4,null_model);
    fprintf('%s: n = %d electrodes\n',plot_labels{p}, ...
        num_electrodes_curves(p));
end

% Include rank zero so the source data and plotted curves share their origin.
x = (0:4)';
null_plot = [0;null_mean];
cumhit_plot = [zeros(1,num_curves);cumhit_curves];

figure; hold on;
h_null = plot(x,null_plot*100,'k--','LineWidth',1.5);
h_observed = gobjects(num_curves,1);
for p = 1:num_curves
    h_observed(p) = plot(x,cumhit_plot(:,p)*100,'-', ...
        'Color',colors(p,:),'LineWidth',2);
end

% Save the exact values used by the figure.
curve_table = table();
for p = 1:num_curves
    this_curve = table(repmat(string(plot_labels{p}),5,1),x, ...
        cumhit_plot(:,p)*100,repmat(num_electrodes_curves(p),5,1), ...
        'VariableNames',{'Series','ContactRank', ...
        'CumulativeHitRatioPercent','EligibleElectrodes'});
    if p == 1
        curve_table = this_curve;
    else
        curve_table = [curve_table;this_curve]; %#ok<AGROW>
    end
end
null_table = table(repmat("Null Mean",5,1),x,null_plot*100, ...
    repmat(num_electrodes_curves(3),5,1), ...
    'VariableNames',curve_table.Properties.VariableNames);
curve_table = [null_table;curve_table];
writetable(curve_table,[savepath,'CumulativeHitRatio_AllSweetspots.csv']);

% Sweet-spot and LFP mixed-effects summaries from the revised analysis.
for ss = 1:length(sweetspots)
    tbs_pred = results(ss).tbs;
    [coeffs{ss},pvals{ss}] = corr(tbs_pred.Beta, ...
        zscore(tbs_pred.TherapeuticWindow),'Type','spearman','Rows','complete');
    tbs_pred.TherapeuticWindowz = zscore(tbs_pred.TherapeuticWindow);
    tbs_pred.Distancez = zscore(tbs_pred.Distance);
    tbs_pred.ClinicalThresholdz = zscore(tbs_pred.ClinicalThreshold);
    tbs_pred.SpatialCorrelationz = zscore(tbs_pred.SpatialCorrelation);
    tbs_pred.Betaz = tbs_pred.Beta;
    beta_idx = isfinite(tbs_pred.Betaz);
    tbs_pred.Betaz(beta_idx) = zscore(tbs_pred.Beta(beta_idx));

    lmeW{ss} = fitlme(tbs_pred,'TherapeuticWindowz~Distancez+(1|SubID)');
    diary([savepath,'LME_sweetspot',sweetspot_labels{ss},'.txt']);
    disp(lmeW{ss}); diary off;

    if ss == 1
        lmeB = fitlme(tbs_pred,'TherapeuticWindowz~Betaz+(1|SubID)');
        diary([savepath,'LME_LFP_LowBeta.txt']); disp(lmeB); diary off;
        lmeBTT = fitlme(tbs_pred,'ClinicalThresholdz~Betaz+(1|SubID)');
        diary([savepath,'LME_LFP_LowBeta_therthresh.txt']); disp(lmeBTT); diary off;
        lmecomb = fitlme(tbs_pred, ...
            'TherapeuticWindowz~Betaz*SpatialCorrelationz+(1|SubID)');
        diary([savepath,'LME_LFP_LowBeta_comb.txt']); disp(lmecomb); diary off;
    end
end

xlabel('Contact Rank'); ylabel('Cumulative Hit Ratio (%)');
legend([h_null;h_observed],[{'Null Mean'};plot_labels(:)], ...
    'Location','Southeast');
set(gca,'TickDir','out','Box','off','FontWeight','bold','FontSize',14);
set(gcf,'Color','w'); xlim([0 4.5]); xticks(0:4);
ax = gca; ax.Legend.Box = 'off';
f = gcf; f.Position(4) = f.Position(3)+f.Position(3)*0.025;

saveas(gcf,[savepath,'CumulativeHitRatio_AllSweetspots.png']);
set(gcf,'Color','none'); set(gca,'Box','off','Color','none');
export_fig([savepath,'CumulativeHitRatio_AllSweetspots_transp.png'], ...
    '-png','-transparent','-r600');

%% Function for cumulative hit ratio
function [cumhit,null_mean,null_lower,null_upper,num_electrodes] = ...
    compute_hit_ratio(predicted,ids,active,max_rank,null_model)
    all_electrodes = unique(ids);
    valid_electrode = false(numel(all_electrodes),1);
    for i = 1:numel(all_electrodes)
        idx = ids == all_electrodes(i);
        valid_electrode(i) = sum(idx) == max_rank && ...
            all(isfinite(predicted(idx))) && all(isfinite(active(idx))) && ...
            any(active(idx) == 1);
    end
    unique_electrodes = all_electrodes(valid_electrode);
    num_electrodes = numel(unique_electrodes);
    if num_electrodes == 0
        error('No complete electrodes with at least one active contact.');
    end
    cumulative_hits = zeros(max_rank, 1);

    % Observed
    for i = 1:num_electrodes
        idx = ids == unique_electrodes(i);
        pred = predicted(idx); is_act = active(idx) == 1;
        [~, sorted_idx] = sort(pred,'descend','MissingPlacement','last');
        active_pos = find(is_act(sorted_idx),1,'first');
        for r = active_pos:max_rank
            cumulative_hits(r) = cumulative_hits(r) + 1;
        end
    end
    cumhit = cumulative_hits / num_electrodes;

    % Null distribution uses the same eligible-electrode denominator.
    num_shuffles = 10000; null_dist = zeros(max_rank,num_shuffles);
    switch lower(string(null_model))
        case "matchedactivecontacts"
            for s = 1:num_shuffles
                hits = zeros(max_rank,1);
                for i = 1:num_electrodes
                    idx = ids == unique_electrodes(i);
                    is_act = active(idx) == 1;
                    active_pos = find(is_act(randperm(max_rank)),1,'first');
                    hits(active_pos:max_rank) = hits(active_pos:max_rank)+1;
                end
                null_dist(:,s) = hits/num_electrodes;
            end
            null_mean = mean(null_dist,2);
        case "singlecontacttheoretical"
            for s = 1:num_shuffles
                hits = zeros(max_rank,1);
                for i = 1:num_electrodes
                    active_pos = randi(max_rank);
                    hits(active_pos:max_rank) = hits(active_pos:max_rank)+1;
                end
                null_dist(:,s) = hits/num_electrodes;
            end
            null_mean = (1:max_rank)'/max_rank;
        otherwise
            error('Unknown null_model: %s',null_model);
    end
    null_lower = prctile(null_dist,2.5,2);
    null_upper = prctile(null_dist,97.5,2);
end

