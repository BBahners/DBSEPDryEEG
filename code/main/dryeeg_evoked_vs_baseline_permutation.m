%% Test DBS-evoked EEG responses against the pre-stimulation baseline
% Paired tests are followed by subject-level sign-flip permutations, joint
% Benjamini-Hochberg correction, and a minimum temporal-run requirement.
clear variables; close all;
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot,'functions'),'-begin');
dryeeg_setup(repoRoot);

%% Configuration
data_path = dryeeg_data_path(repoRoot,'EPmaps_Discovery');
results_path = dryeeg_result_path(repoRoot,'baseline_statistics');
metadata_file = dryeeg_data_path(repoRoot, ...
    'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp.xlsx');
time_file = dryeeg_data_path(repoRoot,'timevec.mat');
channel_file = dryeeg_data_path(repoRoot,'channels.mat');
exclude_subject_ids = "P027";
epoch_window = [-0.050,0.200];
baseline_window = [-0.050,-0.006];
test_window = [0.010,0.200];
initial_alpha = 0.05;
fdr_alpha = 0.005;
n_permutations = 10000;
minimum_duration = 0.005;
random_seed = 20260722;
permutation_batch_size = 250;

%% Load study inputs
loaded_time = load(time_file);
assert(isfield(loaded_time,'time'),'timevec.mat must contain time.');
time_full = double(loaded_time.time(:)');
metadata = readtable(metadata_file);
subject_ids = string(metadata.SubID);
subject_ids = subject_ids(~ismember(subject_ids,exclude_subject_ids));

epoch_idx = find(time_full>=epoch_window(1) & time_full<=epoch_window(2));
time = time_full(epoch_idx);
baseline_idx = find(time>=baseline_window(1) & time<=baseline_window(2));
test_idx = find(time>=test_window(1) & time<=test_window(2));
assert(~isempty(baseline_idx) && ~isempty(test_idx), ...
    'Baseline and test windows must both contain samples.');
assert(max(baseline_idx)<min(test_idx), ...
    'Baseline and post-stimulation windows must not overlap.');
sampling_interval = median(diff(time));
minimum_samples = max(1,ceil(minimum_duration/sampling_interval));

loaded_channels = load(channel_file);
channels = loaded_channels.channels;
channels(:,33:35) = [];
channel_labels = strings(1,numel(channels));
for channel_idx = 1:numel(channels)
    channel_labels(channel_idx) = string(channels(channel_idx).Name);
end
n_channels = numel(channel_labels);

channel_side = {'right','right','right','right','right','right','right', ...
    'right','right','right','left','left','z','right','z','right','left', ...
    'left','z','right','right','z','left','left','left','left','left', ...
    'left','left','left','left','left'};
channel_flip = {'P7','T7','CP5','FC5','F7','F3','C3','P3','AF3','Fp1', ...
    'Fp2','AF4','Fz','FC1','Cz','CP1','PO4','O2','Oz','O1','PO3','Pz', ...
    'CP2','FC2','P4','C4','F4','F8','FC6','CP6','T8','P8'};
[found,flip] = ismember(channel_flip,cellstr(channel_labels));
assert(all(found),'Channel-flip order does not match channels.mat.');
right_idx = contains(channel_side,'right');
left_idx = contains(channel_side,'left');

%% Rebuild the bilateral EP-map collection
n_requested = numel(subject_ids);
maps = cell(n_requested*2,1);
map_subject = strings(n_requested*2,1);
included = false(n_requested,1);
for subject_idx = 1:n_requested
    subject_id = subject_ids(subject_idx);
    left_files = dir(fullfile(data_path,subject_id+'_left*.mat'));
    right_files = dir(fullfile(data_path,subject_id+'_right*.mat'));
    if numel(left_files)~=1 || numel(right_files)~=1
        warning('Skipping %s: expected one left and one right map.',subject_id);
        continue;
    end
    left = local_load_map(fullfile(left_files(1).folder,left_files(1).name));
    right = local_load_map(fullfile(right_files(1).folder,right_files(1).name));
    assert(size(left,1)==n_channels && size(right,1)==n_channels, ...
        'Unexpected channel count for %s.',subject_id);
    L = left(:,epoch_idx); L(right_idx,:) = NaN;
    R = right(:,epoch_idx); R(left_idx,:) = NaN;
    maps{subject_idx} = R(flip,:);
    maps{n_requested+subject_idx} = L;
    map_subject([subject_idx,n_requested+subject_idx]) = subject_id;
    included(subject_idx) = true;
end

keep_entry = [included;included];
maps = maps(keep_entry);
map_subject = map_subject(keep_entry);
subject_ids = subject_ids(included);
n_subjects = numel(subject_ids);
assert(n_subjects>=2,'At least two subjects with bilateral maps are required.');

map_data = cat(3,maps{:});
subject_data = nan(n_channels,numel(time),n_subjects);
for subject_idx = 1:n_subjects
    subject_data(:,:,subject_idx) = mean( ...
        map_data(:,:,map_subject==subject_ids(subject_idx)),3,'omitnan');
end

baseline_mean = mean(subject_data(:,baseline_idx,:),2,'omitnan');
difference_data = subject_data-baseline_mean;

%% Parametric and permutation inference
[t_full,p_full,n_full] = local_ttest_zero(difference_data);
t_observed = t_full(:,test_idx);
p_parametric = p_full(:,test_idx);
n_valid = n_full(:,test_idx);
initial_mask = p_parametric<initial_alpha;

rng(random_seed,'twister');
test_difference = permute(difference_data(:,test_idx,:),[3,1,2]);
test_difference = reshape(test_difference,n_subjects,[]);
t_vector = t_observed(:)';
valid = isfinite(test_difference);
difference_zero = test_difference;
difference_zero(~valid) = 0;
n_per_test = sum(valid,1);
sum_squares = sum(difference_zero.^2,1);
exceed_count = zeros(1,size(test_difference,2));

for first_permutation = 1:permutation_batch_size:n_permutations
    this_batch = min(permutation_batch_size, ...
        n_permutations-first_permutation+1);
    signs = 2*(rand(this_batch,n_subjects)>=0.5)-1;
    t_permutation = local_t_from_sums(signs*difference_zero, ...
        sum_squares,n_per_test);
    exceed_count = exceed_count+sum(abs(t_permutation)>=abs(t_vector),1);
end

p_permutation = reshape((exceed_count+1)./(n_permutations+1), ...
    size(t_observed));
p_permutation(n_valid<2 | ~isfinite(t_observed)) = NaN;
q_fdr = local_bh_adjust(p_permutation);
combined_mask = initial_mask & q_fdr<=fdr_alpha;
significant_mask = local_minimum_run(combined_mask,minimum_samples);

%% Return full-epoch maps and write machine-readable results
t_map = nan(n_channels,numel(time));
p_parametric_map = nan(n_channels,numel(time));
p_permutation_map = nan(n_channels,numel(time));
q_fdr_map = nan(n_channels,numel(time));
initial_mask_map = false(n_channels,numel(time));
significant_mask_map = false(n_channels,numel(time));
t_map(:,test_idx) = t_observed;
p_parametric_map(:,test_idx) = p_parametric;
p_permutation_map(:,test_idx) = p_permutation;
q_fdr_map(:,test_idx) = q_fdr;
initial_mask_map(:,test_idx) = initial_mask;
significant_mask_map(:,test_idx) = significant_mask;

cluster_table = local_run_table(significant_mask,t_observed,q_fdr, ...
    time(test_idx),channel_labels,sampling_interval);
results = struct();
results.subjectIDs = subject_ids;
results.channelLabels = channel_labels;
results.time = time;
results.baselineIndices = baseline_idx;
results.testIndices = test_idx;
results.nValid = n_valid;
results.baselineMean = baseline_mean;
results.subjectData = subject_data;
results.grandAverage = mean(subject_data,3,'omitnan');
results.tMap = t_map;
results.pParametric = p_parametric_map;
results.pPermutation = p_permutation_map;
results.qFDR = q_fdr_map;
results.initialMask = initial_mask_map;
results.significantMask = significant_mask_map;
results.clusterTable = cluster_table;
results.settings = struct('epochWindow',epoch_window, ...
    'baselineWindow',baseline_window,'testWindow',test_window, ...
    'initialAlpha',initial_alpha,'fdrAlpha',fdr_alpha, ...
    'nPermutations',n_permutations,'minimumDuration',minimum_duration, ...
    'randomSeed',random_seed);
save(fullfile(results_path,'evoked_vs_baseline_statistics.mat'),'results','-v7.3');
writetable(cluster_table, ...
    fullfile(results_path,'evoked_vs_baseline_significant_runs.csv'));

fprintf('Significant bins: %d; temporal runs: %d.\n', ...
    nnz(significant_mask),height(cluster_table));

function map = local_load_map(filename)
loaded = load(filename);
assert(isfield(loaded,'m'),'File %s does not contain m.',filename);
map = double(loaded.m);
end

function [t_value,p_value,n_valid] = local_ttest_zero(data)
valid = isfinite(data); data(~valid) = 0;
n_valid = sum(valid,3);
sum_data = sum(data,3); sum_squares = sum(data.^2,3);
mean_data = sum_data./n_valid;
variance = (sum_squares-(sum_data.^2)./n_valid)./(n_valid-1);
variance = max(variance,0);
standard_error = sqrt(variance./n_valid);
t_value = mean_data./standard_error;
zero_over_zero = standard_error==0 & mean_data==0;
t_value(zero_over_zero) = 0; t_value(n_valid<2) = NaN;
degrees_freedom = n_valid-1;
p_value = nan(size(t_value));
testable = n_valid>=2 & ~isnan(t_value);
argument = degrees_freedom(testable)./(degrees_freedom(testable)+ ...
    t_value(testable).^2);
p_value(testable) = betainc(argument,degrees_freedom(testable)/2,0.5);
p_value(zero_over_zero) = 1;
end

function t_value = local_t_from_sums(signed_sums,sum_squares,n_valid)
mean_value = signed_sums./n_valid;
variance = (sum_squares-(signed_sums.^2)./n_valid)./(n_valid-1);
variance = max(variance,0);
standard_error = sqrt(variance./n_valid);
t_value = mean_value./standard_error;
t_value(standard_error==0 & mean_value==0) = 0;
t_value(:,n_valid<2) = NaN;
end

function q_value = local_bh_adjust(p_value)
q_value = nan(size(p_value));
finite_idx = find(isfinite(p_value));
if isempty(finite_idx), return; end
[p_sorted,sort_order] = sort(p_value(finite_idx));
n_tests = numel(p_sorted);
q_sorted = p_sorted.*n_tests./(1:n_tests)';
for idx = n_tests-1:-1:1
    q_sorted(idx) = min(q_sorted(idx),q_sorted(idx+1));
end
unsorted = nan(n_tests,1); unsorted(sort_order) = min(q_sorted,1);
q_value(finite_idx) = unsorted;
end

function filtered = local_minimum_run(mask,minimum_samples)
filtered = false(size(mask));
for channel_idx = 1:size(mask,1)
    transitions = diff([false,mask(channel_idx,:),false]);
    starts = find(transitions==1); stops = find(transitions==-1)-1;
    for run_idx = find((stops-starts+1)>=minimum_samples)
        filtered(channel_idx,starts(run_idx):stops(run_idx)) = true;
    end
end
end

function output = local_run_table(mask,t_value,q_value,time,labels,dt)
channel_col = strings(0,1); start_col = zeros(0,1); end_col = zeros(0,1);
duration_col = zeros(0,1); peak_time_col = zeros(0,1);
peak_t_col = zeros(0,1); minimum_q_col = zeros(0,1); direction_col = strings(0,1);
for channel_idx = 1:size(mask,1)
    transitions = diff([false,mask(channel_idx,:),false]);
    starts = find(transitions==1); stops = find(transitions==-1)-1;
    for run_idx = 1:numel(starts)
        idx = starts(run_idx):stops(run_idx);
        [~,relative_peak] = max(abs(t_value(channel_idx,idx)));
        peak_idx = idx(relative_peak); peak_t = t_value(channel_idx,peak_idx);
        channel_col(end+1,1) = labels(channel_idx); %#ok<AGROW>
        start_col(end+1,1) = time(idx(1))*1000; %#ok<AGROW>
        end_col(end+1,1) = time(idx(end))*1000; %#ok<AGROW>
        duration_col(end+1,1) = numel(idx)*dt*1000; %#ok<AGROW>
        peak_time_col(end+1,1) = time(peak_idx)*1000; %#ok<AGROW>
        peak_t_col(end+1,1) = peak_t; %#ok<AGROW>
        minimum_q_col(end+1,1) = min(q_value(channel_idx,idx)); %#ok<AGROW>
        if peak_t>=0, direction_col(end+1,1) = "Positive"; ...
        else, direction_col(end+1,1) = "Negative"; end %#ok<AGROW>
    end
end
output = table(channel_col,start_col,end_col,duration_col,peak_time_col, ...
    peak_t_col,minimum_q_col,direction_col,'VariableNames', ...
    {'Channel','Start_ms','End_ms','Duration_ms','PeakTime_ms','Peak_t', ...
    'Minimum_FDR_q','Direction'});
end
