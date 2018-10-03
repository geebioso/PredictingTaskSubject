%% Script to plot a matrix of standard deviation in BOLD activation
function [] = plot_BV_heatmap (whs, truncate, roi_sort_by, fontsize, labels_on, ...
    normalize, grid_on, remove_max, isHPC) 

%% INPUT: 
%     str roi_sort_by: how to sort rois {'Lobe' or 'Hemishpere'}
%     int fontsize
%     bool labels_on: should we have individual row labels
%     int normalize: how should we normalize? % 0 is no normalization, 1 is by subject, 2 is by ROI, 3 is globally
%     bool grid_on = 1; % do we want to plot a grid on figures? 
%     bool remove_max = 0; % do we want to remove the max variance value? 

if whs==0; bad_roi=253; else; bad_roi=[]; end % ROI with low outlying variance (near 0)

if normalize
    cmin = -4;
    cmax = 4;
    colormapping = 'scaled';
else
    cmin = 13.6211;
    cmax = 199.1290;
    colormapping = 'direct';
end

plot_19_subjs_subj_mean = 0;
plot_19_subjs_subj = 1;
plot_19_subjs_task = 1;
plot_ordered = 1;
plot_sanity_check = 0;
K = 9; 

dorepresents      = [ 1 ];

% whrep==1: ROI StdDev
% whrep==2: ROI means


featurelabels = { ...
    'StdDev ROIs' , ...
    'Means  ROIs' , ...
    'ROI FV scaled',...
    'ROI FV',...
    'MSSD', ...
    'SQRT_MSS', ...
    'Corrs  ROIs' , ...
    'Covs   ROIS', ...
    'Covs+Vars ROIs', ...
    'Covs   ROIS scaled'
    };

%% Load the task names
[wd, rd] = set_directories(isHPC); % working directory and results directory

if ~exist( 'tasks' , 'var' )
    load( fullfile(rd, 'ICAresults','tasklist') );
end

task_labels = {'Emotional Pictures', 'Emotional Faces', ...
    'Encoding', 'Go/No-go', 'Rest', 'Retrieval', 'Monetary Incentive', 'Working Memory', 'Theory of Mind'}; 

%% Load Data 
try
    filenm = fullfile( rd, 'features', sprintf('whs%d_truncate%d.mat', whs, truncate));
    fprintf('Loading %s\n', filenm);
    load( filenm );  
%                     NT: 9
%                     NR: 269
%                  combs: [1827×3 double]
%                     NC: 1827
%     subsinbothsessions: [19×1 double]
%               features: [1×1 struct]

% combs columns are [subject, task, session] 
catch
    error('Features are not computed for dataset %d', whs);
end

% compute extra features
roi_vars_scale = features.stds.^2;
features.vars_scale = (roi_vars_scale - mean(roi_vars_scale(:)))/std(roi_vars_scale(:));
features.covs_scale = (features.covs - mean(features.covs(:)))/(std(features.covs(:)));

% compute number of subjects in both sessions 
unique_subs = unique(combs(:,1));
NS = length(unique_subs); 
nspersubj = accumarray( combs(:,1) , combs(:,3) , [ NS 1 ] , @mean );
subsinbothsessions = find( nspersubj == 1.5 );
NSboth = length(subsinbothsessions);

% remove bad_roi 
features.stds(:,bad_roi) = []; 

%% Normalization 
if min(features.stds(:)) > 0
    fprintf('No Rois have variance of 0\nmin(features.stds(:))=%2.2f', min(features.stds(:)));
end

if normalize
    session1_idx = ismember(combs(:,3),1);
    
    
    if normalize == 1
        
        for i = 1:NS
            s = unique_subs(i);
            sub_idx = ismember(combs(:,1),s);
            
            % session 1
            idx = and(sub_idx, session1_idx);
            features.stds_norm(idx,:) = bsxfun( @minus, features.stds(idx,:),mean(features.stds(idx,:)));
            features.stds_norm(idx,:) = bsxfun( @rdivide, features.stds_norm(idx,:), std(features.stds(idx,:)));
            
            features.means_norm(idx,:) = bsxfun( @minus, features.means(idx,:),mean(features.means(idx,:)));
            features.means_norm(idx,:) = bsxfun( @rdivide, features.means_norm(idx,:), std(features.means(idx,:)));
            
            % session 2
            idx = and(sub_idx, ~session1_idx);
            features.stds_norm(idx,:) = bsxfun( @minus, features.stds(idx,:),mean(features.stds(idx,:)));
            features.stds_norm(idx,:) = bsxfun( @rdivide, features.stds_norm(idx,:), std(features.stds(idx,:)));
            
            features.means_norm(idx,:) = bsxfun( @minus, features.means(idx,:),mean(features.means(idx,:)));
            features.means_norm(idx,:) = bsxfun( @rdivide, features.means_norm(idx,:), std(features.means(idx,:)));
            
            
        end
    elseif normalize == 2
        for r = 1:NR
            % sessoin 1 
            idx = session1_idx;
            features.stds_norm(idx, r) = bsxfun( @minus, features.stds(idx, r),mean(features.stds(idx, r)));
            features.stds_norm(idx, r) = bsxfun( @rdivide, features.stds_norm(idx, r), std(features.stds(idx, r)));
            
            features.means_norm(idx, r) = bsxfun( @minus, features.means(idx, r),mean(features.means(idx, r)));
            features.means_norm(idx, r) = bsxfun( @rdivide, features.means_norm(idx, r), std(features.means(idx, r)));
            
            % session 2 
            idx = ~session1_idx;
            features.stds_norm(idx, r) = bsxfun( @minus, features.stds(idx, r),mean(features.stds(idx, r)));
            features.stds_norm(idx, r) = bsxfun( @rdivide, features.stds(idx, r), std(features.stds(idx, r)));
            
            features.means_norm(idx, r) = bsxfun( @minus, features.means(idx, r),mean(features.means(idx, r)));
            features.means_norm(idx, r) = bsxfun( @rdivide, features.means_norm(idx, r), std(features.means(idx, r)));
        end
    elseif normalize == 3
        
        % session 1
        idx = session1_idx;
        features.stds_norm(idx,:) = bsxfun( @minus, features.stds(idx,:),mean(features.stds(idx,:)));
        features.stds_norm(idx,:) = bsxfun( @rdivide, features.stds_norm(idx,:), std(features.stds(idx,:)));
        
        features.means_norm(idx,:) = bsxfun( @minus, features.means(idx,:),mean(features.means(idx,:)));
        features.means_norm(idx,:) = bsxfun( @rdivide, features.means_norm(idx,:), std(features.means(idx,:)));
        
        % session 2
        idx = ~session1_idx;
        features.stds_norm(idx,:) = bsxfun( @minus, features.stds(idx,:),mean(features.stds(idx,:)));
        features.stds_norm(idx,:) = bsxfun( @rdivide, features.stds_norm(idx,:), std(features.stds(idx,:)));
        
        features.means_norm(idx,:) = bsxfun( @minus, features.means(idx,:),mean(features.means(idx,:)));
        features.means_norm(idx,:) = bsxfun( @rdivide, features.means_norm(idx,:), std(features.means(idx,:)));
    end
end

%% Get new roi indices into old rois
[new_rois_old_idx] = get_new_roi_index_into_old_rois(rd);

%% Compute ROI Ordering
filename = fullfile(rd, '..','FMRI','BrainVisualization','Network Classification','ROI_299_list.csv');
T = readtable(filename);
T = table(T.Index, T.aal_brodmannLabel, T.notes, 'VariableNames', {'ROI', 'label', 'notes'});
T = T(new_rois_old_idx,:);
regionnames = T.label;

[ lobes, lobenames, primaryareas, primarylobes, regionareas, regionlobes, regionareas2, regionlobes2, lobenames2, primaryareas2, primarylobes2 ] = getLobeSets( regionnames );

% get hemisphere
left = ~cellfun(@isempty, strfind(regionnames, '_L_'));
right = ~cellfun(@isempty, strfind(regionnames, '_R_'));
center = and(~left, ~right);

hemisphere = cell(length(left), 1);
hemisphere(left) = {'left'};
hemisphere(right) = {'right'};
hemisphere(center) = {'center'};

% sort
regionlobes(bad_roi) = []; 
hemisphere(bad_roi) = []; 
[~, lobe_idx] = sort(regionlobes);
[~, hemi_idx] = sort(hemisphere);

if strcmp(roi_sort_by, 'Lobe')
    sort_idx = lobe_idx;
elseif strcmp(roi_sort_by, 'Hemisphere')
    sort_idx = hemi_idx;
end

[sorted_tasks, task_sort_idx] = sort(combs(:,2));
%% Plots

if plot_ordered
    % Get ROI sort labels etc.
    
    if strcmp(roi_sort_by, 'Lobe')
        [unique_lobes,first_lobe_occurances,~] = unique(regionlobes(sort_idx), 'first'); % get first occurance of each unique element
    elseif strcmp(roi_sort_by, 'Hemisphere')
        [unique_lobes,first_lobe_occurances,~] = unique(hemisphere(sort_idx), 'first'); % get first occurance of each unique element
    end
    h = figure(); clf;
    h.Position =[92 0 1638 665];
    h1 = subplot(1,2,1);
    
    [unique_subjs,first_subj_occurances,~] = unique(combs(:,1), 'first'); % get first occurance of each unique element
    
    
    image(features.stds(:,sort_idx), 'CDataMapping', colormapping); colorbar; hold on;
    h1.YTick = first_subj_occurances(1:5:end) + K/2;
    h1.YTickLabels = 1:5:length(unique_subjs);
    
    diffs = diff([first_lobe_occurances; NR])/2;
    xticks = diffs + first_lobe_occurances;
    h1.XTick = xticks; % first_lobe_occurances;
    h1.XTickLabels = unique_lobes;
    h1.XTickLabelRotation = 45;
    h1.FontSize = fontsize;
    
    %title(sprintf('FV Ordered by Subject'));
    title('A');
    if labels_on
        xlabel(sprintf('ROI (ordered by %s)', lower(roi_sort_by)));
        ylabel('Scan (ordered by subject)');
    end
    
    % plot gridlines
    if grid_on
        for i = 1:length(first_lobe_occurances)
            x = first_lobe_occurances(i);
            plot([x,x], [1 size(features.stds,1)], 'k', 'LineWidth', 1); hold on;
        end
    end
    
    [sorted_tasks, task_sort_idx] = sort(combs(:,2));
    
    h2 = subplot(1,2,2);
    
    [unique_tasks,first_task_occurances,~] = unique(sorted_tasks, 'first'); % get first occurance of each unique element
    image(features.stds(task_sort_idx, sort_idx), 'CDataMapping', colormapping); colorbar; hold on;
    h2.YTick = first_task_occurances + NSboth/2;
    h2.YTickLabels = task_labels(unique_tasks);
    
    diffs = diff([first_lobe_occurances; NR])/2;
    xticks = diffs + first_lobe_occurances;
    h2.XTick = xticks; % first_lobe_occurances;
    h2.XTickLabels = unique_lobes;
    h2.XTickLabelRotation = 45;
    h2.FontSize = fontsize;
    
    %title(sprintf('FV Ordered by Task'));
    title('B');
    if labels_on
        xlabel(sprintf('ROI (ordered by %s)', lower(roi_sort_by)));
        ylabel('Scan (ordered by subject)');
    end
    
    % plot gridlines
    if grid_on
        for i = 1:length(first_lobe_occurances)
            x = first_lobe_occurances(i);
            plot([x,x], [1 length(task_sort_idx)], 'k', 'LineWidth', 1); hold on;
        end
        
        for i = 1:length(first_task_occurances)
            y = first_task_occurances(i);
            plot([1 length(task_sort_idx)], [y,y], 'k', 'LineWidth', 1); hold on;
        end
    end
    
    filename = fullfile(rd, 'images', sprintf('functional_variance_subject_vs_task_heatmap_%s_normalize=%d_truncate%d', ...
        roi_sort_by, normalize, truncate));
    print(filename, '-depsc');
    
end
%% Plot For only 19 subjects in both sessions
if plot_19_subjs_subj
    
    bothsess = ismember(combs(:,1), subsinbothsessions);
    session1 = and(combs(:,3) == 1, bothsess);
    session2 = and(combs(:,3) == 2, bothsess);
    
    features.stds_sess1 = features.stds(session1, :);
    features.stds_sess2 = features.stds(session2, :);
    
    h = figure(); clf;
    h.Position =[92 0 1638 665];
    
    % session 1
    h1 = subplot(1,2,1);
    
    [sess1_subjs, sort_idx] = sort(combs(session1, 1));
    [unique_subjs,first_subj_occurances,~] = unique(sess1_subjs, 'first'); % get first occurance of each unique element
    
    image(features.stds_sess1(sort_idx, lobe_idx), 'CDataMapping', colormapping); colorbar;
    if normalize == 0
        caxis([cmin, cmax]);
    end
    hold on;
    
    h1.YTick = first_subj_occurances + K/2;
    h1.YTickLabels = 1:NSboth; % subsinbothsessions;
    h1.TickLength = [0 0]; 
    
    diffs = diff([first_lobe_occurances; NR])/2;
    xticks = diffs + first_lobe_occurances - 0.5;
    h1.XTick = xticks; % first_lobe_occurances;
    h1.XTickLabels = unique_lobes;
    h1.XTickLabelRotation = 45;
    h1.FontSize = fontsize;
    
    %title(sprintf('Session 1'));
    title('A');
    if labels_on
        xlabel(sprintf('ROI (ordered by %s)', lower(roi_sort_by)));
        ylabel('Scan (ordered by subject)');
    end
    
    % plot gridlines
    if grid_on
        for i = 1:length(first_lobe_occurances)
            x = first_lobe_occurances(i) - 0.5;
            plot([x,x], [1 length(task_sort_idx)], 'k', 'LineWidth', 1); hold on;
        end
        
        for i = 1:length(first_subj_occurances)
            y = first_subj_occurances(i) - 0.5;
            plot([1 length(task_sort_idx)], [y,y], 'k', 'LineWidth', 1); hold on;
        end
    end
    % Session 2
    h2 = subplot(1,2,2);
    
    [sess2_subjs, sort_idx] = sort(combs(session2, 1));
    [unique_subjs,first_subj_occurances,~] = unique(sess2_subjs, 'first'); % get first occurance of each unique element
    
    if remove_max
        [~,ii] = max(features.stds_sess2(:));
        features.stds_sess2(ii) = mean(features.stds_sess2(:));
        image(features.stds_sess2(sort_idx, :), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    else
        image(features.stds_sess2(sort_idx,lobe_idx), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    end
    
    h2.YTick = first_subj_occurances + K/2;
    h2.YTickLabels = 1:NSboth; % subsinbothsessions;
    h2.TickLength = [0 0]; 
    
    diffs = diff([first_lobe_occurances; NR])/2;
    xticks = diffs + first_lobe_occurances;
    h2.XTick = xticks - 0.5; % first_lobe_occurances;
    h2.XTickLabels = unique_lobes;
    h2.XTickLabelRotation = 45;
    h2.FontSize = fontsize;
    
    %title(sprintf('Session 2'));
    title('B');
    if labels_on
        xlabel(sprintf('ROI (ordered by %s)', lower(roi_sort_by)));
        ylabel('Scan (ordered by subject)');
    end
    
    % plot gridlines
    if grid_on
        for i = 1:length(first_lobe_occurances)
            x = first_lobe_occurances(i) - 0.5;
            plot([x,x], [1 length(task_sort_idx)], 'k', 'LineWidth', 1); hold on;
        end
        
        for i = 1:length(first_subj_occurances)
            y = first_subj_occurances(i)- 0.5;
            plot([1 length(task_sort_idx)], [y,y], 'k', 'LineWidth', 1); hold on;
        end
    end
    filename = fullfile(rd, 'images', sprintf('FV_19_subjs_heatmap_subjorder_%s_normalize=%d_whs%d_truncate%d',...
        roi_sort_by, normalize, whs, truncate));
    print(filename, '-depsc');
    
    %% Subject statistical testing 
    BV = features.stds_sess2(sort_idx,lobe_idx); 
    subj_idx = repmat(1:NSboth, NT, 1);
    subj_idx = subj_idx(:); 
    
    % manova1(BV, subj_idx)
    
    [coeff,score,latent,tsquared,explained,mu] = pca(BV); 
    i = find( cumsum(explained) > 99.5); 
    
    t = score(:, 1:i(1))*coeff(:, 1:i(1))' + repmat(mu,length(subj_idx),1); 
    [d,p,stats] = manova1(score(:,1:i(1)), subj_idx)
    
    %% Task Statistical Testing 
    [sess2_subjs, sort_idx] = sort(combs(session2, 2));
    taskBV = features.stds_sess2(sort_idx,lobe_idx); 
    task_idx = repmat(1:NT, NSboth, 1);
    task_idx = task_idx(:); 
    
    [coeff,score,latent,tsquared,explained,mu] = pca(taskBV); 
    i = find( cumsum(explained) > 100); 
    
    t = score(:, 1:i(1))*coeff(:, 1:i(1))' + repmat(mu,length(task_idx),1); 
    [d,p,stats] = manova1(score(:,1:i(1)), task_idx)
    
    idx = kmeans(taskBV,d)
    
end


%% Plot For only 19 subjects in both sessions (TASK ordering) 
if plot_19_subjs_task
    
    
    bothsess = ismember(combs(:,1), subsinbothsessions);
    session1 = and(combs(:,3) == 1, bothsess);
    session2 = and(combs(:,3) == 2, bothsess);
    
    features.stds_sess1 = features.stds(session1, :);
    features.stds_sess2 = features.stds(session2, :);
    
    taskid = combs(:,2);
    [sess1_tasks, sort_idx] = sort(taskid(session1));
    
    if strcmp(roi_sort_by, 'Lobe')
        [unique_lobes,first_lobe_occurances,~] = unique(regionlobes(lobe_idx), 'first'); % get first occurance of each unique element
    elseif strcmp(roi_sort_by, 'Hemisphere')
        [unique_lobes,first_lobe_occurances,~] = unique(hemisphere(hem_idx), 'first'); % get first occurance of each unique element
    end
    
    h = figure(); clf;
    h.Position =[92 0 1638 665];
    
    % session 1
    h1 = subplot(1,2,1);
    
    [unique_tasks,first_task_occurances,~] = unique(sess1_tasks, 'first'); % get first occurance of each unique element
    image(features.stds_sess1(sort_idx,lobe_idx), 'CDataMapping', colormapping); colorbar;
    
    K = length(unique_tasks);
    
    if normalize == 0
        caxis([cmin, cmax]);
    end
    hold on;
    
    h1.YTick = first_task_occurances + NSboth/2;
    h1.YTickLabels = task_labels; % subsinbothsessions;
    h1.TickLength = [0 0]; 
    
    diffs = diff([first_lobe_occurances; NR])/2;
    xticks = diffs + first_lobe_occurances;
    h1.XTick = xticks - 0.5; % first_lobe_occurances;
    h1.XTickLabels = unique_lobes;
    h1.XTickLabelRotation = 45;
    h1.FontSize = fontsize;
    
    %title(sprintf('Session 1'));
    title('A');
    if labels_on
        xlabel(sprintf('ROI (ordered by %s)', lower(roi_sort_by)));
        ylabel('Scan (ordered by subject)');
    end
    
    % plot gridlines
    if grid_on
        for i = 1:length(first_lobe_occurances)
            x = first_lobe_occurances(i) - 0.5;
            plot([x,x], [1 length(sort_idx)], 'k', 'LineWidth', 1); hold on;
        end
        
        for i = 1:length(first_task_occurances)
            y = first_task_occurances(i) - 0.5;
            plot([1 length(lobe_idx)], [y,y], 'k', 'LineWidth', 1); hold on;
        end
    end
    
    % Session 2
    h2 = subplot(1,2,2);
    
    [sess2_tasks, sort_idx] = sort(taskid(session2));
    [unique_tasks,first_task_occurances,~] = unique(sess2_tasks, 'first'); % get first occurance of each unique element
    
    
    
    if strcmp(roi_sort_by, 'Lobe')
        [unique_lobes,first_lobe_occurances,~] = unique(regionlobes(lobe_idx), 'first'); % get first occurance of each unique element
    elseif strcmp(roi_sort_by, 'Hemisphere')
        [unique_lobes,first_lobe_occurances,~] = unique(hemisphere(hem_idx), 'first'); % get first occurance of each unique element
    end
    
    if remove_max
        [~,ii] = max(features.stds_sess2(:));
        features.stds_sess2(ii) = mean(features.stds_sess2(:));
        image(features.stds_sess2(sort_idx,lobe_idx), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    else
        image(features.stds_sess2(sort_idx,lobe_idx), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    end
    
    h2.YTick = first_task_occurances + NSboth/2; 
    h2.YTickLabels = task_labels; % subsinbothsessions;
    h2.TickLength = [0 0]; 
    
    diffs = diff([first_lobe_occurances; NR])/2;
    xticks = diffs + first_lobe_occurances;
    h2.XTick = xticks - 0.5; % first_lobe_occurances;
    h2.XTickLabels = unique_lobes;
    h2.XTickLabelRotation = 45;
    h2.FontSize = fontsize;
    
    %title(sprintf('Session 2'));
    title('B');
    if labels_on
        xlabel(sprintf('ROI (ordered by %s)', lower(roi_sort_by)));
        ylabel('Scan (ordered by subject)');
    end
    
    % plot gridlines
    if grid_on
        for i = 1:length(first_lobe_occurances)
            x = first_lobe_occurances(i) - 0.5;
            plot([x,x], [1 length(sort_idx)], 'k', 'LineWidth', 1); hold on;
        end
        
        for i = 1:length(first_task_occurances)
            y = first_task_occurances(i) - 0.5;
            plot([1 length(lobe_idx)], [y,y], 'k', 'LineWidth', 1); hold on;
        end
    end
    
    filename = fullfile(rd, 'images', sprintf('FV_19_subjs_heatmap_taskorder_%s_normalize=%d_whs%d_truncate%d', ...
        roi_sort_by, normalize, whs, truncate));
    print(filename, '-depsc');
    
end


%% Plot FV ordered by highest variance for each subject

if plot_sanity_check
    
    % sort each subject by highest var
    features.stds_sort = zeros(size(features.stds));
    for s = 1:size(features.stds,1)
        features.stds_sort(s,:) = sort(features.stds(s,:));
    end
    
    % make some random data
    features.stds_rand_sort = zeros(size(features.stds));
    for s = 1:size(features.stds,1)
        indices = datasample(1:size(features.stds,1), NR); % sample a subject to take value from for each region
        
        rand_subj = zeros(1,NR);
        for r = 1:NR
            rand_subj(r) = features.stds(indices(r), r); % choose ROI from subject we sampled
        end
        features.stds_rand_sort(s,:) = rand_subj;
    end
    
    for s = 1:size(features.stds,1)
        features.stds_rand_sort(s,:) = sort(features.stds_rand_sort(s,:));
    end
    
    h = figure(); clf;
    h.Position =[92 0 600 600];
    
    % session 1
    subplot(1,1,1);
    
    h1 = plot(features.stds_sort(combs(:,3) == 1, :)', 'Color', [0 0 1 0.5]); hold on;
    
    title(sprintf('Session 1 Sorted FVs'));
    xlabel('ROI (sorted by FV)');
    ylabel('FV');
    
    %h2 = subplot(1,2,2);
    
    h2 = plot(features.stds_rand_sort(combs(:,3) == 1, :)', 'Color', [1 0 0 0.5]);
    
    title(sprintf('Session 1 Sorted FVs'));
    xlabel('ROI (sorted by FV)');
    ylabel('FV');
    
    h = legend([h1(1), h2(1)], {'sorted subjects', 'sorted random'});
    
    xlim([1,NR]);
    filename = fullfile(rd, 'images', sprintf('functional_variance_sorted_vs_random_%s_truncate%d', ...
        roi_sort_by, truncate));
    print(filename, '-depsc');
end


%% Means
if plot_19_subjs_subj_mean
    
    colormapping = 'scaled';
    
    bothsess = ismember(combs(:,1), subsinbothsessions);
    session1 = and(combs(:,3) == 1, bothsess);
    session2 = and(combs(:,3) == 2, bothsess);
    
    features.means_norm_sess1 = features.means_norm(session1, :);
    features.means_norm_sess2 = features.means_norm(session2, :);
    
    h = figure(); clf;
    h.Position =[92 0 1638 665];
    
    % session 1
    h1 = subplot(1,2,1);
    
    [sess1_subjs, sort_idx] = sort(combs(session1, 1));
    [unique_subjs,first_subj_occurances,~] = unique(sess1_subjs, 'first'); % get first occurance of each unique element
    
    image(features.means_norm_sess1(:,sort_idx), 'CDataMapping', colormapping); colorbar;
    if normalize == 0
        caxis([cmin, cmax]);
    end
    hold on;
    
    h1.YTick = first_subj_occurances - 0.5;
    h1.YTickLabels = 1:NSboth; % subsinbothsessions;
    
    diffs = diff([first_lobe_occurances; NR])/2;
    xticks = diffs + first_lobe_occurances;
    h1.XTick = xticks - 0.5; % first_lobe_occurances;
    h1.XTickLabels = unique_lobes;
    h1.XTickLabelRotation = 45;
    h1.FontSize = fontsize;
    
    %title(sprintf('Session 1'));
    title('A');
    if labels_on
        xlabel(sprintf('ROI (ordered by %s)', lower(roi_sort_by)));
        ylabel('Scan (ordered by subject)');
    end
    
    % plot gridlines
    if grid_on
        for i = 1:length(first_lobe_occurances)
            x = first_lobe_occurances(i) - 0.5;
            plot([x,x], [1 length(task_sort_idx)], 'k', 'LineWidth', 1); hold on;
        end
        
        for i = 1:length(first_subj_occurances)
            y = first_subj_occurances(i) - 0.5;
            plot([1 length(task_sort_idx)], [y,y], 'k', 'LineWidth', 1); hold on;
        end
    end
    % Session 2
    h2 = subplot(1,2,2);
    
    [sess2_subjs, sort_idx] = sort(combs(session2, 1));
    [unique_subjs,first_subj_occurances,~] = unique(sess2_subjs, 'first'); % get first occurance of each unique element
    
    if remove_max
        [~,ii] = max(features.means_norm_sess2(:));
        features.means_norm_sess2(ii) = mean(features.means_norm_sess2(:));
        image(features.means_norm_sess2(:,sort_idx), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    else
        image(features.means_norm_sess2(:,sort_idx), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    end
    
    h2.YTick = first_subj_occurances - 0.5;
    h2.YTickLabels = 1:NSboth; % subsinbothsessions;
    
    diffs = diff([first_lobe_occurances; NR])/2;
    xticks = diffs + first_lobe_occurances;
    h2.XTick = xticks - 0.5; % first_lobe_occurances;
    h2.XTickLabels = unique_lobes;
    h2.XTickLabelRotation = 45;
    h2.FontSize = fontsize;
    
    %title(sprintf('Session 2'));
    title('B');
    if labels_on
        xlabel(sprintf('ROI (ordered by %s)', lower(roi_sort_by)));
        ylabel('Scan (ordered by subject)');
    end
    
    % plot gridlines
    if grid_on
        for i = 1:length(first_lobe_occurances)
            x = first_lobe_occurances(i) - 0.5;
            plot([x,x], [1 length(task_sort_idx)], 'k', 'LineWidth', 1); hold on;
        end
        
        for i = 1:length(first_subj_occurances)
            y = first_subj_occurances(i)- 0.5;
            plot([1 length(task_sort_idx)], [y,y], 'k', 'LineWidth', 1); hold on;
        end
    end
    filename = fullfile(rd, 'images', sprintf('M_19_subjs_heatmap_subjorder_%s_normalize=%d_truncate%d', ...
        roi_sort_by, normalize, truncate));
    print(filename, '-dpng');
    
end


%% Get Mean FD and sort 





mFD1 = load(fullfile(rd, '..', 'FMRI', 'restingstatedata', 'meanFD_174subj.mat')); 
mFD2 = load(fullfile(rd, '..', 'FMRI', 'restingstatedata', 'batch2_meanFD.mat'));


%% Load Head Motion
load( fullfile( rd, '..', 'FMRI', 'restingstatedata', 'meanFD_174subj.mat'));
analysis_subjs = subjs;
% x = load( fullfile( '..', 'FMRI', 'restingstatedata', 'batch2_meanFD.mat'));
% meanFD_2 = x.meanFD;
load( fullfile( rd, '..', 'FMRI', 'restingstatedata', 'motion_203Subj.mat'));
motion = motion( ismember( subjs, analysis_subjs), :);

[NS, T] = size(motion);







