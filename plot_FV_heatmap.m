%% Script to plot a matrix of standard deviation in BOLD activation

addpath(fullfile('..','FMRI','BrainVisualization','Network Classification'));
whs = 8;     % which data set? 8 is the most recent with edge voxels removed
remove_outlier = 0; % Remove outlying variance from session2 data
roi_sort_by = 'Lobe'; % 'Lobe'; % 'Lobe' or 'Hemishpere'


fontsize = 20;
labels_on = 0;
normalize = 0; % 0 is no normalization, 1 is by subject, 2 is by ROI, 3 is globally
grid_on = 1;

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


if (whs==8)
    roifile    = ['ICAresults' filesep 'tc_csf_wm_motion_out_globalMask5000_203subj_conv'];
    icfile     = '';
end

dorepresents      = [ 1 ];

% whrep==1: ROI StdDev
% whrep==2: ROI means


featurelabels = { ...
    'StdDev ROIs' , ...
    'Means  ROIs'
    };

%% Load the task names
if ~exist( 'tasks' , 'var' )
    load( fullfile('ICAresults','tasklist') );
end

task_labels = {'Emotional Pictures', 'Emotional Faces', ...
    'Encoding', 'Go/No-go', 'Rest', 'Retrieval', 'Monetary Incentive', 'Working Memory', 'Theory of Mind'}; 

%% Load the original ROI timecourse data
if ~exist( 'D' , 'var' )
    load( roifile );
end

%% Constants
NS = max( D.SubjectIndex2 ); % number of subjects

NT = length( tasks ); % number of tasks
NR = size( D.ROI , 2 ); % number of regions
%% Get combinations of Task x Subject x Session
[ combs , ~ , cindex ] = unique( [ D.SubjectIndex2 D.TaskIndex D.SessionIndex ], 'rows' );
NC = size( combs , 1 );

% Find subjects that are in session 1 and session 2
nspersubj = accumarray( combs(:,1) , combs(:,3) , [ NS 1 ] , @mean );
subsinbothsessions = find( nspersubj == 1.5 );
NSboth = length(subsinbothsessions);
%% Calculate Correlations for ROI across Task x Subject x Session
%if ~exist( 'roi_corrs' , 'var' )
ROI = D.ROI;

if ~exist( 'roi_means' , 'var' )
    % Cap the FFT at 120
    L = 120;
    L2 = L/2 + 1;
    
    NRP = NR * (NR-1) / 2; % off diagonal entries of FC matrix
    roi_covs  = zeros( NC , NRP );
    roi_corrs = zeros( NC , NRP );
    roi_means = zeros( NC , NR  );
    roi_stds  = zeros( NC , NR  );
    roi_fft   = zeros( NC , L2 * NR );
    
    warning off;
    triidx = ~triu(true(NR));
    for j=1:NC
        if mod(j,100) == 0
            fprintf( 'Working on combination %d of %d\n' , j , NC );
        end
        wh = find( cindex == j );
        ROI_NOW = ROI( wh , : );
        
        roi_means( j , : ) = mean( ROI_NOW , 1 );
        roi_stds( j , : )   = std( ROI_NOW , [] , 1 );
        rhos  = 1 - pdist( ROI_NOW' , 'correlation' );
        rhos2 = cov(ROI_NOW);
        rhos2 = rhos2(triidx)';
        roi_corrs( j , : ) = rhos;
        roi_covs( j, : ) = rhos2;
        
        
        yfft = fft( ROI_NOW( 1:L , : ) ,[],1)';
        P2 = abs( yfft/L );
        P1 = P2(:,1:L/2+1);
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        %roi_fft( j , : ) = mean( P1 , 1 );
        roi_fft( j , : ) = P1(:);
        
    end
    warning on;
    
    % Remove NaNs
    roi_corrs( isnan( roi_corrs )) = 0;
    roi_covs( isnan( roi_covs )) = 0;
    
    roi_stds_old = roi_stds;
    roi_means_old = roi_means;
    
    save( 'Results/roi_stds.mat', 'roi_stds'); 
end

if min(roi_stds(:)) > 0
    fprintf('No Rois have variance of 0\nmin(roi_stds(:))=%2.2f', min(roi_stds(:)));
end

if normalize
    unique_subs = unique(combs(:,1));
    session1_idx = ismember(combs(:,3),1);
    if normalize == 1
        
        for i = 1:NS
            s = unique_subs(i);
            sub_idx = ismember(combs(:,1),s);
            
            % session 1
            idx = and(sub_idx, session1_idx);
            roi_stds(idx,:) = bsxfun( @minus, roi_stds_old(idx,:),mean(roi_stds_old(idx,:)));
            roi_stds(idx,:) = bsxfun( @rdivide, roi_stds(idx,:), std(roi_stds_old(idx,:)));
            
            roi_means(idx,:) = bsxfun( @minus, roi_means_old(idx,:),mean(roi_means_old(idx,:)));
            roi_means(idx,:) = bsxfun( @rdivide, roi_means(idx,:), std(roi_means_old(idx,:)));
            
            % session 2
            idx = and(sub_idx, ~session1_idx);
            roi_stds(idx,:) = bsxfun( @minus, roi_stds_old(idx,:),mean(roi_stds_old(idx,:)));
            roi_stds(idx,:) = bsxfun( @rdivide, roi_stds(idx,:), std(roi_stds_old(idx,:)));
            
            roi_means(idx,:) = bsxfun( @minus, roi_means_old(idx,:),mean(roi_means_old(idx,:)));
            roi_means(idx,:) = bsxfun( @rdivide, roi_means(idx,:), std(roi_means_old(idx,:)));
            
            
        end
    elseif normalize == 2
        for r = 1:NR
            % sessoin 1 
            idx = session1_idx;
            roi_stds(idx, r) = bsxfun( @minus, roi_stds_old(idx, r),mean(roi_stds_old(idx, r)));
            roi_stds(idx, r) = bsxfun( @rdivide, roi_stds(idx, r), std(roi_stds_old(idx, r)));
            
            roi_means(idx, r) = bsxfun( @minus, roi_means_old(idx, r),mean(roi_means_old(idx, r)));
            roi_means(idx, r) = bsxfun( @rdivide, roi_means(idx, r), std(roi_means_old(idx, r)));
            
            % session 2 
            idx = ~session1_idx;
            roi_stds(idx, r) = bsxfun( @minus, roi_stds_old(idx, r),mean(roi_stds_old(idx, r)));
            roi_stds(idx, r) = bsxfun( @rdivide, roi_stds(idx, r), std(roi_stds_old(idx, r)));
            
            roi_means(idx, r) = bsxfun( @minus, roi_means_old(idx, r),mean(roi_means_old(idx, r)));
            roi_means(idx, r) = bsxfun( @rdivide, roi_means(idx, r), std(roi_means_old(idx, r)));
        end
    elseif normalize == 3
        
        % session 1
        idx = session1_idx;
        roi_stds(idx,:) = bsxfun( @minus, roi_stds_old(idx,:),mean(roi_stds_old(idx,:)));
        roi_stds(idx,:) = bsxfun( @rdivide, roi_stds(idx,:), std(roi_stds_old(idx,:)));
        
        roi_means(idx,:) = bsxfun( @minus, roi_means_old(idx,:),mean(roi_means_old(idx,:)));
        roi_means(idx,:) = bsxfun( @rdivide, roi_means(idx,:), std(roi_means_old(idx,:)));
        
        % session 2
        idx = ~session1_idx;
        roi_stds(idx,:) = bsxfun( @minus, roi_stds_old(idx,:),mean(roi_stds_old(idx,:)));
        roi_stds(idx,:) = bsxfun( @rdivide, roi_stds(idx,:), std(roi_stds_old(idx,:)));
        
        roi_means(idx,:) = bsxfun( @minus, roi_means_old(idx,:),mean(roi_means_old(idx,:)));
        roi_means(idx,:) = bsxfun( @rdivide, roi_means(idx,:), std(roi_means_old(idx,:)));
    end
else
    roi_stds = roi_stds_old;
    roi_means = roi_means_old;
end

%% Get new roi indices into old rois
[new_rois_old_idx] = get_new_roi_index_into_old_rois();

%% Compute ROI Ordering
filename = fullfile('..','FMRI','BrainVisualization','Network Classification','ROI_299_list.csv');
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
    h = figure(1); clf;
    h.Position =[92 0 1638 665];
    h1 = subplot(1,2,1);
    
    [unique_subjs,first_subj_occurances,~] = unique(combs(:,1), 'first'); % get first occurance of each unique element
    
    
    image(roi_stds(:,sort_idx), 'CDataMapping', colormapping); colorbar; hold on;
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
            plot([x,x], [1 size(roi_stds,1)], 'k', 'LineWidth', 1); hold on;
        end
    end
    
    [sorted_tasks, task_sort_idx] = sort(combs(:,2));
    
    h2 = subplot(1,2,2);
    
    [unique_tasks,first_task_occurances,~] = unique(sorted_tasks, 'first'); % get first occurance of each unique element
    image(roi_stds(task_sort_idx, sort_idx), 'CDataMapping', colormapping); colorbar; hold on;
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
    
    filename = fullfile('images', sprintf('functional_variance_subject_vs_task_heatmap_%s_normalize=%d', roi_sort_by, normalize));
    print(filename, '-depsc');
    
end
%% Plot For only 19 subjects in both sessions

if plot_19_subjs_subj
    
    bothsess = ismember(combs(:,1), subsinbothsessions);
    session1 = and(combs(:,3) == 1, bothsess);
    session2 = and(combs(:,3) == 2, bothsess);
    
    roi_stds_sess1 = roi_stds(session1, :);
    roi_stds_sess2 = roi_stds(session2, :);
    
    h = figure(2); clf;
    h.Position =[92 0 1638 665];
    
    % session 1
    h1 = subplot(1,2,1);
    
    [sess1_subjs, sort_idx] = sort(combs(session1, 1));
    [unique_subjs,first_subj_occurances,~] = unique(sess1_subjs, 'first'); % get first occurance of each unique element
    
    image(roi_stds_sess1(sort_idx, lobe_idx), 'CDataMapping', colormapping); colorbar;
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
    
    if remove_outlier
        [~,ii] = max(roi_stds_sess2(:));
        roi_stds_sess2(ii) = mean(roi_stds_sess2(:));
        image(roi_stds_sess2(sort_idx, :), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    else
        image(roi_stds_sess2(sort_idx,lobe_idx), 'CDataMapping', colormapping); colorbar;
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
    filename = fullfile('images', sprintf('FV_19_subjs_heatmap_subjorder_%s_normalize=%d', roi_sort_by, normalize));
    print(filename, '-depsc');
    
end


%% Plot For only 19 subjects in both sessions

if plot_19_subjs_task
    
    
    bothsess = ismember(combs(:,1), subsinbothsessions);
    session1 = and(combs(:,3) == 1, bothsess);
    session2 = and(combs(:,3) == 2, bothsess);
    
    roi_stds_sess1 = roi_stds(session1, :);
    roi_stds_sess2 = roi_stds(session2, :);
    
    taskid = combs(:,2);
    [sess1_tasks, sort_idx] = sort(taskid(session1));
    
    if strcmp(roi_sort_by, 'Lobe')
        [unique_lobes,first_lobe_occurances,~] = unique(regionlobes(lobe_idx), 'first'); % get first occurance of each unique element
    elseif strcmp(roi_sort_by, 'Hemisphere')
        [unique_lobes,first_lobe_occurances,~] = unique(hemisphere(hem_idx), 'first'); % get first occurance of each unique element
    end
    
    h = figure(3); clf;
    h.Position =[92 0 1638 665];
    
    % session 1
    h1 = subplot(1,2,1);
    
    [unique_tasks,first_task_occurances,~] = unique(sess1_tasks, 'first'); % get first occurance of each unique element
    image(roi_stds_sess1(sort_idx,lobe_idx), 'CDataMapping', colormapping); colorbar;
    
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
    
    if remove_outlier
        [~,ii] = max(roi_stds_sess2(:));
        roi_stds_sess2(ii) = mean(roi_stds_sess2(:));
        image(roi_stds_sess2(sort_idx,lobe_idx), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    else
        image(roi_stds_sess2(sort_idx,lobe_idx), 'CDataMapping', colormapping); colorbar;
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
    
    filename = fullfile('images', sprintf('FV_19_subjs_heatmap_taskorder_%s_normalize=%d', roi_sort_by, normalize));
    print(filename, '-depsc');
    
end


%% Plot FV ordered by highest variance for each subject

if plot_sanity_check
    
    % sort each subject by highest var
    roi_stds_sort = zeros(size(roi_stds));
    for s = 1:size(roi_stds,1)
        roi_stds_sort(s,:) = sort(roi_stds(s,:));
    end
    
    % make some random data
    roi_stds_rand_sort = zeros(size(roi_stds));
    for s = 1:size(roi_stds,1)
        indices = datasample(1:size(roi_stds,1), NR); % sample a subject to take value from for each region
        
        rand_subj = zeros(1,NR);
        for r = 1:NR
            rand_subj(r) = roi_stds(indices(r), r); % choose ROI from subject we sampled
        end
        roi_stds_rand_sort(s,:) = rand_subj;
    end
    
    for s = 1:size(roi_stds,1)
        roi_stds_rand_sort(s,:) = sort(roi_stds_rand_sort(s,:));
    end
    
    h = figure(3); clf;
    h.Position =[92 0 600 600];
    
    % session 1
    subplot(1,1,1);
    
    h1 = plot(roi_stds_sort(combs(:,3) == 1, :)', 'Color', [0 0 1 0.5]); hold on;
    
    title(sprintf('Session 1 Sorted FVs'));
    xlabel('ROI (sorted by FV)');
    ylabel('FV');
    
    %h2 = subplot(1,2,2);
    
    h2 = plot(roi_stds_rand_sort(combs(:,3) == 1, :)', 'Color', [1 0 0 0.5]);
    
    title(sprintf('Session 1 Sorted FVs'));
    xlabel('ROI (sorted by FV)');
    ylabel('FV');
    
    h = legend([h1(1), h2(1)], {'sorted subjects', 'sorted random'});
    
    xlim([1,NR]);
    filename = fullfile('images', sprintf('functional_variance_sorted_vs_random_%s', roi_sort_by));
    print(filename, '-depsc');
end


%% Means
if plot_19_subjs_subj_mean
    
    colormapping = 'scaled';
    
    bothsess = ismember(combs(:,1), subsinbothsessions);
    session1 = and(combs(:,3) == 1, bothsess);
    session2 = and(combs(:,3) == 2, bothsess);
    
    roi_means_sess1 = roi_means(session1, :);
    roi_means_sess2 = roi_means(session2, :);
    
    h = figure(4); clf;
    h.Position =[92 0 1638 665];
    
    % session 1
    h1 = subplot(1,2,1);
    
    [sess1_subjs, sort_idx] = sort(combs(session1, 1));
    [unique_subjs,first_subj_occurances,~] = unique(sess1_subjs, 'first'); % get first occurance of each unique element
    
    image(roi_means_sess1(:,sort_idx), 'CDataMapping', colormapping); colorbar;
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
    
    if remove_outlier
        [~,ii] = max(roi_means_sess2(:));
        roi_means_sess2(ii) = mean(roi_means_sess2(:));
        image(roi_means_sess2(:,sort_idx), 'CDataMapping', colormapping); colorbar;
        if normalize == 0
            caxis([cmin, cmax]);
        end
        hold on;
    else
        image(roi_means_sess2(:,sort_idx), 'CDataMapping', colormapping); colorbar;
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
    filename = fullfile('images', sprintf('M_19_subjs_heatmap_subjorder_%s_normalize=%d', roi_sort_by, normalize));
    print(filename, '-dpng');
    
end


















