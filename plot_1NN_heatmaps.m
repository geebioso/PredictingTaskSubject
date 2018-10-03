function [] = plot_1NN_heatmaps( isHPC, whs, truncate )
%% Subject Classification Plots

[wd, rd] = set_directories(isHPC); % working directory and results directory

%% Options
classifiertype = 8;
FC_and_FV = 1;
subj_plot = 1;

%% Load the original ROI timecourse data

if ~exist( 'tasks' , 'var' )
    load( fullfile(rd, 'ICAresults','tasklist') );
end


filenm = fullfile( rd, 'features', sprintf('whs%d_truncate%d.mat', whs, truncate));
fprintf('Loading %s\n', filenm);
load( filenm, 'combs' );

% constants
NS = length( unique( combs(:, 1) ) ); % number of subjects
NT = length( unique( combs(:, 2)) ); % number of tasks

% Find subjects that are in session 1 and session 2
nspersubj = accumarray( combs(:,1) , combs(:,3) , [ NS 1 ] , @mean );
subsinbothsessions = find( nspersubj == 1.5 );
NSboth = length(subsinbothsessions);

%% Load Results

filenm = fullfile( rd, 'Results', sprintf('AC_c%d_truncate%d', classifiertype, truncate));
load(filenm); % loads AC, CR, CORRECT


%% Plot both on same Figure

whreps = {'BVSD', 'FCP'};
dopredtasks = {'between', 'within'};
count = 1;


fontsize = 18;
labelbounds = 1:9;
lim = [1,9];
nanmat = zeros(NT);
nanmat(1:size(nanmat,1)+1:end) = nan;
h0 = figure(); clf;
h0.Position = [139 5 905 700];

ha = tight_subplot(2,2,[.05 .05],[.27 .05],[.2 .01]);
cmin = 0; cmax = 1;

for whrep = whreps
    for pt = dopredtasks
        
        corrects = CORRECT.(sprintf('whs%s', genvarname(num2str(whs)))).(whrep{:}).(pt{:});
        
        if strcmp( pt{:}, 'between')
            task_accuracies = mean( corrects, 3);
            accuracy = mean( corrects(:) );
            
        elseif strcmp( pt{:}, 'within')
            accuracy = mean( corrects, 3 );
            idx = eye( size( accuracy, 1));
            accuracy = accuracy( ~idx );
            accuracy = mean( accuracy(:) );
            task_accuracies = mean( corrects, 3);
            task_accuracies(eye(size(task_accuracies,1)) == 1) = nan; 
        end
        
        fprintf('Mean Accuracy %s %s: %2.4f\n', whrep{:}, pt{:}, accuracy);
        
        axes(ha(count));
        h = ha(count);
        %h = subplot(2,2,1);
        imagesc(task_accuracies); colorbar;
        caxis([ cmin cmax ]);
        set(gca, 'YTick', labelbounds);
        set(gca, 'YTickLabels', tasks);
        h.FontSize = fontsize;
        
        if count == 1
            ylabel('Training task');
            h.TickLength = [0 0];
            h.XTick = [];
        elseif count == 2
            h.XTick = [];
            h.YTick = [];
        elseif count == 3
            ylabel('Training task');
            xlabel('Test task');
            h.XTickLabels = tasks;
            h.XTick =  labelbounds;
            h.XTickLabelRotation = 90;
            h.TickLength = [0 0];
        elseif count == 4
            xlabel('Test task');
            h.XTickLabels = tasks;
            h.XTick =  labelbounds;
            h.XTickLabelRotation = 90;
            h.TickLength = [0 0];
            h.YTick = [];
        end
        
        title(sprintf('%s %s-session accuracy', whrep{:}, pt{:}));
        count = count + 1;
        
    end
end


filenm = fullfile(rd, 'images', sprintf( 'finn_et_al_replication_newsplit_pt2-3_whs%d_whr1&4_truncate%d', whs, truncate));
print(filenm, '-depsc', '-r300');




















