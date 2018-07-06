%% Subject Classification Plots


%% Options

whs = 8;
classifiertype = 8;
whrep = 1;
pt = 3;
FC_and_FV = 1; 
subj_plot = 1; 
task_plot = 1; 

newtasklabels = {'Emotional Pictures', 'Emotional Faces', ...
    'Encoding', 'Rest', 'Retrieval', 'Go/No-go', 'Monetary Incentive', 'Working Memory', 'Theory of Mind'}; 

%% Load the original ROI timecourse data

if (whs==8)
    roifile    = ['ICAresults' filesep 'tc_csf_wm_motion_out_globalMask5000_203subj_conv'];
    icfile     = '';
end

if ~exist( 'D' , 'var' )
    load( roifile );
end

if ~exist( 'tasks' , 'var' )
    load( fullfile('ICAresults','tasklist') );
end

tasks = newtasklabels; 

NS = max( D.SubjectIndex2 ); % number of subjects
NT = length( tasks ); % number of tasks
NR = size( D.ROI , 2 ); % number of regions

%% Get combinations of Task x Subject x Session
[ combs , ~ , cindex ] = unique( [ D.SubjectIndex2 D.TaskIndex D.SessionIndex ], 'rows' );
NC = size( combs , 1 );

% Find subjects that are in session 1 and session 2
nspersubj = accumarray( combs(:,1) , combs(:,3) , [ NS 1 ] , @mean );
subsinbothsessions = find( nspersubj == 1.5 );


if subj_plot

filename = fullfile('Results', sprintf('PREDIDs_pt%d_whs%d_whrep%d_classifiertype%d.mat', ...
    pt, whs, whrep, classifiertype));
load(filename);

W_CORRECT = CORRECT;

pt = 4;
filename = fullfile('Results', sprintf('PREDIDs_pt%d_whs%d_whrep%d_classifiertype%d.mat', ...
    pt, whs, whrep, classifiertype));
load(filename);
B_CORRECT = CORRECT;


%% Calculate accuracies
b_accuracy = mean(B_CORRECT, 3);
w_accuracy = mean(W_CORRECT, 3);

NT = 9;
NSboth = 19;

idx = eye(NT,NT);
%Y = (1-idx).*X;
w_accuracy_nondiagonal = w_accuracy(~idx);

b_accuracy_overall = mean(b_accuracy(:));
w_accuracy_overall = mean(w_accuracy_nondiagonal);

fprintf('Mean Accuracy is\n\tBetween: %2.4f\n\tWithin: %2.4f\n', mean(b_accuracy(:)), w_accuracy_overall);

%% Plots
fontsize = 18;
labelbounds = 1:9;
lim = [1,9];
nanmat = zeros(NT);
nanmat(1:size(nanmat,1)+1:end) = nan;
h0 = figure(5); clf;
h0.Position = [94 497 1.8*1183 2*445];

cmin = 0; cmax = 1;
h = subplot(1,2,1);
imagesc(b_accuracy); colorbar;
caxis([ cmin cmax ]);
set(gca, 'YTick', labelbounds);
set(gca, 'YTickLabels', tasks);
ylabel('training task');
xlabel('test task');
h.XTickLabels = tasks;
h.XTick =  labelbounds;
h.XTickLabelRotation = 90;
h.TickLength = [0 0];
h.FontSize = fontsize;
%h.XLim = lim;
%title(sprintf('Between Session Subject\n Identity Classification Accuracy'));
% title('A'); 
title('FC between-session accuracy');

h = subplot(1,2,2);
imagesc(w_accuracy); colorbar;
%patchimg((taskaccuracies{1,2} + taskaccuracies{2,2})/2 + nanmat, 'none'); colorbar;
caxis([ cmin cmax ]);
set(gca, 'YTick',  labelbounds);
set(gca, 'YTickLabels', tasks);
ylabel('training task');
xlabel('test task');
h.XTickLabels = tasks;
h.XTick =  labelbounds;
h.XTickLabelRotation = 90;
h.TickLength = [0 0];
h.FontSize = fontsize;
%h.XLim = lim;
%title(sprintf('Within Session Subject\n Identity Classification Accuracy'));
%title('B');
title('FC within-session accuracy');

filenm = fullfile('images', sprintf( 'finn_et_al_replication_newsplit_pt2-3_whs%d_whr%d', whs , whrep ));
print(filenm, '-depsc', '-r300');


end

%% Task Classification Plots

if task_plot 
whs = 8;
classifiertype = 8;
whrep = 4;

pt = 1;
filename = fullfile('Results', sprintf('PREDIDs_pt%d_whs%d_whrep%d_classifiertype%d.mat', ...
    pt, whs, whrep, classifiertype));
load(filename);

W_CORRECT = CORRECT;

pt = 2;
filename = fullfile('Results', sprintf('PREDIDs_pt%d_whs%d_whrep%d_classifiertype%d.mat', ...
    pt, whs, whrep, classifiertype));
load(filename);
B_CORRECT = CORRECT;


%% Calculate accuracies
b_accuracy = mean(B_CORRECT, 3);
w_accuracy = mean(W_CORRECT, 3);

idx = eye(NT,NT);
%Y = (1-idx).*X;
w_accuracy_nondiagonal = w_accuracy(~idx);

b_accuracy_overall = mean(b_accuracy(:));
w_accuracy_overall = mean(w_accuracy_nondiagonal);

fprintf('Mean Accuracy is\n\tBetween: %2.4f\n\tWithin: %2.4f\n', mean(b_accuracy(:)), w_accuracy_overall);

%% Plots

fontsize = 18;
labelbounds = 1:NSboth;
lim = [1,9];
nanmat = zeros(NT);
nanmat(1:size(nanmat,1)+1:end) = nan;
h0 = figure(6); clf;
h0.Position = [94 497 1.8*1183 2*445];

cmin = 0; cmax = 1;
h = subplot(1,2,1);
imagesc(b_accuracy); colorbar;
caxis([ cmin cmax ]);
set(gca, 'YTick', labelbounds);
set(gca, 'YTickLabels', subsinbothsessions);
ylabel('subject id known');
xlabel('subject id unknown');
h.XTickLabels = subsinbothsessions;
h.XTick =  labelbounds;
h.XTickLabelRotation = 90;
h.TickLength = [0 0];
h.FontSize = fontsize;
%h.XLim = lim;
%title(sprintf('Between Session\n Task Classification Accuracy'));
title('A');

h = subplot(1,2,2);
imagesc(w_accuracy); colorbar;
%patchimg((taskaccuracies{1,2} + taskaccuracies{2,2})/2 + nanmat, 'none'); colorbar;
caxis([ cmin cmax ]);
set(gca, 'YTick',  labelbounds);
set(gca, 'YTickLabels', subsinbothsessions);
ylabel('subject id known');
xlabel('subject id unknown');
h.XTickLabels = subsinbothsessions;
h.XTick =  labelbounds;
h.XTickLabelRotation = 90;
h.TickLength = [0 0];
h.FontSize = fontsize;
%h.XLim = lim;
%title(sprintf('Within Session\n Task Classification Accuracy'));
title('B');

filenm = fullfile('images', sprintf( 'finn_et_al_replication_newsplit_pt1-2_whs%d_whr%d', whs , whrep ));
%print(filenm, '-dpng', '-r0');
end



%% Plot both on same Figure
whrep = 4;
pt = 3;
filename = fullfile('Results', sprintf('PREDIDs_pt%d_whs%d_whrep%d_classifiertype%d.mat', ...
    pt, whs, whrep, classifiertype));
load(filename);

W_CORRECT_FC = CORRECT;

pt = 4;
filename = fullfile('Results', sprintf('PREDIDs_pt%d_whs%d_whrep%d_classifiertype%d.mat', ...
    pt, whs, whrep, classifiertype));
load(filename);
B_CORRECT_FC = CORRECT;

whrep = 1;
pt = 3;
filename = fullfile('Results', sprintf('PREDIDs_pt%d_whs%d_whrep%d_classifiertype%d.mat', ...
    pt, whs, whrep, classifiertype));
load(filename);

W_CORRECT_FV = CORRECT;

pt = 4;
filename = fullfile('Results', sprintf('PREDIDs_pt%d_whs%d_whrep%d_classifiertype%d.mat', ...
    pt, whs, whrep, classifiertype));
load(filename);
B_CORRECT_FV = CORRECT;


%% Calculate accuracies
b_accuracy_fc = mean(B_CORRECT_FC, 3);
w_accuracy_fc = mean(W_CORRECT_FC, 3);

b_accuracy_fv = mean(B_CORRECT_FV, 3);
w_accuracy_fv = mean(W_CORRECT_FV, 3);

NT = 9;
NSboth = 19;

idx = eye(NT,NT);
%Y = (1-idx).*X;
w_accuracy_nondiagonal_fc = w_accuracy_fc(~idx);

b_accuracy_overall_fc = mean(b_accuracy_fc(:));
w_accuracy_overall_fc = mean(w_accuracy_nondiagonal_fc);

w_accuracy_nondiagonal_fv = w_accuracy_fv(~idx);

b_accuracy_overall_fv = mean(b_accuracy_fv(:));
w_accuracy_overall_fv = mean(w_accuracy_nondiagonal_fv);

fprintf('Mean FC Accuracy is\n\tBetween: %2.4f\n\tWithin: %2.4f\n', mean(b_accuracy_fc(:)), w_accuracy_overall_fc);
fprintf('Mean FV Accuracy is\n\tBetween: %2.4f\n\tWithin: %2.4f\n', mean(b_accuracy_fv(:)), w_accuracy_overall_fv);

%% Plots
switch FC_and_FV
    case 1
        
      

        fontsize = 18;
        labelbounds = 1:9;
        lim = [1,9];
        nanmat = zeros(NT);
        nanmat(1:size(nanmat,1)+1:end) = nan;
        h0 = figure(7); clf;
        h0.Position = [139 5 905 700];
        
        ha = tight_subplot(2,2,[.05 .05],[.27 .05],[.2 .01]);
        
        cmin = 0; cmax = 1;
        axes(ha(1)); 
        h = ha(1); 
        %h = subplot(2,2,1);
        imagesc(b_accuracy_fc); colorbar;
        caxis([ cmin cmax ]);
        set(gca, 'YTick', labelbounds);
        set(gca, 'YTickLabels', tasks);
        ylabel('training task');
        h.TickLength = [0 0];
        h.XTick = []; 
        h.FontSize = fontsize;
        %h.XLim = lim;
        %title(sprintf('Between Session Subject\n Identity Classification Accuracy'));
        %title('A');
        title('FC between-session accuracy');
        
        axes(ha(2)); 
        h = ha(2); 
        imagesc(w_accuracy_fc); colorbar;
        %patchimg((taskaccuracies{1,2} + taskaccuracies{2,2})/2 + nanmat, 'none'); colorbar;
        caxis([ cmin cmax ]);
        set(gca, 'YTick',  labelbounds);
        set(gca, 'YTickLabels', tasks);
        h.XTick = [];
        h.YTick = []; 
        h.FontSize = fontsize;
        %h.XLim = lim;
        %title(sprintf('Within Session Subject\n Identity Classification Accuracy'));
        %title('B');
        title('FC within-session accuracy');
        
        
        axes(ha(3)); 
        h = ha(3); 
        imagesc(b_accuracy_fv); colorbar;
        caxis([ cmin cmax ]);
        set(gca, 'YTick', labelbounds);
        set(gca, 'YTickLabels', tasks);
        ylabel('training task');
        xlabel('test task');
        h.XTickLabels = tasks;
        h.XTick =  labelbounds;
        h.XTickLabelRotation = 90;
        h.TickLength = [0 0];
        h.FontSize = fontsize;
        %h.XLim = lim;
        %title(sprintf('Between Session Subject\n Identity Classification Accuracy'));
        %title('C');
        title('BV between-session accuracy');
        
        axes(ha(4)); 
        h = ha(4); 
        imagesc(w_accuracy_fv); colorbar;
        %patchimg((taskaccuracies{1,2} + taskaccuracies{2,2})/2 + nanmat, 'none'); colorbar;
        caxis([ cmin cmax ]);
        set(gca, 'YTick',  labelbounds);
        set(gca, 'YTickLabels', tasks);
        %ylabel('training task');
        xlabel('test task');
        h.XTickLabels = tasks;
        h.XTick =  labelbounds;
        h.XTickLabelRotation = 90;
        h.TickLength = [0 0];
        h.FontSize = fontsize;
        h.YTick = []; 
        %h.XLim = lim;
        %title(sprintf('Within Session Subject\n Identity Classification Accuracy'));
        %title('D');
        title('BV within-session accuracy');
        
        
        filenm = fullfile('images', sprintf( 'finn_et_al_replication_newsplit_pt2-3_whs%d_whr1&4', whs));
        print(filenm, '-depsc', '-r300');
        
        
end



















