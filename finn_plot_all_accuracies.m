
% NOTE: THE ORDER OF DATASETS IN AC IS IMPORTANT!!!!

plot_type = 'finn'; % 'finn' or 'LR'

if strcmp(plot_type, 'LR')
   pred_mode = 'task'; 
else
   pred_mode = 'subj';  
end


% Default font size
set(0,'defaultAxesFontSize',12);
addpath('../../MATLAButils/')

% Get AC by running print_result_table_regularization_optimized.m
% types of feature

feat_types = struct();
feat_types.BVSD = 'BV';
feat_types.BVM = 'M';
feat_types.BVVS = 'BV';
feat_types.BVV = 'BV';
feat_types.MSSD = 'BV';
feat_types.SQRT_MSSD = 'BV';
feat_types.FCP = 'FC';
feat_types.FCC = 'FC';
feat_types.FCCV = 'FC';
feat_types.FCCS = 'FC';

set_descriptions0 = struct();

% PUT BASELINE FIRST!!!!!
set_descriptions0.whset26 = {'200s Filter','ICA','With Experimental Design','GM Regressed Out'};


set_descriptions0.whset8 = {'90s Filter' 'No ICA' 'With Experimental Design' 'CSF+WM Regressed Out'};
set_descriptions.whset20 = {'200s Filter','ICA','Experimental Design Regressed Out','CSF+WM Regressed Out'};
set_descriptions0.whset21 = {'200s Filter','No ICA','With Experimental Design','CSF+WM Regressed Out'};
set_descriptions0.whset22 = {'200s Filter','ICA','With Experimental Design','CSF+WM Regressed Out'};
set_descriptions0.whset23 = {'200s Filter','No ICA','With Experimental Design','GM Regressed Out'};
set_descriptions0.whset24 = {'200s Filter','No ICA','Experimental Design Regressed Out','CSF+WM Regressed Out'};
set_descriptions0.whset27 = {'200s Filter','ICA','Experimental Design Regressed Out','GM Regressed Out'};

if strcmp(plot_type, 'finn') % get from finn_et_al_analysis.m
    classifiertype = 8;
elseif strcmp(plot_type, 'LR') 
    classifiertype = 0; 
end

filenm = fullfile( rd, 'Results', sprintf('AC_c%d', classifiertype)); 
load(filenm); 

set_fields = fieldnames(AC);
rep_fields = fieldnames(AC.(set_fields{1}));
NPT = size(AC.whset26.BVSD,1);
NS  = length(set_fields);
NR  = length(rep_fields);
if strcmp(plot_type, 'LR') 
    pred_tasks = {'Task Within', 'Task Between', 'Subject Within', 'Subject Between'};
elseif strcmp(plot_type, 'finn')
    pred_tasks = {'Within-Session', 'Between-Session', 'Task-Same Task', 'Task-Different Task', 'Rest-Rest'};
end

F = fieldnames(set_descriptions0);
D = struct2cell(set_descriptions0);
keep = ismember( F, set_fields);
set_descriptions = cell2struct(D(keep), F(keep));

if ~isequal( fieldnames(set_descriptions), set_fields)
    error( 'fieldnames of accuracies and set_descriptions don''t match');
end

leg_labels = struct2cell(set_descriptions);
leg_labels_together = cell(1);
for i = 1:size(leg_labels,1)
    temp = leg_labels{i,:};
    temp = cellfun( @(x) [x ','], temp, 'UniformOutput', 0); temp{end} = temp{end}(1:end-1);
    leg_labels_together{i} = [temp{:}];
end

% populate matrix with data
A = zeros( NPT, length(set_fields), length(rep_fields));
C = zeros( NPT, length(set_fields), length(rep_fields), 2);
COR = cell( NPT,length(set_fields), length(rep_fields));
for i = 1:NS
    set_field = set_fields{i};
    for j = 1:NR
        rep_field = rep_fields{j};
        A(:,i,j) = AC.(set_field).(rep_field);
        C(:,i,j,:) = CR.(set_field).(rep_field); 
        
        for t = 1:NPT
         COR{t,i,j} = CORRECT.(set_field).(rep_field){t}; 
        end
    end
end

%A(A==0) = NaN;

if strcmp(plot_type, 'finn')
    ncol = 3;
    nrow = 2;
    subj_ymin = min(min(min(A(:,:,:))));
    subj_ymax = max(max(max(A(:,:,:))));
elseif strcmp(plot_type, 'LR')
    ncol = 2;
    nrow = 2;
    task_ymin = min(min(min(A(1:2,:,:))));
    task_ymax = max(max(max(A(1:2,:,:))));
    subj_ymin = min(min(min(A(3:4,:,:))));
    subj_ymax = max(max(max(A(3:4,:,:))));
end
if strcmp(pred_mode, 'task')
    NPT = 2;
    A = A(1:2,:,:);
    
    nrow = 1;
end

%% Plot By Prediction Task
markers = '+o*xsd^';
map = brewermap(length(set_fields), 'Set1');
clear lines

h0 = figure(1); clf;
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
if classifiertype ~= 8
    h0.Position = [152 377 998 328];
else
    h0.Position = [152 41 998 664];
end
%index = reshape(1:6, 3, 2).';
for i = 1:NPT
    %h = subplot(2,3,index(i));
    if and( strcmp(plot_type, 'finn'), i>2)
        h = subplot(nrow,ncol,i+1);
    else
        h = subplot(nrow,ncol,i);
    end
    for j = 1:NS
        for k = 1:NR
            if feat_types.(rep_fields{k}) == 'FC'
                markernow = '+';
            elseif feat_types.(rep_fields{k}) == 'BV'
                markernow = 'o';
            else
                markernow = '^';
            end
            plot( k, squeeze( A(i,j,k)), markernow, 'Color', map(j,:)); hold on;
        end
        
        p = plot( 1:NR, squeeze( A(i,j,:)), '-', 'Color', map(j,:)); hold on;
        lines(j) = p(1);
        
    end
    h.XTick = 1:NR;
    h.XTickLabel = strrep( rep_fields, '_', ' ');
    h.XTickLabelRotation = 0;
    if strcmp(plot_type, 'LR')
        if i < 3
            ylim([task_ymin, task_ymax]);
        else
            ylim([subj_ymin, subj_ymax]);
        end
    else
        ylim([subj_ymin, subj_ymax]);
    end
    title(pred_tasks{i});
    if or( and( i == 1, strcmp(plot_type, 'finn')), strcmp(plot_type, 'LR'))
        legend(lines, leg_labels_together, 'Location', 'Best'); % , set_fields);
    end
    ylabel('accuracy')
    %end
end

if strcmp(plot_type, 'finn')
    filenm = fullfile('..', 'images', 'finn_accuracies_all_datasets');
elseif strcmp(plot_type, 'LR') 
    filenm = fullfile('..', 'images', 'accuracies_all_datasets');
end
print(filenm, '-dtiff');

%% Pairwise Subplot Bar

allpairs = nchoosek(set_fields,2);
tokeep = [];
for p = 1:size(allpairs,1)
    if strcmp( allpairs{p,1}, set_fields{1})
        tokeep = [tokeep, p];
    end
end
allpairs = allpairs(tokeep,:);

NP = size(allpairs,1);

for p = 1:NP
    h0 = figure(p+NP+1); clf;
    set(h0, 'DefaultAxesFontName', 'Arial');
    set(h0, 'DefaultTextFontName', 'Arial');
    if classifiertype ~=8 
        h0.Position = [371 356 670 296];
    else
        h0.Position = [90 160 836 545];
    end
    pair_now = allpairs(p,:);
    
    set1_idx = ismember( set_fields, pair_now{1});
    set2_idx = ismember( set_fields, pair_now{2});
    
    for t = 1:NPT
        if and( strcmp(plot_type, 'finn'), t>2)
            h = subplot(nrow,ncol,t+1);
        else
            h = subplot(nrow,ncol,t);
        end
        
        % paired accuracies 
        pairA = [squeeze( A(t,set1_idx,:) ) squeeze( A(t,set2_idx,:) )]; 
        
        % paired credible intervals 
        C1 =  squeeze( C(t,set1_idx,:,:) );
        C2 =  squeeze( C(t,set2_idx,:,:));
        
        % compute significant differences 
        SIG = zeros(NR,1); 
        for n = 1:NR
           if pairA(n,1) > pairA(n,2)
              if C1(n,1) > C2(n,2) 
                 SIG(n) = 1;  
              end
           else
               if C1(n,2) < C2(n,1) 
                 SIG(n) = 1;  
              end
           end
        end
        
        h = bar( [squeeze( A(t,set1_idx,:) ) squeeze( A(t,set2_idx,:) )]);
        set(gca, 'XTickLabel', strrep( rep_fields, '_', ' '));
        set(gca, 'XTickLabelRotation', 0);
        
        title(pred_tasks{t})
        
    end
    leg_idx = ~strcmp( leg_labels{set1_idx}, leg_labels{set2_idx});
    
    suptitle(sprintf('%s (Blue) vs %s (Red)', leg_labels{set1_idx}{leg_idx}, ...
        leg_labels{set2_idx}{leg_idx})  );
    
    if strcmp( plot_type, 'finn')
        filenm = fullfile( '..', 'images', sprintf('finn_accuracies_bar_%s_vs_%s', ...
            pair_now{1}, pair_now{2}));
    elseif strcmp(plot_type, 'LR') 
        filenm = fullfile( '..', 'images', sprintf('accuracies_bar_%s_vs_%s', ...
            pair_now{1}, pair_now{2}));
    end
    print(filenm, '-depsc');
    
end

% %% Pairwise scatter
% if ~task_only
%     allpairs = nchoosek(set_fields,2);
%     tokeep = [];
%     for p = 1:size(allpairs,1)
%         if strcmp( allpairs{p,1}, set_fields{1})
%             tokeep = [tokeep, p];
%         end
%     end
%     allpairs = allpairs(tokeep,:);
%     NP = size(allpairs,1);
%     
%     labels = repmat(rep_fields, 4, 1);
%     
%     markers = '^sxo';
%     
%     for p = 1:NP
%         h0 = figure(p +3*NP+1); clf;
%         h0.Position = [1 193 670 512];
%         pair_now = allpairs(p,:);
%         
%         set1_idx = ismember( set_fields, pair_now{1});
%         set2_idx = ismember( set_fields, pair_now{2});
%         
%         h = plot(A(4,set1_idx,1), 1 - 0.1, 'r-');  hold on;
%         line_handles(1) = h(1);
%         h = plot(A(4,set2_idx,1), 1 - 0.1, 'b-');  hold on;
%         line_handles(2) = h(1);
%         
%         for t = 1:NPT
%             markernow = markers(t);
%             
%             for r = 1:NR
%                 
%                 if r == 1
%                     h = plot( A(t,set1_idx,r), r - 0.1, [ 'k' markernow], 'MarkerSize', 9);  hold on;
%                     marker_handles(t) = h(1);
%                 end
%                 plot( A(t,set1_idx,r), r - 0.1, [ 'r' markernow], 'MarkerSize', 9);  hold on;
%                 plot( A(t,set2_idx,r), r + 0.1, [ 'b' markernow], 'MarkerSize', 9);  hold on;
%             end
%         end
%         
%         set(gca, 'YTick', 1:NR);
%         set(gca, 'YTickLabel', strrep( rep_fields, '_', ' '));
%         
%         leg_idx = ~strcmp( leg_labels{set1_idx}, leg_labels{set2_idx});
%         
%         title(sprintf('%s vs %s', leg_labels{set1_idx}{leg_idx}, ...
%             leg_labels{set2_idx}{leg_idx})  );
%         
%         xlabel('accuracy');
%         
%         l = legend(marker_handles, pred_tasks, 'Location', 'NorthWest');
%         a=axes('position',get(gca,'position'),'visible','off');
%         leg_idx = ~strcmp( leg_labels{set1_idx}, leg_labels{set2_idx});
%         l2 = legend(a, line_handles, { leg_labels{set1_idx}{leg_idx} leg_labels{set2_idx}{leg_idx}} );
%         % l2.Position = [0.1482    0.7833-0.1179/1.6    0.1955    0.1179/2];
%         l2.Position = [0.1448 0.7930-0.1123/1.6 0.1843 0.1123/2];
%         
%         filenm = fullfile( '..', 'images', sprintf('accuracies_scatter_%s_vs_%s', ...
%             pair_now{1}, pair_now{2}));
%         print(filenm, '-dpng');
%         
%     end
%     
% end


