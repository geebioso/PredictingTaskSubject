function [] = plot_LR_all_accuacies(isHPC, classifiertype, truncate) 

[wd, rd] = set_directories(isHPC); % working directory and results directory

task_only = 1;
baseline = sprintf( 'whs%s', genvarname('0')); 

if ~ismember( classifiertype , [0,6])
   error( 'classifier type not supported for this function')  
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

set_descriptions0.(sprintf('whs%s', genvarname('-1'))) = {'90s Filter' 'No ICA' 'With Experimental Design' 'CSF+WM Regressed Out'};
set_descriptions0.(sprintf('whs%s', genvarname('0'))) = {'200s Filter','ICA','With Experimental Design','GM Regressed Out'};
set_descriptions0.(sprintf('whs%s', genvarname('1'))) = {'200s Filter','ICA','With Experimental Design','CSF+WM Regressed Out'};
set_descriptions0.(sprintf('whs%s', genvarname('2'))) = {'200s Filter','No ICA','With Experimental Design','GM Regressed Out'};
set_descriptions0.(sprintf('whs%s', genvarname('3'))) = {'200s Filter','ICA','Experimental Design Regressed Out','GM Regressed Out'};

%% Load Data 
filenm = fullfile( rd, 'Results', sprintf('AC_c%d_truncate%d', classifiertype, truncate));
load(filenm); 

baseline_idx = find( strcmp(fieldnames(AC), baseline)); 

if isempty(baseline_idx)
   error(['Accuracy not computed for baseline dataset. Either respecify', ...
       ' baseline dataset or run models for baseline dataset']); 
end

set_fields = fieldnames(AC);
rep_fields = fieldnames(AC.(set_fields{baseline_idx}));
NPT = size(AC.(baseline).BVSD,1); % number of prediction tasks 
NS  = length(set_fields); % number of datasets 
NR  = length(rep_fields); % number of feature representations 
pred_tasks = {'Task Classification\nWithin Session', 'Task Classification\nBetween Session', ...
    'Subject Classification\nWithin Session', 'Subject Classification\nBetween Session'};
pred_tasks = {'Within Session', 'Between Session', 'Within Session', 'Between Session'};


F = fieldnames(set_descriptions0);
D = struct2cell(set_descriptions0);
keep = ismember( F, set_fields);
set_descriptions = cell2struct(D(keep), F(keep));

% set_fields = sort(set_fields); 
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
COR = cell(NPT,length(set_fields), length(rep_fields));
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

A(A==0) = NaN;

task_ymin = min(min(min(A(1:2,:,:))));
task_ymax = max(max(max(A(1:2,:,:))));
subj_ymin = min(min(min(A(3:4,:,:))));
subj_ymax = max(max(max(A(3:4,:,:))));

ncol = 2;
nrow = 2;
if task_only
    NPT = 2;
    A = A(1:2,:,:);
    C = C(1:2,:,:,:); 
    COR = COR(1:2,:,:);
   
    nrow = 1;
end

%% Compute Overall Credible intervals between datasets (over all feature types) 
% credible interval size 
p = 0.95; 

% hypterparameters for mean at chance
alpha = 1.25;
beta = 3;
                
JOINTCOR = cell( NPT, NS);
JOINTA = zeros(  NPT, NS);
JOINTCRED = zeros(  NPT, NS, 2);
for t = 1:NPT
   for s = 1:NS
      JOINTCOR{t,s} = vertcat( COR{t,:,s}); 
      JOINTA(t,s) = mean(JOINTCOR{t,s}); 
      
      [ cred, pbayesmean, alpha_post, beta_post] = beta_binomial_cred_interval(alpha, beta, JOINTCOR{t,s}, p);
      
      JOINTCRED(t,s,:) = cred; 
   end
end



%% Plot By Prediction Task
markers = '+o*xsd^';
map = brewermap(length(set_fields), 'Set1');
clear lines

h0 = figure(); clf;
if task_only
    h0.Position = [152 377 998 328];
else
    h0.Position = [152 41 998 664];
end
%index = reshape(1:6, 3, 2).';
for i = 1:NPT
    %h = subplot(2,3,index(i));
    h = subplot(nrow,ncol,i);
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
    h.XTickLabelRotation = 45;
    if i < 3
        ylim([task_ymin, task_ymax]);
    else
        ylim([subj_ymin, subj_ymax]);
    end
    title(pred_tasks{i});
    %if i == 3
    legend(lines, leg_labels_together, 'Location', 'SouthEast'); % , set_fields);
    ylabel('accuracy') 
    %end
end

filenm = fullfile('..', 'images', 'accuracies_all_datasets');
print(filenm, '-depsc');

%% Pairwise Subplot Bar

delta = 0.2; 
allpairs = nchoosek(set_fields,2);
tokeep = [];
for p = 1:size(allpairs,1)
    if strcmp( allpairs{p,1}, baseline)
        tokeep = [tokeep, p];
    end
end
allpairs = allpairs(tokeep,:);

NP = size(allpairs,1);

for p = 1:NP
    h0 = figure(); clf;
    if task_only
         h0.Position = [371 356 670 296];
    else
        h0.Position = [1 193 670 512];
    end
    pair_now = allpairs(p,:);
    
    % find indices into accuracy for the datasets we are comparing
    set1_idx = ismember( set_fields, pair_now{1});
    set2_idx = ismember( set_fields, pair_now{2});
    
    % for each prediction task 
    for t = 1:NPT
        subplot(nrow,ncol,t);
        
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
        
        % plot bar 
        h = bar( pairA );
        
        % plot significance 
        SIGidx = find(SIG); 
        if ~isempty(SIGidx)
           left = SIGidx - delta; 
           right = SIGidx + delta; 
           
           sigstar_input = cell(length(SIGidx),1); 
           for i = 1:length(SIGidx,1)
                sigstar_input{i}(1) = left(i);
                sigstar_input{i}(2) = right(i);
           end
           
           sigstar(sigstar_input); 
        end
        
        set(gca, 'XTickLabel', strrep( rep_fields, '_', ' '));
        set(gca, 'XTickLabelRotation', 0);
        
        
        title(pred_tasks{t})
        
    end
    leg_idx = ~strcmp( leg_labels{set1_idx}, leg_labels{set2_idx});
    
    suptitle(sprintf('%s (Blue) vs %s (Red)', leg_labels{set1_idx}{leg_idx}, ...
        leg_labels{set2_idx}{leg_idx})  );
    
    filenm = fullfile( '..', 'images', sprintf('accuracies_bar_%s_vs_%s', ...
        pair_now{1}, pair_now{2}));
    print(filenm, '-depsc', '-r300');
    
end

%% Pairwise scatter
if ~task_only
    allpairs = nchoosek(set_fields,2);
    tokeep = [];
    for p = 1:size(allpairs,1)
        if strcmp( allpairs{p,1}, set_fields{1})
            tokeep = [tokeep, p];
        end
    end
    allpairs = allpairs(tokeep,:);
    NP = size(allpairs,1);
    
    labels = repmat(rep_fields, 4, 1);
    
    markers = '^sxo';
    
    for p = 1:NP
        h0 = figure(); clf;
        h0.Position = [1 193 670 512];
        pair_now = allpairs(p,:);
        
        set1_idx = ismember( set_fields, pair_now{1});
        set2_idx = ismember( set_fields, pair_now{2});
        
        h = plot(A(4,set1_idx,1), 1 - 0.1, 'r-');  hold on;
        line_handles(1) = h(1);
        h = plot(A(4,set2_idx,1), 1 - 0.1, 'b-');  hold on;
        line_handles(2) = h(1);
        
        for t = 1:NPT
            markernow = markers(t);
            
            for r = 1:NR
                
                if r == 1
                    h = plot( A(t,set1_idx,r), r - 0.1, [ 'k' markernow], 'MarkerSize', 9);  hold on;
                    marker_handles(t) = h(1);
                end
                plot( A(t,set1_idx,r), r - 0.1, [ 'r' markernow], 'MarkerSize', 9);  hold on;
                plot( A(t,set2_idx,r), r + 0.1, [ 'b' markernow], 'MarkerSize', 9);  hold on;
            end
        end
        
        set(gca, 'YTick', 1:NR);
        set(gca, 'YTickLabel', strrep( rep_fields, '_', ' '));
        
        leg_idx = ~strcmp( leg_labels{set1_idx}, leg_labels{set2_idx});
        
        title(sprintf('%s vs %s', leg_labels{set1_idx}{leg_idx}, ...
            leg_labels{set2_idx}{leg_idx})  );
        
        xlabel('accuracy');
        
        l = legend(marker_handles, pred_tasks, 'Location', 'NorthWest');
        a=axes('position',get(gca,'position'),'visible','off');
        leg_idx = ~strcmp( leg_labels{set1_idx}, leg_labels{set2_idx});
        l2 = legend(a, line_handles, { leg_labels{set1_idx}{leg_idx} leg_labels{set2_idx}{leg_idx}} );
        % l2.Position = [0.1482    0.7833-0.1179/1.6    0.1955    0.1179/2];
        l2.Position = [0.1448 0.7930-0.1123/1.6 0.1843 0.1123/2];
        
        filenm = fullfile( '..', 'images', sprintf('accuracies_scatter_%s_vs_%s', ...
            pair_now{1}, pair_now{2}));
        print(filenm, '-depsc', '-r300');
        
    end
    
end


