function [ AC , CR, CORRECT, best_costs, tab ] = print_result_table_regularization_optimized( classifiertype, whsets, isHPC, truncate)

[wd, rd] = set_directories(isHPC);

if classifiertype ~= 8
    pred_task = 'task'; % both
else
    error( 'DO NOT USE classifiertype=8. DEPRECATED. USE finn_et_al_analysis.m INSTEAD');
end

whreps = [1,4, 7:9, 11]; % 10]; % [1,2,4,9,10]
whreps_labels = { 'BVSD', 'BVM', 'BVVS', 'BVV', 'MSSD', 'SQRT_MSSD', 'FCP', 'FCC', 'FCCV', 'FCCS', 'BVMAD'};
whreps_type = {'BV', 'BV', 'BV', 'BV', 'BV', 'BV', 'FC', 'FC', 'FC', 'FC', 'BV' };

% whrep==1: ROI StdDev
% whrep==2: ROI means
% whrep==3: ROI FV scaled
% whrep==4: ROI FV
% whrep==5: MSSD
% whrep==6: SQRT MSSD
% whrep==7: ROI FC correlation
% whrep==8: ROI FC covariance
% whrep==9: ROI FC covariance and variance
% whrep==10: ROI FCcov scaled

costs = [ 1e-10 1e-9 1e-8 1e-7 1e-6 0.00001, 0.0001, 0.001, 0.01, .1, 1, 10, 100, 1000];

print_dataset = 1;

latex_filename = fullfile(getenv('HOME'), 'Dropbox', 'anatomical', 'Results',...
    sprintf('LR_all_latex_table_regularization_optimized_HCP_truncate%d.txt', truncate));

[ datasetnames] = load_dataset_names();

if strcmp(pred_task, 'task')
    dopredtasks = [1:2];
elseif strcmp(pred_task, 'subj')
    dopredtasks = [3:4];
elseif strcmp(pred_task, 'both')
    dopredtasks=1:4;
end

% credible interval range
p = 0.95;

% accuracy, credible interval, storage
AC = struct();
CR = struct();
CORRECT = struct();

%% Print Headers
f = fopen(latex_filename, 'w');

if strcmp( pred_task, 'task')
    fprintf(f, '\\begin{tabular}{ l l l l l } \n \\hline \\\\[-1ex] ');
    fprintf(f, '& & & \\multicolumn{2}{l}{\\bf Task  }  \\\\ \\cline{4-7} \\\\[-2ex] \n');
    fprintf(f, '{\\bf Type} & {\\bf Feature}  & {\\bf \\# Features} & {\\bf Within } & {\\bf Between }  \\\\ \\\\[-2ex] \n');
elseif strcmp( pred_task, 'subj')
    fprintf(f, '\\begin{tabular}{ l l l l l } \n \\hline \\\\[-1ex] ');
    fprintf(f, '& & & \\multicolumn{2}{l}{\\bf Subject  }  \\\\ \\cline{4-7} \\\\[-2ex] \n');
    fprintf(f, '{\\bf Type} & {\\bf Feature}  & {\\bf \\# Features} & {\\bf Within } & {\\bf Between }  \\\\ \\\\[-2ex] \n');
elseif strcmp( pred_task, 'both')
    fprintf(f, '\\begin{tabular}{ l l l l l l l } \n \\hline \\\\[-1ex] ');
    fprintf(f, '& & & \\multicolumn{2}{l}{\\bf Task  } & \\multicolumn{2}{l}{\\bf Subject Identity } \\\\ \\cline{4-7} \\\\[-2ex] \n');
    fprintf(f, '{\\bf Type} & {\\bf Feature}  & {\\bf \\# Features} & {\\bf Within } & {\\bf Between } &{\\bf Within } & {\\bf Between} \\\\ \\\\[-2ex] \n');
end

for whs = whsets
    
    if print_dataset
        data_set_idx = find([datasetnames{:,2}]== whs);
        if any(strcmp( pred_task, {'task', 'subj'}))
            fprintf(f, '\n \\multicolumn{5}{l}{Dataset %d: %s }\\\\ \\cline{1-5} \\\\ \\\\[-2ex] \n ', ...
                whs, strrep( datasetnames{data_set_idx, 1}, '_', '\_'));
        else
            fprintf(f, '\n \\multicolumn{7}{l}{Dataset %d: %s }\\\\ \\cline{1-7} \\\\ \\\\[-2ex] \n ', ...
                whs, strrep( datasetnames{data_set_idx, 1}, '_', '\_'));
        end
    end
    
    load(fullfile(rd, 'features', sprintf('whs%d_truncate%d.mat', whs, truncate)), 'NT');
    
    %% Print Body
    
    filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_truncate%d_valid' ,...
        1 , whs , 0, 1, truncate) ];
    
    x = load(filenm);
    modelnow = x.modelnow;
    
    % set feature set dimensions
    if size(modelnow.w, 2) == 269
        whrep_feat_no = {'269', '269', '269', '269', '269', '269', '$\frac{269*268}{2}$',...
            '$\frac{269*268}{2}$', '$\frac{269*269}{2}$', '$\frac{269*268}{2}$', '269'};
    elseif size(modelnow.w, 2) == 305
        whrep_feat_no = {'305', '305', '305', '305', '305', '305', '$\frac{305*304}{2}$',...
            '$\frac{305*304}{2}$', '$\frac{305*305}{2}$', '$\frac{305*304}{2}$', '305'};
    end
    
    
    for i = 1:length(whreps)
        whrep = whreps(i);
        
        fprintf(f, '%s & %10.10s & %20.20s', whreps_type{whrep}, strrep( whreps_labels{whrep}, '_', '\_'), whrep_feat_no{whrep});
        
        accuracies = zeros(4,1);
        creds = zeros(4,2);
        corrects = cell(4,1);
        for pt = dopredtasks
            
            
            filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_truncate%d_valid' ,...
                pt , whs , classifiertype , whrep, truncate ) ];
            
            
            % Beta Binomial Hyperparameters set for mode at chance level
            % performance
            if pt < 3
                alpha = 1.25;
                beta = 3;
            else
                alpha = 1.0116;
                beta = 1;
            end
            
            try
                load(filenm);
                
                
                costlabel = strrep( num2str(best_cost), '.', '_');  % costlabels{ bestcostidx( pt, i )};
                fprintf('Loading %s, Best cost = %s\n', filenm, costlabel);
                correct = Y_pred(Y_pred~=0) == YTRUE(Y_pred ~=0);
                accuracy = mean(correct);
                [ cred, pbayesmean, alpha_post, beta_post] = beta_binomial_cred_interval(alpha, beta, correct, p);
                
                
                % convert to percent
                cred = cred*100;
                accuracy = accuracy*100;
                
                fprintf(f, ' & %3.0f (%d, %d)', accuracy, round(cred(1)), round(cred(2)));
                accuracies(pt) = accuracy;
                creds(pt, :) = cred;
                corrects{pt} = correct;
                
            catch
                fprintf('\t FILE NOT FOUND: %s\n', filenm);
                fprintf(f, ' & ~ ');
            end
        end
        
        AC.(sprintf('whs%s', genvarname(num2str(whs)))).(whreps_labels{whrep}) = accuracies;
        CR.(sprintf('whs%s', genvarname(num2str(whs)))).(whreps_labels{whrep}) = creds;
        CORRECT.(sprintf('whs%s', genvarname(num2str(whs)))).(whreps_labels{whrep}) = corrects;
        
        
        fprintf(f, ' \\\\ \n ');
    end
    
    fprintf(f, ' \\\\[-1.5ex]');
    
    
    
    
end

fprintf(f, '\\\\[-1.5ex] \\hline \n \\end{tabular}');

% Save AC

filenm = fullfile( rd, 'Results', sprintf('AC_c%d_truncate%d', classifiertype, truncate));
save(filenm, 'AC', 'CR', 'CORRECT' );



%% Print Best Costs

best_costs = zeros( length(whsets), length(whreps), length(dopredtasks));
for j = 1:length(whsets)
    whs = whsets(j);
    
    for i = 1:length(whreps)
        whrep = whreps(i);
        
        for pt = dopredtasks
            
            
            %                  filenm = [ 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_cost%s' ,...
            %                     pt , whs , classifiertype , whrep, costlabel ) ];
            filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_truncate%d_valid' ,...
                pt , whs , classifiertype , whrep, truncate ) ];
            
            
            try
                load(filenm);
                costlabel = strrep( num2str(best_cost), '.', '_');  % costlabels{ bestcostidx( pt, i )};
                fprintf('Loading %s, Best cost = %s\n', filenm, costlabel);
                
                best_costs(whs, i, pt) = best_cost;
            catch
                continue;
            end
            
        end
    end
end

% print out latex table 
f = fopen(latex_filename, 'r'); 
line = fgets(f);
while line ~= -1
    fprintf('%s', line); 
    line = fgets(f); 
end

%% Plot best costs for each representation

py = zeros( length(whreps) , length(costs) );
count = 1;
for i = 1:length(whreps)
    temp = squeeze( best_costs(:, i, :) );
    ps = hist( temp(:), costs );
    py( count , : ) = ps * 100 / sum( ps );
    count = count + 1;
end


% print table
ii = sum(py, 1) > 0;
py = [costs; py];
py = py(:,ii);
tab = table( ['costs'; whreps_labels(whreps)'], 'VariableName', {'whrep'});
tab = [tab, array2table( py)]

end

