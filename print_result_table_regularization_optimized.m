function [ AC , CR, CORRECT, best_costs, tab ] = print_result_table_regularization_optimized( classifiertype, whsets, isHPC)
[wd, rd] = set_directories(isHPC);


if classifiertype ~= 8
pred_task = 'both'; 
else 
   pred_task = 'task';  
end

whreps = [1,4,7:9]; % 10]; % [1,2,4,9,10]
whreps_labels = { 'BVSD', 'BVM', 'BVVS', 'BVV', 'MSSD', 'SQRT_MSSD', 'FCP', 'FCC', 'FCCV', 'FCCS'};
whreps_type = {'BV', 'BV', 'BV', 'BV', 'BV', 'BV', 'FC', 'FC', 'FC', 'FC' };

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

% FOR 1NN
%classifiertypes = 8;
%filename = 'Results/nearest_neighbor_latex_table.txt';

% FOR LR
% classifiertypes = 0;

if classifiertype == 8
    filename = fullfile(getenv('HOME'), 'Dropbox', 'anatomical', 'Results',...
        '1NN_all_latex_table_regularization_optimized_HCP.txt');
else
    filename = fullfile(getenv('HOME'), 'Dropbox', 'anatomical', 'Results',...
        'LR_all_latex_table_regularization_optimized_HCP.txt');
end

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
f = fopen(filename, 'w');

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
        if any(strcmp( pred_task, {'task', 'subj'}))
            fprintf(f, '\n \\multicolumn{5}{l}{Dataset %d: %s }\\\\ \\cline{1-5} \\\\ \\\\[-2ex] \n', whs, strrep( datasetnames{whs}, '_', '\_'));
        else
            fprintf(f, '\n \\multicolumn{7}{l}{Dataset %d: %s }\\\\ \\cline{1-7} \\\\ \\\\[-2ex] \n', whs, strrep( datasetnames{whs}, '_', '\_'));
        end
    end
    
    load(fullfile(rd, 'features', sprintf('whs%d.mat', whs)), 'NT');
    
    %% Print Body
    
    
    filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_valid' ,...
        1 , whs , 0, 1) ];
    
    x = load(filenm);
    modelnow = x.modelnow;
    
    % set feature set dimensions
    if size(modelnow.w, 2) == 269
        whrep_feat_no = {'269', '269', '269', '269', '269', '269', '$\frac{269*268}{2}$',...
            '$\frac{269*268}{2}$', '$\frac{269*269}{2}$', '$\frac{269*268}{2}$', };
    elseif size(modelnow.w, 2) == 305
        whrep_feat_no = {'305', '305', '305', '305', '305', '305', '$\frac{305*304}{2}$',...
            '$\frac{305*304}{2}$', '$\frac{305*305}{2}$', '$\frac{305*304}{2}$', };
    end
    
    
    for i = 1:length(whreps)
        whrep = whreps(i);
        
        fprintf(f, '%s & %s & %s', whreps_type{whrep}, strrep( whreps_labels{whrep}, '_', '\_'), whrep_feat_no{whrep});
        
        accuracies = zeros(4,1);
        creds = zeros(4,2); 
        corrects = cell(4,1); 
        for pt = dopredtasks
            
            if classifiertype == 8
                filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d' ,...
                    pt , whs , classifiertype , whrep ) ];
            else
                filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_valid' ,...
                    pt , whs , classifiertype , whrep ) ];
            end
            
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
                
                if classifiertype == 8
                    %% overall accuracies
                    if pt == 1 || pt ==2
                        error('Task classification not supported');
                    elseif pt == 3
                        W = W(eye(NT) == 0);
                        W = [W{:}];
                        correct = W(:);
                        accuracy = mean(correct);
                        [ cred, pbayesmean, alpha_post, beta_post] = beta_binomial_cred_interval(alpha, beta, correct, p);
                        cred = cred;
                    elseif pt == 4
                        B2 = [B{:}];
                        correct = B2(:);
                        accuracy = mean(correct);
                        [ cred, pbayesmean, alpha_post, beta_post] = beta_binomial_cred_interval(alpha, beta, correct, p);
                        cred = cred;
                    end
                    
                else
                    costlabel = strrep( num2str(best_cost), '.', '_');  % costlabels{ bestcostidx( pt, i )};
                    fprintf('Loading %s, Best cost = %s\n', filenm, costlabel);
                    correct = Y_pred(Y_pred~=0) == YTRUE(Y_pred ~=0);
                    accuracy = mean(correct);
                    [ cred, pbayesmean, alpha_post, beta_post] = beta_binomial_cred_interval(alpha, beta, correct, p);
                end
                
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
        
        AC.(sprintf('whset%d', whs)).(whreps_labels{whrep}) = accuracies;
        CR.(sprintf('whset%d', whs)).(whreps_labels{whrep}) = creds;
        CORRECT.(sprintf('whset%d', whs)).(whreps_labels{whrep}) = corrects;
            
        
        fprintf(f, ' \\\\ \n ');
    end
    
    fprintf(f, ' \\\\[-1.5ex]');
    
    
    
    
end

fprintf(f, '\\\\[-1.5ex] \\hline \n \\end{tabular}');

% Save AC 
if classifiertype ~= 8
    filenm = fullfile( rd, 'Results', sprintf('AC_c%d', classifiertype)); 
    save(filenm, 'AC'); 
end



%% Print Best Costs
if classifiertype ~= 8
    
    best_costs = zeros( length(whsets), length(whreps), length(dopredtasks));
    for j = 1:length(whsets)
        whs = whsets(j);
        
        for i = 1:length(whreps)
            whrep = whreps(i);
            
            for pt = dopredtasks
                
                
                %                  filenm = [ 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_cost%s' ,...
                %                     pt , whs , classifiertype , whrep, costlabel ) ];
                filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_valid' ,...
                    pt , whs , classifiertype , whrep ) ];
                
                
                try
                    load(filenm);
                    costlabel = strrep( num2str(best_cost), '.', '_');  % costlabels{ bestcostidx( pt, i )};
                    fprintf('Loading %s, Best cost = %s\n', filenm, costlabel);
                    
                    best_costs(whs, whrep, pt) = best_cost;
                catch
                    continue;
                end
                
            end
        end
    end
else
    best_costs = NaN; tab =NaN;
end

%% Plot best costs for each representation

if classifiertype ~=8
    
    py = zeros( length(whreps) , length(costs) );
    count = 1;
    for whrep = whreps
        temp = squeeze( best_costs(:, whrep, :) );
        temp = temp( whsets, :);
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
end

