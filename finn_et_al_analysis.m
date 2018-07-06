
classifiertype = 8;
print_dataset = 1; 
isHPC = 0; 
whs = 26; 

[wd, rd] = set_directories(isHPC); % working directory and results directory

whsets = [26, 22, 23, 27]; 
whreps = [1,4, 7:9]; % 10]; % [1,2,4,9,10]
whreps_labels = { 'BVSD', 'BVM', 'BVVS', 'BVV', 'MSSD', 'SQRT_MSSD', 'FCP', 'FCC', 'FCCV', 'FCCS'};
whreps_type = {'BV', 'BV', 'BV', 'BV', 'BV', 'BV', 'FC', 'FC', 'FC', 'FC' };

load(fullfile(rd, 'features', sprintf('whs%d.mat', whs)), 'NR');
load( fullfile(rd, 'ICAresults','tasklist') );
NT = length(tasks); 

if NR== 269
    whrep_feat_no = {'269', '269', '269', '269', '269', '269', '$\frac{269*268}{2}$',...
        '$\frac{269*268}{2}$', '$\frac{269*269}{2}$', '$\frac{269*268}{2}$', };
elseif size(modelnow.w, 2) == 305
    whrep_feat_no = {'305', '305', '305', '305', '305', '305', '$\frac{305*304}{2}$',...
        '$\frac{305*304}{2}$', '$\frac{305*305}{2}$', '$\frac{305*304}{2}$', };
end

[ datasetnames] = load_dataset_names();

AC = struct(); 
CR = struct(); 
CORRECT = struct(); 
accuracy_labels = {'Within', 'Between', 'Task-Same Task', 'Task-Different Task', 'Rest-Rest'}; 

%% PRINT LATEX TABLE
% beta-binomial hyper parameters
alpha = 1.0116; % mode is 1/174, but with pretty smooth density
beta = 1;

p = 0.95;

% prediction tasks
dopredtasks = [3,4];

% load cost info
filename = fullfile(rd, 'Results', 'finn_et_al_table.txt');

f = fopen(filename, 'w');

% Print Headers
fprintf(f, '\\begin{tabular}{ l l l l l l l l } \n \\hline \\\\[-1ex] ');
fprintf(f, '{\\bf Type} & {\\bf Method}  & {\\bf \\# Feat} & {\\bf Within } & \\multicolumn{4}{c}{\\bf Between } \\\\ \\\\[-2ex] \\cline{5-8}\n');
fprintf(f, ' & & & & {\\bf All} & {\\bf Same Task} & {\\bf Diff Task} & {\\bf Rest} \\\\ \\\\[-2ex] \n');

% Print Body

for whs = whsets
    
    if print_dataset
            fprintf(f, '\n \\multicolumn{8}{l}{Dataset %d: %s }\\\\ \\cline{1-7} \\\\ \\\\[-2ex] \n', whs, strrep( datasetnames{whs}, '_', '\_'));
    end
    
    for i = 1:length(whreps)
        whrep = whreps(i);
        fprintf(f, '%s & %s & %s ', whreps_type{i}, strrep(whreps_labels{whrep}, '_', ' '), whrep_feat_no{i});
        
        accuracies = zeros(5,1);
        creds = zeros(5,2); 
        corrects = cell(5,1); 
        for pt = dopredtasks
            
            filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d' , pt , whs , classifiertype , whrep ) ];
            
            load(filenm); % loads W and B
            
            %% Calculate accuracies
            
            % overall accuracies
            if pt == 1 || pt ==2
                error('Task classification not supported');
            elseif pt == 3
                W = W(eye(NT) == 0);
                W = [W{:}];
                correct = W(:);
                accuracy = 100*mean(correct);
                [ cred, pbayesmean, alpha_post, beta_post] = beta_binomial_cred_interval(alpha, beta, correct, p);
                cred = cred*100; 
                
                accuracies(1) = accuracy; 
                creds(1,:) = cred; 
                corrects{1} = correct; 
            elseif pt == 4
                B2 = [B{:}];
                correct = B2(:); 
                accuracy = 100*mean(correct); 
                [ cred, pbayesmean, alpha_post, beta_post] = beta_binomial_cred_interval(alpha, beta, correct, p);
                cred = 100*cred; 
                
                accuracies(2) = accuracy; 
                creds(2,:) = cred; 
                corrects{2} = correct; 
            end
            
            
            if pt == 4
                % credible interval for type of pairs of task scans
                K = size(W, 1);  % number of tasks
                
                rest_idx = strcmp(tasks , 'Resting');
                r2r = B(rest_idx, rest_idx);
                t2st= B(eye(NT)==1);
                t2st(rest_idx) = []; 
                t2dt = B( ~eye(NT));
                
                correct = struct();
                correct.r2r = r2r{:};
                correct.t2st = [ t2st{:}];
                correct.t2dt = [ t2dt{:}];
                
                fun = @(x) round( beta_binomial_cred_interval(alpha,beta,x,0.95)*100, 0);
                acc_cred  = structfun( @(x) fun(x(:)), correct, 'UniformOutput', false);
                acc_mean  = structfun( @(x) round(100*mean(x(:)),0), correct, 'UniformOutput', false);
                
                
                fprintf(f, ' & %3.0f (%d, %d) & %3.0f (%d, %d) & %3.0f (%d, %d) & %3.0f (%d, %d)', ...
                    accuracy, round(cred(1)), round(cred(2)), ...
                    acc_mean.t2st, acc_cred.t2st(1), acc_cred.t2st(2), ...
                    acc_mean.t2dt, acc_cred.t2dt(1), acc_cred.t2dt(2), ...
                    acc_mean.r2r, acc_cred.r2r(1), acc_cred.r2r(2) );
                
                accuracies(3) = acc_mean.t2st; 
                accuracies(4) = acc_mean.t2dt; 
                accuracies(5) = acc_mean.r2r; 
                
                creds(3,:) = acc_cred.t2st;
                creds(4,:) = acc_cred.t2dt;
                creds(5,:) = acc_cred.r2r;
                
                corrects{3} = correct.t2st; 
                corrects{4} = correct.t2dt; 
                corrects{5} = correct.r2r; 
            else
                fprintf(f, ' & %3.0f (%d, %d)', ...
                    accuracy, round(cred(1)), round(cred(2)));
            end
       
        end
        
        AC.(sprintf('whset%d', whs)).(whreps_labels{whrep}) = accuracies;
        CR.(sprintf('whset%d', whs)).(whreps_labels{whrep}) = creds;
        CORRECT.(sprintf('whset%d', whs)).(whreps_labels{whrep}) = corrects;
        
        fprintf(f, ' \\\\ \n ');
    end
    
end

fprintf(f, '\\\\[-1.5ex] \\hline \n \\end{tabular}');

filenm = fullfile( rd, 'Results', sprintf('AC_c%d', classifiertype)); 
save( filenm, 'AC'); 





