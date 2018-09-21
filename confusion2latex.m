function [ ] = confusion2latex( C, accuracies, filename, tasknames, include_zeros, zero_char, whrep_labels)
% puts the confusion matrix into format for latex code

N = length(C);  % number of confusion matrices
T = length(tasknames);

for i = 1:N
    if i > 1
        f = fopen(filename, 'a');
    else
        f = fopen(filename, 'w');
    end
    
    fprintf(f, '\\multicolumn{8}{l}{\\bf %s } \\\\ \n', whrep_labels{i}); 
    for t = 1:T
        fprintf(f, '& {\\bf %s}', char(tasknames{t})); 
    end
    fprintf(f,'\\\\\n'); % prints \\ newline 
    
    
    for t = 1:T
        fprintf(f, '{\\bf %s} ', char(tasknames{t}));
        for j = 1:T
            if include_zeros
                fprintf(f, '& %d ', C{i}(t,j));
            else
                if C{i}(t,j) == 0
                    fprintf(f, '&  %s', zero_char);
                else
                   fprintf(f, '& %d ', C{i}(t,j)); 
                end
            end
        end
        if t == T
            fprintf(f, '\\\\ \\\\[-2ex] \\hline \\\\ \n\n');  
        else
            fprintf(f,'\\\\ \n');
        end
    end
    
    fprintf(f,'{ \\bf Accuracy} ');
    
    for t = 1:T
        fprintf(f, '& %2.2f ', accuracies{i}(t));
    end
    fprintf( f, '\\\\ \\\\[-2ex] \\hline \\\\ \n\n'); 
    
    fclose(f);
    
end

end
