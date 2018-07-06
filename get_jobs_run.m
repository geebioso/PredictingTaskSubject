% function [] = get_jobs_run( isHPC ) 

[wd, rd] = set_directories(isHPC); % working directory and results directory 

files = dir([rd filesep 'discrweights']); 

NF = length(files);
filenames = cell(NF,1);

for f = 1:NF
    filenames{f} = files(f).name;
end

% get rid of hidden files 
for f = 1:NF
    filenames{f} = files(f).name; 
    isvalid = strfind(filenames{f}, 'valid'); 
    if filenames{f}(1) == '.'
      filenames{f} = [];  
    elseif filenames{f}(1:2) ~= 'pt' % take only liblinear classifier results 
        filenames{f} = [];  
    elseif isempty(isvalid)% make sure we only take the validated models 
        filenames{f} = []; 
    end
end

ii = cellfun( @(x) ~isempty(x), filenames, 'UniformOutput', true);
filenames = filenames(ii); 
NF = length(filenames); 

%% Split and get information to store in table 
info = zeros(NF, 4); 
for f = 1:NF
    filenm = filenames{f};
    name_vec = split(filenm, '_');
    
    %     {'pt101'    }
    %     {'whs10'    }
    %     {'c0'       }
    %     {'whr10'    }
    %     {'valid.mat'}
    
    pt = regexp(name_vec{1},'\d*','Match');
    pt = str2num( pt{1} );
    
    whs = regexp(name_vec{2},'\d*','Match');
    whs = str2num( whs{1} );
    
    c = regexp(name_vec{3},'\d*','Match');
    c = str2num( c{1} );
    
    whrep = regexp(name_vec{4},'\d*','Match');
    whrep = str2num( whrep{1} );
    
    info(f,:) = [pt, whs, c, whrep]; 
    
end

tab = table( info(:,2), info(:,1), info(:,4), info(:,3), 'VariableNames', {'whs', 'pt', 'whrep', 'c'}); 
ii = tab.pt < 11; 
tab = tab(ii,:); 
tab = sortrows(tab); 

