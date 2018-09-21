function [] = plot_feature_selection_acc(isHPC, whs, whreps, line_smoother)

% line_smoother = 1
[wd, rd] = set_directories(isHPC); % working directory and results directory
classifiertype=6; %6=L1
fontsize = 12;

featurelabels = { ...
    'BVSD' , ...
    'M' , ...
    'BVVS',...
    'BV',...
    'MSSD', ...
    'SQRT_MSSD', ...
    'FCP' , ...
    'FCC', ...
    'FCCV', ...
    'FCCVS'
    };

map = brewermap( length(whreps), 'Set1');

h0 = figure(1); clf;
h0.Position = [146 414 1284 564]

for pt = 1:2
    s = subplot(1,2,pt);
    min_max_nonzeros = 0;
    for i = 1:length(whreps)
        whrep = whreps(i);
        filenm = [ rd filesep 'feature_selection' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_valid' ,...
            pt , whs , classifiertype , whrep ) ];
        
        try 
            dat = load(filenm); % acc and nonzeros are of dimension [Nfold x Ncost]
        catch 
            error( sprintf(['Analysis not run for: \n\tpt=%d \n\twhs=%d \n\tclassifiertype=%d', ...
                '\n\twhrep=%d\nRun the file run_all_feature_selection.m with', ...;
                ' the appropriate options'], ...
                pt, whs, classifiertype, whrep));  
        end
        acc = dat.acc;
        nonzeros = dat.nonzeros;
        
        fprintf( 
       
        
        
        mean_nonzeros = cellfun( @(x) sum(x), nonzeros);
        m = max(mean_nonzeros(:));
        if min_max_nonzeros==0
            min_max_nonzeros = m;
        elseif m < min_max_nonzeros
            min_max_nonzers=m;
        end
        
        % average over folds
        x = mean(mean_nonzeros,1);
        %x = log(mean(mean_nonzeros,1));
        y = mean( acc, 1 );
        scatter( x, y, [], map(i,:) ); hold on;
        
        % lowess smoother
        % this loop will run a lowess smoother with the smallest smoothing
        % window possible that doesn't lead to NaN results or estimation
        % problems 
        if line_smoother
            f = 0.15;
            count = 0; 
            while ~count
                try
                    f = f + 0.01;
                    datain = [x;y]';
                    [dataout lowerLimit upperLimit xy] = lowess(datain,f,0);
                    plot(dataout(:,1), dataout(:,3), '-', 'Color', map(i,:) ); hold on;
                    if any(isnan( dataout(:,3)))
                       error('nan values');  
                    end
                    count = 1;
                    fprintf('f = %2.2f\n', f); 
                catch
                    warning(sprintf( 'problem with lowess smoother: f = %d', f));
                end
            end
        end
        % max accuracy
        max_acc = max(y);
        h(i) = plot( [min(x), max(x)], [max_acc, max_acc], '--', 'Color', map(i,:) ); hold on;
        
    end
    xlabel('# nonzeros');
    ylabel('accuracy');
    xlim( [0, min_max_nonzeros]);
    ylim( [0, 100]);
    s.FontSize = fontsize;
    
    if pt==1
        title(sprintf('Within-Session Task Classification\nNumber of nonzero features vs. accuracy'));
    else
        title(sprintf('Between-Session Task Classification\nNumber of nonzero features vs. accuracy'));
    end
    
    legend( h, featurelabels(whreps), 'Location', 'best' );
    
end

rep_label = strjoin( featurelabels(whreps), '_');
rep_label = strrep( rep_label, ' ROIs', '');
rep_label = strip(rep_label);
filenm = fullfile(rd, 'images', sprintf( 'feature_selection_%s', rep_label));

print(filenm, '-depsc');
