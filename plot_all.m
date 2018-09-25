% This file will create plots for each type of analysis. Analysis can be
% set using the plot_mode variable

%% Options
plot_mode = "processing_robustness"; % feature_selection, data_exploration, task_classification, subject_classification, processing_robustness
truncate = 0;
isHPC = 0;

% subject_classification options 

whs = 0;
whreps = [1,7];

% whs==-1: old dataset (old OSU preprocessing pipeline)
% whs== 0: baseline ( HCP preprocessing pipeline )
% whs== 1: dataset 1 ( HCP preprocessing pipeline )
% whs== 2: dataset 2 ( HCP preprocessing pipeline )
% whs== 3: dataset 3 ( HCP preprocessing pipeline )

% whrep==1: BVSD
% whrep==2: M
% whrep==3: BVVS
% whrep==4: BVV
% whrep==5: MSSD
% whrep==6: SQRT MSSD
% whrep==7: FCP
% whrep==8: FCC
% whrep==9: FCCV
% whrep==10: FCCS

if strcmp( plot_mode, 'feature_selection')
    
    % options 
    line_smoother = 1; % do lowess smoothing for the feature selection accuracy plot?
    
    plot_feature_selection_acc(isHPC, whs, whreps, line_smoother)
    
elseif strcmp( plot_mode, 'data_exploration')
    
    % options 
    roi_sort_by = 'Lobe'; % 'Lobe'; % 'Lobe' or 'Hemishpere'
    fontsize = 20;
    labels_on = 0;
    normalize = 0; % 0 is no normalization, 1 is by subject, 2 is by ROI, 3 is globally
    grid_on = 1;
    remove_max = 0;
    
    plot_BV_heatmap (whs, truncate, roi_sort_by, fontsize, labels_on, ...
        normalize, grid_on, remove_max) 
    
    fontsize = 18; 
    plot_BV_rest_vs_task( whs, truncate, fontsize); 
    
elseif strcmp( plot_mode, 'task_classification')

     % options 
    fontsize = 15;
    fontcolor = 'k';
    usepatchimg = 0;
    include_zeros = 0;
    zero_char = '-';
    facealpha = 1;
    classifiertype = 0; % same as liblinear 0=L2, 6=L1, we added 8 to represent 1NN 
    
    plot_confusion_matrices( isHPC, whreps, whs, truncate, classifiertype, ...
    fontsize, fontcolor, usepatchimg, include_zeros, zero_char, facealpha)
    
elseif strcmp( plot_mode, 'subject_classifications')
    
    plot_1NN_heatmaps( isHPC, whs )
    
elseif strcmp( plot_mode, 'processing_robustness')
    
     % task classification plots 
     
     % run aggregation of results 
     [ AC , CR, CORRECT, best_costs, tab ] = print_result_table_regularization_optimized( classifiertype, whsets, isHPC, truncate); 
     
     % plot 
     plot_LR_all_accuacies(isHPC, classifiertype, truncate); 
     
     % subject identification  
     
     % run aggregation of results 
     finn_et_al_analysis(isHPC, whs, truncate);  
     
     % plot
     plot_finn_all_accuracies( isHPC, truncate ); 
     
end


