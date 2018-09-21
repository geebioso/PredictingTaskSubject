%% Run Classifier on ICA coefficients V4
% Only train on session 1 data

%%% This file will is the same as runclassifier_ICAdata_v5.m except that
%%% for within session classification, it will test on only the 19 subjects
%%% that were scanned in both sessions

% Feb 2 2017
clear D;
isHPC = 0;
truncate = 1; 

[wd, rd] = set_directories(isHPC); % working directory and results directory

whsets = [-1:3];     % which data set? 
runIC = 0;   % run models with ICA results?

dopredtasks    = 3:4;%[1:5];
% pt == 1, Predict Task Session 1
% pt == 2, Predict Task Session 2
% pt == 3, Predict Subject Session 1
% pt == 4, Predict Subject Session 2
% pt == 5, Predict Task Session 1 All Subject

K                 = 5; % number of folds, use 5 so that stratified cross validation tests points from all 19 subjects

dorepresents      = [1:10 ];% [ 1,2,4,9:12 ];

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

featurelabels = { ...
    'StdDev ROIs' , ...
    'Means  ROIs' , ...
    'ROI FV scaled',...
    'ROI FV',...
    'MSSD', ...
    'SQRT_MSS', ...
    'Corrs  ROIs' , ...
    'Covs   ROIS', ...
    'Covs+Vars ROIs', ...
    'Covs   ROIS scaled'
    };

classifiertype=8;
seed              = 1;
rng( seed );

%%
predtasklabels = { ...
    'Task Classifier Session 1' , ...
    'Task Classifier Session 2' , ...
    'Subject Identity Classifier Session 1' , ...
    'Subject Identity Classifier Session 2', ...
    'Task Classifier Session 1 All Subjects'};

%% Load the task names
if ~exist( 'tasks' , 'var' )
    load( fullfile(rd, 'ICAresults','tasklist') );
end

for whs = whsets
    
    fprintf('Dataset %d\n', whs); 
    %% Load the Features
    load( fullfile(rd, 'tasklist.mat') );
    
    try
        filenm = fullfile( rd, 'features', sprintf('whs%d_truncate%d.mat', whs, truncate));
        fprintf('\tLoading %s\n', filenm);
        load( filenm );
        %                     NT: 9
        %                     NR: 269
        %                  combs: [1827×3 double]
        %                     NC: 1827
        %     subsinbothsessions: [19×1 double]
        %               features: [1×1 struct]
        
        % combs columns are [subject, task, session]
    catch
        error('Features are not computed for dataset %d', whs);
    end
    
    % compute extra features
    roi_vars_scale = features.stds.^2;
    features.vars_scale = (roi_vars_scale - mean(roi_vars_scale(:)))/std(roi_vars_scale(:));
    features.covs_scale = (features.covs - mean(features.covs(:)))/(std(features.covs(:)));
    clear roi_vars_scale
    
    
    %% Classification loop
    for pt=dopredtasks
        fprintf( '\t**************************************************************\n' );
        fprintf( '\tPrediction Task: %s\n' , predtasklabels{ pt });
        fprintf( '\t**************************************************************\n' );
         
        for whrep=dorepresents
            fprintf( '\tFeatures:  %s\n' , featurelabels{ whrep } );
            
            if whrep==1
                X  = features.stds;
            end
            
            if whrep==2
                X  = features.means;
                
                % Standardize?
                X = X -  repmat( mean( X , 1 ) , size(X,1) , 1 );
                X = X ./ repmat( std( X , [] , 1 ) , size(X,1) , 1 );
            end
            
            if whrep==3
                X = features.vars_scale;
            end
            
            if whrep==4
                X = features.stds.^2;
            end
            
            if whrep==5
                X = features.mssd;
            end
            
            if whrep==6
                X = sqrt(features.mssd);
            end
            
            if whrep==7
                X  = features.corrs;
            end
            
            if whrep==8
                X = features.covs;
            end
            
            if whrep==9
                X = [features.covs, features.stds.^2];
            end
            
            if whrep==10
                X = features.covs_scale;
            end
            
            X = double(X);
            % number of features
            NF = size( X , 2 );           
            
%             %% Parallel loop
%             Y_pred_cell = cell( nsets , 1 );
%             YP_pred_cell = cell( nsets , 1 );
%             models = cell( nsets , 1 );
%             %parfor j=1:K
            
            %% Train/Test
            
            % precompute all similarities 
            C = round(corr( X', X' ), 6);
            
            % storage 
            W = cell(NT); 
            B = cell(NT); 
            
            % session indices 
            sess1_idx = combs(:,3) == 1; 
            sess2_idx = and( combs(:,3) == 2, ismember( combs( : , 1) , subsinbothsessions )); 
            
            % Classification loop 
            for t1 = 1:NT
                t1_idx = combs(:,2) == t1; 
                for t2 = 1:NT
                    t2_idx = combs(:,2) == t2; 
                    
                    %% Within 
                    % train/test indices 
                    test_idx = and(sess1_idx, t1_idx);
                    train_idx = and(sess1_idx, t2_idx);
                    
                    % subject labels 
                    test_labels = combs(test_idx,1); 
                    train_labels = combs(train_idx,1); 
                    
                    % reduce similarity matrix and find most similar 
                    Cnow = C( test_idx, train_idx); 
                    [~, most_similar] = max(Cnow, [], 2);
                    
                    W{t1,t2} = test_labels == train_labels(most_similar);
                    
                    %% Between 
                    % train/test indices
                    test_idx = and(sess2_idx, t1_idx);
                    train_idx = and(sess1_idx, t2_idx);
                    
                    % subject labels
                    test_labels = combs(test_idx,1);
                    train_labels = combs(train_idx,1);
                    
                    % reduce similarity matrix and find most similar
                    Cnow = C( test_idx, train_idx);
                    [~, most_similar] = max(Cnow, [], 2);
                    
                    B{t1,t2} = test_labels == train_labels(most_similar);
                    
                end
            end
            
            % And save the resulting model trained on all of session 1 data
            filenm = [ rd filesep 'discrweights' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_truncate%d' ,...
                pt , whs , classifiertype , whrep, truncate ) ];
            save( filenm , 'W' , 'B', 'whs' , 'classifiertype', 'whrep' , 'featurelabels'); 
            
        end
    end
end
