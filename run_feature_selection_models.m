function [] = run_feature_selection_models(isHPC, whs, whrep, pt)
%%% This file will is the same as ../runclassifier_ICAdata_v7.m except that
%%% is will run on the UCI HPC
% Must do L1 regularization for variable selection 

% isHPC = 0; whs = 26; whrep = 1;

% March 26 2018

% 8GB for features of size 269
%

if ~isa(isHPC, 'numeric')
    isHPC = str2num(isHPC);
    whs = str2num(whs);
    whrep = str2num(whrep);
    pt = str2num(pt); 
end

fprintf('isHPC=%d, whs=%d, whrep=%d, pt=%d\n', isHPC, whs, whrep, pt);

%% Options
dosave = 1; % save models?
K                 = 2; % number of folds, use 5 so that stratified cross validation tests points from all 19 subjects
classifiertype    = 6; % 0: L2 regression, 6:L1 regression or 7: L2 regression dual

%costs = [ 1e-10 1e-9 1e-8 1e-7 1e-6 0.00001, 0.0001, 0.001, 0.01, .1, 1, 10, 100, 1000];
%costs = 2.^(-20:10);
costs = 1.1.^(-200:73); 

NCO = length(costs);
costlabels = cell(NCO,1);
for c = 1:NCO
    costlabels{c} = strrep( num2str(costs(c)), '.', '_');
end
seed              = 1;

[wd, rd] = set_directories(isHPC); % working directory and results directory

% Classifier Type
% 	 0 -- L2-regularized logistic regression (primal)
% 	 1 -- L2-regularized L2-loss support vector classification (dual)
% 	 2 -- L2-regularized L2-loss support vector classification (primal)
% 	 3 -- L2-regularized L1-loss support vector classification (dual)
% 	 4 -- support vector classification by Crammer and Singer
% 	 5 -- L1-regularized L2-loss support vector classification
% 	 6 -- L1-regularized logistic regression
% 	 7 -- L2-regularized logistic regression (dual)

if isHPC
    fprintf('Using mex files in current directory\n');
else
    if ispc
        addpath(fullfile(getenv('HOME'), 'Dropbox','FMRI','Projects','liblinear-2.11', 'windows'));
    else
        addpath(fullfile(getenv('HOME'), 'Dropbox','FMRI','Projects','liblinear-2.11', 'matlab'));
    end
end
%%

% pt == 1, Predict Task Session 1
% pt == 2, Predict Task Session 2
% pt == 3, Predict Subject Session 1
% pt == 4, Predict Subject Session 2
% pt == 5, Predict Task Session 1 All Subject

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

rng( seed );

%%
predtasklabels = { ...
    'Task Classifier Session 1' , ...
    'Task Classifier Session 2' , ...
    'Subject Identity Classifier Session 1' , ...
    'Subject Identity Classifier Session 2', ...
    'Task Classifier Session 1 All Subjects'};

%% Load the Features
load( fullfile(rd, 'tasklist.mat') );

try
    filenm = fullfile( rd, 'features', sprintf('whs%d.mat', whs));
    fprintf('Loading %s\n', filenm);
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
fprintf( '**************************************************************\n' );
fprintf( 'Prediction Task: %s\n' , predtasklabels{ pt });
fprintf( '**************************************************************\n' );

if pt==1
    % Predict Task Session 1
    whok = find( combs( : , 3 ) == 1 ); % session 1 subjects
    Y =  combs( whok , 2 ); % task ids
elseif pt==2
    % Predict Task Session 2
    % Include all session 1 data and only session 2 data from subjects who are in session 1
    whok = find( combs( : , 3 ) == 1 | ismember( combs( : , 1 ) , subsinbothsessions ));
    Y =  combs( whok , 2 );
elseif pt==3
    % Predict Subject Session 1
    whok = find( combs( : , 3 ) == 1 ); % session 1 subjects
    Y =  combs( whok , 1 );
    % Remap subject indices to contiguous set
    OLDY = Y;
    [ temp , ~ , Y ] = unique( Y );
elseif pt==4
    % Predict Subject Session 2
    % Include all session 1 data and only session 2 data from subjects who are in session 1
    whok = find( combs( : , 3 ) == 1 | ismember( combs( : , 1 ) , subsinbothsessions ));
    Y =  combs( whok , 1 );
    % Remap subject indices to contiguous set
    OLDY = Y;
    [ temp , ~ , Y ] = unique( Y );
elseif pt==5
    % Predict Task Session 1
    whok = find( combs( : , 3 ) == 1 ); % session 1 subjects
    Y =  combs( whok , 2 ); % task ids
end
NY = length( Y );

if pt==1 | pt==3
    % Stratified CV partition but only on session 1 data from 19 subjects who participated
    % in both sessions
    insession2 = ismember(combs(whok,1), subsinbothsessions); % session 1 subjects that are also in session 2
    always_train = and(combs(whok,3) == 1, ~insession2) ;
    insession2 = find(insession2); % convert to numerical indices to index into training split later
    nsets = K;
    whtrain_sets = cell( nsets , 1 );
    whtest_sets  = cell( nsets , 1 );
    cv = cvpartition( combs(insession2, 1),'kfold',K);
    for j = 1:nsets
        if (j<=K) % test on stratified CV partition of 19 subjects
            whtrain = always_train;
            whtrain(insession2(cv.training(j))) = 1;
            whtest  = false(size(always_train));
            whtest(insession2(cv.test( j ))) = 1;
        else % test on training data
            whtrain = 1:length(Y);
            whtest  = 1:length(Y);
        end
        whtrain_sets{ j } = whtrain;
        whtest_sets{ j } = whtest;
    end
    
elseif pt==5
    % CV partitions when generalizing to session 1 data
    nsets = K;
    whtrain_sets = cell( nsets , 1 );
    whtest_sets  = cell( nsets , 1 );
    cv = cvpartition( length( Y ) ,'kfold',K);
    for j=1:nsets+1
        if (j<=K)
            whtrain = find( cv.training( j ));
            whtest  = find( cv.test( j ));
        else
            whtrain = 1:length(Y);
            whtest  = 1:length(Y);
        end
        whtrain_sets{ j } = whtrain;
        whtest_sets{ j } = whtest;
    end
elseif pt == 4
    % Training and Test set when generalizing to session 2 data
    nsets = 1;
    whtrain_sets = cell( nsets , 1 );
    whtest_sets  = cell( nsets , 1 );
    
    % Train on all session 1 data
    whtrain = find( combs( whok , 3 ) == 1 );
    % Test on all session 2 data
    whtest = find( combs( whok , 3 ) == 2 );
    
    whtrain_sets{ 1 } = whtrain;
    whtest_sets{ 1 } = whtest;
elseif pt == 2
    % Training and Test set when generalizing to session 2 data
    nsets = 1;
    whtrain_sets = cell( nsets , 1 );
    whtest_sets  = cell( nsets , 1 );
    
    % Train on all session 1 data
    whtrain = find( and( combs( whok , 3 ) == 1, ~ismember(combs(whok,1), subsinbothsessions ) ));
    %whtrain = find( combs( whok , 3 ) == 1); 
    % Test on all session 2 data
    whtest = find( combs( whok , 3 ) == 2 );
    
    whtrain_sets{ 1 } = whtrain;
    whtest_sets{ 1 } = whtest;
end


fprintf( 'whrep = %d\n', whrep);
fprintf( 'Features:  %s,\n' , featurelabels{ whrep } );

if whrep==1
    X  = features.stds( whok , : );
end

if whrep==2
    X  = features.means( whok , : );
    
    % Standardize?
    X = X -  repmat( mean( X , 1 ) , NY , 1 );
    X = X ./ repmat( std( X , [] , 1 ) , NY , 1 );
end

if whrep==3
    X = features.vars_scale(whok,:);
end

if whrep==4
    X = features.stds(whok,:).^2;
end

if whrep==5
    X = features.mssd(whok,:);
end

if whrep==6
    X = sqrt(features.mssd(whok,:));
end

if whrep==7
    X  = features.corrs( whok , : );
end

if whrep==8
    X = features.covs(whok,:);
end

if whrep==9
    X = [features.covs(whok,:) , features.stds(whok,:).^2];
end

if whrep==10
    X = features.covs_scale(whok,:);
end

X = double(X);
% number of features
NF = size( X , 2 );

%% Parallel loop
Y_pred_cell = cell( nsets , 1 );
YP_pred_cell = cell( nsets , 1 );
models = cell( nsets , 1 );

acc = zeros( nsets, NCO );
nonzeros = cell(nsets, NCO); 
for j=1:nsets
    fprintf('Fold=%d\n\tcost num=', j); 
    for c = 1:NCO
        
        if mod(c, 10) == 0
           fprintf('%d,', c);  
        end
        cost = costs(c); % make sure
        % costlabel = strrep( num2str(cost), '.', '_');
        costlabel = costlabels{c};
        whtrain = whtrain_sets{ j };
        whtest  = whtest_sets{ j };
        
        %% Train model using LIBLINEAR's logistic regression model
        Y_train = Y( whtrain );
        X_train = sparse( X( whtrain , : ));
        
        Y_test = Y( whtest );
        X_test = sparse( X( whtest , : ));
        
        % compute validation accuracy using the training data from
        % this fold
        options = sprintf( '-q -s %d -c %3.20f' , classifiertype , cost );
        model = train( Y_train , X_train, options );
        models{ j } = model;
        
        % Test model using code from LIBLINEAR
        warning off;
        test_options = '-b 0 -q'; % -b 1 outputs probability estimates
        [predict_label_train, accuracy_train, dec_values_train] = predict( Y_train, X_train, model , test_options );
        [predict_label_test, accuracy_test, dec_values_test] = predict( Y_test, X_test, model , test_options );
        
        % accuracy is [accuracy, MSE, squared correlation coefficient] 
        warning on;
        
        acc(j,c) = accuracy_test(1); 
        nonzeros{j,c} = sum(model.w~=0,2); 
        
        
    end
end

% Save info 
filenm = [ rd filesep 'feature_selection' filesep sprintf( 'pt%d_whs%d_c%d_whr%d_valid' ,...
    pt , whs , classifiertype , whrep ) ];

if dosave
    save( filenm , 'acc' , 'nonzeros', 'pt' , 'whs' , 'classifiertype' , 'whrep' ,...
        'featurelabels');
end

end

