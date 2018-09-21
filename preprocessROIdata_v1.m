%% Read and prepare ROI data

% whs: which dataset to use
%   -1 = non-hcp processing (old pipeline) 
%   0 = baseline
%   1 = dataset 1 
%   2 = dataset 2
%   3 = dataset 3 

whs = 1; 
data_directory = fullfile(getenv('HOME'),'Dropbox','FMRI','restingstatedata', 'OSF_data_dump'); 

% Original dataset (original processing pipeline) 
if (whs==-1) 
    removeoutliers = 0;
    mainstudy   = fullfile( data_directory, 'tc_csf_wm_motion_out_globalMask5000_203subj.mat');
    secondstudy = '';
    savefile    = ['ICAresults' filesep 'tc_csf_wm_motion_out_globalMask5000_203subj_conv.mat'];
    load( fullfile(data_directory, 'tc_csf_wm_motion_out_globalMask5000_203subj'), 'subj_ind_in_250' );  
end 

% Paper datasets (variations on HCP pipeline) 
if (whs==0) 
    removeoutliers = 0;
    mainstudy   = fullfile( data_directory, 'tc_ica_wm_csf_out_180429.mat');
    secondstudy = '';
    savefile    = ['ICAresults' filesep 'tc_ica_wm_csf_out_180429_conv.mat'];
    load( fullfile(data_directory,'tc_csf_wm_motion_out_globalMask5000_203subj'), 'subj_ind_in_250' );  
end

if (whs==1) 
    removeoutliers = 0;
    mainstudy   =  fullfile( data_directory, 'tc_gm_out_180430.mat');
    secondstudy = '';
    savefile    = ['ICAresults' filesep 'tc_gm_out_180430_conv.mat'];
    load( fullfile(data_directory,'tc_csf_wm_motion_out_globalMask5000_203subj'), 'subj_ind_in_250' );  
end

if (whs==2) 
    removeoutliers = 0;
    mainstudy   =  fullfile( data_directory, 'tc_ica_gm_out_180521.mat');
    secondstudy = '';
    savefile    = ['ICAresults' filesep 'tc_ica_gm_out_180521_conv.mat'];
    load( fullfile(data_directory,'tc_csf_wm_motion_out_globalMask5000_203subj'), 'subj_ind_in_250' );  
end

if (whs==3) 
    removeoutliers = 0;
    mainstudy   =  fullfile( data_directory, 'tc_ica_gm_task_out_180523.mat');
    secondstudy = '';
    savefile    = ['ICAresults' filesep 'tc_ica_gm_task_out_180523_conv.mat'];
    load( fullfile(data_directory,'tc_csf_wm_motion_out_globalMask5000_203subj'), 'subj_ind_in_250' );  
end

% Load the task names
load( 'tasklist.mat');

% Load 250 subject names
%load( 'icaresults\subjlist250' ); % subjs;
%fullsubs = subjs;


%% Load data
load( mainstudy );
% [ D , NS1 , NT , NR , N1 ] = preparedatastruct2( tc , tasks );
[ D , NS1 , NT , NR , N1 ] = preparedatastruct2( tc , tasks, subjs, session_idx);

% % Need to change sessionindex and subject index
% wh = D.SubjectIndex > 174;
% D.SessionIndex( wh ) = 2;
% 
% % Create a new index
% D.SubjectIndex = subj_ind_in_250( D.SubjectIndex )';
% [ subs , ~ , D.SubjectIndex2 ] = unique( D.SubjectIndex );

%% Show histogram
figure( 1 ); clf;
subplot( 2,1,1 );
histogram( D.ROI(:) , 500 );
title( 'Original ' );

%% Remove outliers?
if removeoutliers == 1
    ROI = D.ROI;
    %premove = 0.05;
    %lw = prctile( ROI(:) , premove );
    
    ROI( ROI <= 500 ) = NaN;
    ROI( ROI >= 15000 ) = NaN;
    ROI2 = round( ROI , -2 );
    
    globalmean = mode( ROI2( : ));
    ROI( isnan( ROI )) = globalmean;
    
    
    subplot( 2,1,2 );
    histogram( ROI(:) , 500 );
    title( 'After removing outliers' );
    
    D.ROI = ROI;
end

%% Save
save( savefile , 'D' );
