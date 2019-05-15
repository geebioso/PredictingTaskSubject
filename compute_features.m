function [] = compute_features(whsets, isHPC, truncate)

[wd, rd] = set_directories(isHPC);

for whs = whsets
    
    fprintf( 'whset = %d\n', whs);
    
    % Original dataset (original processing pipeline) 
    if (whs==-1)
        roifile = [rd filesep 'ICAresults' filesep 'tc_csf_wm_motion_out_globalMask5000_203subj_conv.mat'];
    end
    
    % Paper Datasets (variations on HCP pipeline) 
    if (whs==0)
        roifile = [rd filesep 'ICAresults' filesep 'tc_ica_wm_csf_out_180429_conv.mat'];
    end
    
    if (whs==1)
        roifile = [rd filesep 'ICAresults' filesep 'tc_gm_out_180430_conv.mat'];
    end
    
    if (whs==2)
        roifile = [rd filesep 'ICAresults' filesep 'tc_ica_gm_out_180521_conv.mat'];
    end
    
    if (whs==3)
        roifile = [rd filesep 'ICAresults' filesep 'tc_ica_gm_task_out_180523_conv.mat'];
    end
    
    % Debugging Datasets 
    if (whs==7)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_csf_wm_motion_task_out_globalMask5000_203subj_conv'];
        icfile     = '';
    end
    
    if (whs==8)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_csf_wm_motion_out_globalMask5000_203subj_conv'];
        icfile     = '';
    end
    
    if (whs==8)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_csf_wm_motion_out_globalMask5000_203subj_conv'];
        icfile     = '';
    end
    
    if (whs==9)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_motion_task_out_ROI305_203subj_conv.mat'];
        icfile     = '';
    end
    
    if (whs==10)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_motion_wm_csf_task_out_ROI305_conv.mat'];
        icfile     = '';
    end
    
    if (whs==11)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_miniIca_motion_wm_csf_task_out_ROI305_conv.mat'];
        icfile     = '';
    end
    
    if (whs==12)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_motion_wm_csf_task_out_ROI305_noICA_conv.mat'];
        icfile     = '';
    end
    
    if (whs==13)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_miniIca_motion_wm_csf_task_out_ROI269_conv.mat'];
        icfile     = '';
    end
    
    if (whs==14)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_miniIca_motion_wm_csf_out_ROI305_conv.mat'];
        icfile     = '';
    end
    
    if (whs==15)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_miniIca_motion_wm_csf_task_out_ROI305_180409_conv.mat'];
        icfile     = '';
    end
    
    if (whs==16)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_compare_to_dataset8_180412_conv.mat'];
        icfile     = '';
    end
    
    if (whs==17)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_compare_to_dataset8_noICA_180416_conv.mat'];
        icfile     = '';
    end
    
    if (whs==18)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_dataset8_repeat_180417_conv.mat'];
        icfile     = '';
    end
    
    if (whs==19)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_hp200_180419_conv.mat'];
        icfile     = '';
    end
    
    if (whs==20)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_ica_wm_csf_task_out_180426_conv.mat'];
        icfile     = '';
    end
    
    if (whs==21)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_wm_csf_out_180427_conv.mat'];
        icfile     = '';
    end
    
    if (whs==22)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_ica_wm_csf_out_180429_conv.mat'];
        icfile     = '';
    end
    
    if (whs==23)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_gm_out_180430_conv.mat'];
        icfile     = '';
    end
    
    if (whs==24)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_wm_csf_task_out_180430_conv.mat'];
        icfile     = '';
    end
    
    
    if (whs==25)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_wm_csf_out_ROI305_1800501_conv.mat'];
        icfile     = '';
    end
    
    
    if (whs==26)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_ica_gm_out_180521_conv.mat'];
        icfile     = '';
    end
    
    if (whs==27)
        roifile    = [rd filesep 'ICAresults' filesep 'tc_ica_gm_task_out_180523_conv.mat'];
        icfile     = '';
    end
    
    
    %% Load Data (loads D)
    load(roifile);
    load( fullfile(rd, 'ICAresults','tasklist') );
    
    %% Constants
    NS = max( D.SubjectIndex2 ); % number of subjects
    NT = length( tasks ); % number of tasks
    NR = size( D.ROI , 2 ); % number of regions
    
    %% Get combinations
    [ combs , ~ , cindex ] = unique( [ D.SubjectIndex2 D.TaskIndex D.SessionIndex ], 'rows' );
    NC = size( combs , 1 );
    
    % Find subjects that are in session 1 and session 2
    nspersubj = accumarray( combs(:,1) , combs(:,3) , [ NS 1 ] , @mean );
    subsinbothsessions = find( nspersubj == 1.5 );
    
    %% Get Length of Shortest Task
    num_time_points = zeros(NT,1);
    for t = 1:NT
        task_now_idx = D.TaskIndex == t;
        num_time_points(t) = max(D.TrialIndex(task_now_idx));
    end
    scan_lengths = num_time_points*2 + 12;          % tasks have a TR of 2s 
    scan_lengths(5) = num_time_points(5)*2.5 + 12;  % resting state has a TR of 2.5s 
    
    min_scan_length = min(scan_lengths);
    
    %% Compute Features
    
    ROI = D.ROI;
    features = struct();    
    
    NRP = NR * (NR-1) / 2; % off diagonal entries of FC matrix
    roi_covs  = zeros( NC , NRP );
    roi_corrs = zeros( NC , NRP );
    roi_means = zeros( NC , NR  );
    roi_stds  = zeros( NC , NR  );
    mssd  = zeros( NC , NR  );
    mad  = zeros( NC , NR  );
    
    warning off;
    triidx = ~triu(true(NR));
    for j=1:NC
        if mod(j,100) == 0
            fprintf( '\tcombination %d of %d\n' , j , NC );
        end
        wh = find( cindex == j );
        ROI_NOW = ROI( wh , : );
        task_now = unique(D.TaskIndex(wh)); 
        
        % truncate based on minimum scan length and TR 
        if truncate
            if task_now == 5 % rest 
                TR = 2.5;
            else % all other tasks 
                TR = 2; 
            end
            time_point_cutoff = round( (min_scan_length - 12)/TR ); % 12 seconds burn-in 
            ROI_NOW = ROI_NOW( 1:time_point_cutoff, :);  
        end
        
        roi_means( j , : ) = mean( ROI_NOW , 1 );
        roi_stds( j , : )   = std( ROI_NOW , [] , 1 );
        mssd( j, : ) = compute_mssd(ROI_NOW);
        mad( j, : ) = compute_mad( ROI_NOW); 
        
        
        rhos  = 1 - pdist( ROI_NOW' , 'correlation' );
        rhos2 = cov(ROI_NOW);
        rhos2 = rhos2(triidx)';
        
        roi_corrs( j , : ) = rhos;
        roi_covs( j, : ) = rhos2;
        
    end
    warning on;
    
    %     Do this computation elsewhere in runclassifier_fully_parallelized to save space 
    %       roi_vars_scale = roi_stds.^2;
    %       roi_vars_scale = (roi_vars_scale - mean(roi_vars_scale(:)))/std(roi_vars_scale(:));
    %       roi_covs_scale = (roi_covs - mean(roi_covs(:)))/(std(roi_covs(:)));
    
    % Remove NaNs
    roi_corrs( isnan( roi_corrs )) = 0;
    roi_covs( isnan( roi_covs )) = 0;
    
    
    %% Save Features
    features.corrs = single(roi_corrs);
    features.covs = single(roi_covs);
    features.means = single(roi_means);
    features.stds = single(roi_stds);
    features.mssd = single(mssd);
    features.mad = single(mad); 
    
    filenm = fullfile( rd, 'features', sprintf( 'whs%d_truncate%d.mat', whs, truncate));
    fprintf('Saving %s\n', filenm); 
    save( filenm, 'features', 'combs', 'NC', 'NR', 'NT', 'subsinbothsessions');
    
end