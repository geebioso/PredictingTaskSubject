function [ datasetnames] = load_dataset_names()

datasetnames = cell(1,2);

datasetnames(1,:) = {'tc_ica_wm_csf_out_180429_conv.mat', 0};
datasetnames(2,:) = {'tc_gm_out_180430_conv.mat', 1};
datasetnames(3,:) = {'tc_ica_gm_out_180521_conv.mat', 2}; 
datasetnames(4,:) = {'tc_ica_gm_task_out_180523_conv.mat', 3};
datasetnames(5,:) = {'tc_ica_wm_csf_out_180429_conv.mat', 26};
datasetnames(6,:) = {'tc_gm_out_180430_conv.mat', 22}; 
datasetnames(7,:) = {'tc_ica_gm_out_180521_conv.mat', 23};
datasetnames(8,:) = {'tc_ica_gm_task_out_180523_conv.mat', 27};
datasetnames(8,:) = {'tc_csf_wm_motion_out_globalMask5000_203subj_conv.mat', -1};




