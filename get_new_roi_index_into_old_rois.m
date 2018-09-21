function [new_rois_old_idx] = get_new_roi_index_into_old_rois(rd) 
    addpath(fullfile(rd, '..','FMRI','hTFA_Toolbox','code','dependencies','NIfTI_Toolbox'));
    nii269 = load_nii( fullfile(rd, '..','FMRI', 'restingstatedata', 'ROI_269_3mm.nii.gz' ));
    nii299 = load_nii( fullfile(rd, 'ICAresults','ROI_299_3mm.nii.gz' ));
    load(fullfile(rd, '..','FMRI', 'restingstatedata', 'nVoxels_ROI269.mat')); 

    % compare images 
    old_rois = 1:299;
    new_rois = zeros(size(old_rois)); 
    for i =1:length(old_rois)
        old_roi = old_rois(i);

        ii = nii299.img == old_roi;

        possible_rois = unique(nii269.img(ii)); 
        if length(possible_rois) > 1
            new_rois(i) = possible_rois(possible_rois ~= 0); 
        else
            new_rois(i) = possible_rois;
        end
    end

    [old_rois', new_rois']; 

    new_rois_old_idx = old_rois(new_rois ~=0); 

end