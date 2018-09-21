

[wd, rd] = set_directories(isHPC); % working directory and results directory

addpath(fullfile(rd, '..','FMRI','BrainVisualization','Network Classification'));
addpath(fullfile(rd, '..','FMRI','BrainVisualization','BrainNetViewer'));

filenm = fullfile( rd, '..', 'FMRI', 'restingstatedata', 'ROI_269_3mm.nii.gz' ); 

nii = load_nii(filenm); 

view_nii(nii)

