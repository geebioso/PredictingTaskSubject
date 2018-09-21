function [] = plot_BV_rest_vs_task( isHPC, whs, truncate, fontsize) 

% This function produces a plot of rest BV versus task BV 
%% Options 
if whs==0
    bad_roi = 253; % ROI with really low variance 
else
    bad_roi = []; 
end

%% Load data 
[wd, rd] = set_directories(isHPC); % working directory and results directory
filenm = fullfile( rd, 'features', sprintf('whs%d_truncate%d.mat', whs, truncate));
load(filenm); 

roi_stds = features.stds; 
roi_stds(:, bad_roi) = []; 

%% Compute average rest and task BV 
restidx = combs(:,2) == 5;
taskidx = combs(:,2) ~= 5;

restvar = mean(roi_stds(restidx,:));
taskvar = mean(roi_stds(taskidx,:));

lim_min = min( [restvar(:), taskvar(:)]);
lim_max = max( [restvar(:), taskvar(:)]);
lims = [lim_min, lim_max];

figure(); clf;
scatter(restvar, taskvar, 'bo'); hold on;
plot(lims, lims, 'k');
xlabel( 'Average Resting State BV');
ylabel( 'Average Task BV');
set(gca, 'FontSize', fontsize)

filename = fullfile( rd, 'images', sprintf('avg_rest_std_dev_vs_ave_task_std_dev_whs%d_truncate%d', whs, truncate));
print(filename, '-depsc');


