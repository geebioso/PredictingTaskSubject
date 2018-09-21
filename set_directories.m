function [wd, rd] = set_directories(isHPC)

% working directory (code) 
if isHPC
    wd = '/data/users/ggaut/BV/'; 
else
   wd = fullfile( getenv('HOME'), 'Dropbox',  'anatomical', 'pred_task_subj'); 
end

% results directory (also main diarectory in this case) 
if isHPC 
    rd = '/pub/ggaut/BV/'; 
else
   rd = fullfile( getenv('HOME'), 'Dropbox', 'anatomical' );  
end