function [wd, rd] = set_directories(isHPC)

% working directory 
if isHPC
    wd = '/data/users/ggaut/BV/'; 
else
   wd = fullfile( getenv('HOME'), 'Dropbox',  'anatomical'); 
end

% results directory 
if isHPC 
    rd = '/pub/ggaut/BV/'; 
else
   rd = wd;  
end