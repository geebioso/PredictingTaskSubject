function [] = plot_confusion_matrices( isHPC, whreps, whs, truncate, classifiertype, ...
    fontsize, fontcolor, usepatchimg, include_zeros, zero_char, facealpha)

% function to compute confusion matrices for each model 
% and output latex formatted tables 

% INPUT: 
%   bool isHPC: are we running on the UCI HPC? 
%   vector(int) whreps: which representations are we running? see below
%   int whs: which dataset? 
%   bool truncate: truncate the time series to shortest task? 
%   int classifiertype: same as liblinear documentation; 0 = L2, 6 = L1
%   int fontsize = 15; 
%   str fontcolor = 'k'; 
%   bool usepatchimg: otherise use imagesc
%   bool include_zeros: in latex file? 
%   char zero_char: character to replace zeros in latex file 
%   num facealpha: opacity of face [0,1] 

% whs==-1: old dataset (old OSU preprocessing pipeline) 
% whs== 0: baseline ( HCP preprocessing pipeline ) 
% whs== 1: dataset 1 ( HCP preprocessing pipeline ) 
% whs== 2: dataset 2 ( HCP preprocessing pipeline ) 
% whs== 3: dataset 3 ( HCP preprocessing pipeline )

% whrep==1: BVSD
% whrep==2: M
% whrep==3: BVVS
% whrep==4: BVV
% whrep==5: MSSD
% whrep==6: SQRT MSSD
% whrep==7: FCP
% whrep==8: FCC
% whrep==9: FCCV
% whrep==10: FCCS

pt = 1; % confusion matrices will be for within-session task classification 
NW = length(whreps); 


[wd, rd] = set_directories(isHPC); % working directory and results directory

load( fullfile(rd, 'ICAresults','tasklist') );

%% Load Results 
% get filenames 
for i = 1:length(whreps)
    whrep = whreps(i); 
    filenm = [sprintf( 'pt%d_whs%d_c%d_whr%d_truncate%d_valid' ,...
        pt , whs , classifiertype , whrep, truncate) ];
    models{i} = filenm; 
end
models = models'; 

% compute confusion and write to latex table 
C = cell(length(models), 1); 
NS = cell(length(models), 1);
accuracies = cell(length(models), 1); 
for i = 1:length(models)
   model = models{i};
   load(fullfile(rd, 'discrweights',model)); 
   
   has_pred = Y_pred ~= 0; 
   
   C{i} = confusionmat( YTRUE(has_pred), Y_pred(has_pred) ); 
   NS{i} = histc(YTRUE,unique(YTRUE));
   accuracies{i} = diag(C{i})./NS{i}; 
   
end

featurelabels = { ...
    'BVSD' , ...
    'M' , ...
    'BVVS',...
    'BV',...
    'MSSD', ...
    'SQRT_MSSD', ...
    'FCP' , ...
    'FCC', ...
    'FCCV', ...
    'FCCVS'
    };

labels = featurelabels(whreps); 
rep_string = strjoin( featurelabels, '_'); 
filename = fullfile( rd, 'confusion', sprintf( 'pt%d_whs%d_c%d_whr%s_latex.txt' , pt , whs , classifiertype , rep_string )); 
confusion2latex( C, accuracies, filename, tasks, include_zeros, zero_char,labels);
 
% Cdiff = {C{2} - C{1}}; 
% accdiff = {accuracies{2} - accuracies{1}}; 
% filename = fullfile( rd, 'confusion', sprintf( 'pt%d_whs%d_c%d_w_diff_whr%d_latex.txt' , pt , whs , classifiertype , rep_str)); 
% confusion2latex( Cdiff, accdiff, filename, tasks, include_zeros, zero_char, featurelabels(whreps));

%% Plots 

h0 = figure(); clf;
h0.Position = [61 285 2*1157 2*420];
h0.Position = [61 285 1157 420];
labelbounds = (1:9);

if usepatchimg
    labelbounds = labelbounds- 0.5;
    xdelta = 0.66;
    ydelta = 0.5;
else
    xdelta = 0.16;
    xdelta = 0.44;
    ydelta = 0.02;
end

h = cell(NW,1);
c = cell(NW,1); 
s = cell(NW,1); 
for i = 1:NW
    s{i} = subplot(1,NW,i);
%     scale = 0.9;
%     pos = get(s{i}, 'Position');
%     pos(2) = pos(2)+scale*pos(4);
%     pos(4) = (1-scale)*pos(4);
%     set(s{i}, 'Position', pos)
    
    if usepatchimg
        h{i} = patchimg(C{i}', 'none'); hold on; 
    else
        h{i} = imagesc(C{i}'); hold on; 
    end
    c{i} = colorbar(); caxis([0,19])
    alpha(h{i}, facealpha); 
    
    %title(sprintf('%s', whrep_labels{i})); 
    title(sprintf('%s', labels{i})); 
    
    s{i}.XTickLabels = tasks;
    s{i}.XTick =  labelbounds;
    s{i}.XTickLabelRotation = 35;
    s{i}.TickLength = [0 0];
    
    s{i}.YTickLabels = tasks;
    s{i}.YTick =  labelbounds;
    s{i}.FontSize = fontsize; 
     
    for k = 1:size(C{i},1)
        for j = 1:size(C{i},2)
            string = num2str(C{i}(k,j));
            xdelta2 = (3 - length(string))*0.06;
            text(k - xdelta + xdelta2, j - ydelta, string, 'Fontsize', 17, 'Color', fontcolor); 
        end
    end
    
 
end

% for i = 1:NW
%     % Adjust Colorbar Alpha
%       alphaVal = 0.5;
%     drawnow; 
%     cdata = c{i}.Face.Texture.CData;
%     cdata(end,:) = uint8(alphaVal * cdata(end,:));
%     c{i}.Face.Texture.ColorType = 'truecoloralpha';
%     c{i}.Face.Texture.CData = cdata;
%     drawnow
%     c{i}.Face.ColorBinding = 'discrete'; 
% end

filename = fullfile(rd, 'images', sprintf('confusion_matrices_pt%d_whs%d_c%d_whr%d_whr%d', ...
    pt , whs , classifiertype , whreps(1), whreps(2) ));
if usepatchimg
    print(filename, '-depsc');
else
    % saveas(figure(4), filename, 'png');
    
    print(filename, '-depsc', '-r0'); 
end