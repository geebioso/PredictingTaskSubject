function [ lobes, lobenames, primaryareas, primarylobes, regionareas, regionlobes, regionareas2, regionlobes2, lobenames2, primaryareas2, primarylobes2] = getLobeSets( names )
% script to return the set of primary areas taken from our list of regions
% that correspond to each of the main lobes of the brain
% INPUT: 
%   names[N x 1]: list of region names as provided in ROI file 

% OUTPUT: 
%   lobes: structure of with lobes as fields and primary areas as values
%   lobenames: names of lobes 
%   primaryares: cell of unique primary areas 
%   primarylobes: cell of lobes corresponding to each primary area
%   regionareas: primary areas corresponding to each region name
%   regionlobes: lobes corresponding to each region name
%   regionares2: brodmann areas+ corresponding to each region name
%   regionlobes2: lobes corresponding to each region name ( decided by brodmann area)
%   lobenames2: unique lobes in regionlobes2
%   primaryareas2: unique brodmann areas 

% main lobes: frontal, parietal, occipital, temporal, limbic, cerebellum, insular

% add more information for regions that just have brodmann number
i = strcmp(names, '_00'); 
names{i} = 'brainstem_00'; 

i = strcmp(names, '_17'); 
names{i} = 'Occipital_17'; 
%% 

lobenames = {'Frontal' 'Parietal' 'Occipital' 'Temporal' 'Limbic' 'Cerebellum' 'Insular' 'Misc'};

lobes = struct();

lobes.Frontal = {'Frontal' 'Olfactory' 'Paracentral' 'Rectus' 'Supp', 'Precentral'};

lobes.Parietal = {'Angular' 'Parietal' 'Postcentral' 'Precuneus' 'SupraMarginal', 'Rolandic' };

lobes.Occipital = {'Calcarine' 'Cuneus' 'Lingual' 'Occipital' 'optic', '_17'};

lobes.Temporal = {'Fusiform' 'Heschl' 'Temporal', 'closeToAC' };

lobes.Limbic = {'Amygdala' 'Cingulum' 'Hippocampus' 'ParaHippocampal' 'Putamen' 'Thalamus', 'Caudate'};

lobes.Cerebellum = {'BrainStem' 'Cerebelum' 'Vermis' 'brainstem', 'PonTop', '_00'};

lobes.Insular = {'Insula'};

% lobes.Misc = {char.empty(1,0)};


[ regionareas2, regionlobes2, lobenames2, primaryareas2, primarylobes2 ] = getBrodmann(names);

primaryareas = names;
regionlobes = names;
startidx = regexp(names, '_');
for i = 1:length(names)
    if ~isempty(startidx{i})
        primaryareas{i} = names{i}(1:(startidx{i}(1)-1));
    else
        primaryareas{i} = names{i};
    end
    
    regionlobes{i} = regionlobe(primaryareas{i}, lobes, lobenames);
end

regionareas = primaryareas;
primaryareas = unique(primaryareas); % sorted

primarylobes = primaryareas;
for i = 1:length(primaryareas)
    primarylobes{i} = regionlobe(primaryareas{i}, lobes, lobenames);
end

end

function [ lobe ] = regionlobe( name, lobes, lobenames)
F = length(lobenames);
found = 0;
for i = 1:(F-1)
    if sum(ismember(lobes.(lobenames{i}), name))
        found = 1;
        lobe = lobenames{i};
        break
    end
end

if found == 0
    lobe = lobenames{end};
end
end

function [ regionareas2, regionlobes2, lobenames2, primaryareas2, primarylobes2 ] = getBrodmann(name)
[ brodmannlabels, brodmannlobes ] = BrodmannLabels( );
regionareas2 = name;
regionlobes2 = name;
for i = 1:length(name)
    bn = str2double(name{i}((end-1):end));
    if bn >0
        regionareas2{i} = brodmannlabels{bn};
        regionlobes2{i} = brodmannlobes{bn};
    else
        if strcmp(name{i}(1:min(9, end)), 'Cerebelum');
            if regexp(name{i},'Crus1') > 0;
                regionareas2{i} = 'Crus1';
            elseif regexp(name{i}, 'Crus2') > 0
                regionareas2{i} = 'Crus2';
            elseif regexp(name{i}, 'Mid') > 0
                regionareas2{i} = 'Mid Cerebellum';
            else
                st = regexp(name{i}, '_');
                idx = (st(1) + 1):( st(2) -1);
                regionareas2{i} = [name{i}(idx) ' Cerebellum'];
            end
            regionlobes2{i} = 'Cerebellum';
        elseif strcmp(name{i}(1:min(6, end)), 'Vermis')
            regionareas2{i} = 'Vermis';
            regionlobes2{i} = 'Cerebellum';
        elseif strcmp(name{i}(1:min(9, end)), 'BrainStem')
            regionareas2{i} = 'NA';
            regionlobes2{i} = 'Brain Stem';
        elseif strcmp(name{i}(1:min(6, end)), 'PonTop')
            regionareas2{i} = 'Pons Top';
            regionlobes2{i} = 'Brain Stem';
        elseif strcmp(name{i}(1:min(8, end)), 'Thalamus')
            regionareas2{i} = 'Thalamus';
            regionlobes2{i} = 'Limbic';
        elseif strcmp(name{i}(1:min(9, end)), 'Precuneus')
            regionareas2{i} = 'Precuneus';
            regionlobes2{i} = 'Parietal';
        elseif strcmp(name{i}(1:min(7, end)), 'Putamen')
            regionareas2{i} = 'Putamen';
            regionlobes2{i} = 'Basal Ganglia';
        elseif strcmp(name{i}(1:min(7, end)), 'Caudate')
            regionareas2{i} = 'Caudate Nucleus';
            regionlobes2{i} = 'Basal Ganglia';
        elseif strcmp(name{i}(1:min(8, end)), 'Cingulum')
            if strcmp(name{i}(10:14), 'Mid_L') > 0
                regionareas2{i} = 'Left Mid Cingulum';
            elseif strcmp(name{i}(10:14), 'Mid_R') > 0
                regionareas2{i} = 'Right Mid Cingulum';
            end
            regionlobes2{i} = 'Limbic';
        elseif strcmp(name{i}(1:min(9, end)), 'brainstem')
            regionareas2{i} = 'Brain Stem';
            regionlobes2{i} = 'Cerebellum';
        elseif regexp(name{i}, 'Rectus') > 0
            regionareas2{i} = 'Rectus';
            regionlobes2{i} = 'Frontal';
        elseif regexp(name{i}, 'ParaHippocampal') > 0
            regionareas2{i} = 'Parahippocampal';
            regionlobes2{i} = 'Limbic';
        elseif regexp(name{i}, 'Insula') > 0
            regionareas2{i} = 'Insula';
            regionlobes2{i} = 'Insular';
        elseif regexp(name{i}, 'Temporal') > 0
            regionareas2{i} = 'Superior Temproal Gyrus';
            regionlobes2{i} = 'Temporal';
         elseif regexp(name{i}, 'optic_chiasm') > 0
            regionareas2{i} = 'NA';
            regionlobes2{i} = 'Optic Chiasm';  
        elseif regexp(name{i}, 'closeToAC') >0
            regionareas2{i} = 'NA';
            regionlobes2{i} = 'Anterior Commissure';  
        else
            regionareas2{i} = 'NA';
            regionlobes2{i} = 'NA';
        end
        regionareas2{i} = [regionareas2{i} '(B00)']; 
    end
end

lobenames2 = unique(regionlobes2); 
[primaryareas2, b, c] = unique(regionareas2); 
primarylobes2 = regionlobes2(b);
end



