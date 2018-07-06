
whsets = [27];
whreps = 7:10;
whpts = 1:2;
isHPC = 0;
Nworkers = min( 4, feature('numCores')); 

%% Start parallel pool

% set number of workers to minimum of number of sims, or number of cores

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end

if ~(poolsize == Nworkers)
    delete(gcp('nocreate'))
    
    parpool(Nworkers);
    fprintf('Parallel using %d cores\n', Nworkers);
end


%% Run classifiers

[p1 p2 p3] = meshgrid(whsets, whreps, whpts);
allpairs = [p1(:), p2(:) p3(:)];

parfor p = 1:size(allpairs,1)
    whs = allpairs(p,1);
    whrep = allpairs(p,2);
    pt = allpairs(p,3); 
    runclassifier_fully_parallelized(isHPC, whs, whrep, pt);
end



