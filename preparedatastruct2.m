function [ D , NS , NT , NR , N ] = preparedatastruct2( tc , tasks, subjs, session_idx)

%% Prepare a new structure for data
NS = size( tc , 1 ); % Number of subjects
NT = length( tasks ); % Number of tasks
NR = size( tc{1,1} , 2 ); % Number of ROIs

[unique_subjs, ~, subj_idx] = unique(subjs); 

% Determine total number of observations
N = 0;
toffset = zeros( NT , 1 );
for j=1:NT
    M = size( tc{1,j} , 1 );
    N = N + M * NS;
    if (j>1)
       toffset( j ) = toffset( j-1 ) + M;
    end
end
fprintf( 'Total number of cases: %d\n' , N );
D = table( zeros(N,1) , zeros(N,1) , zeros(N,1) , zeros(N,1), zeros(N,1), zeros(N,NR),...
           'VariableNames' , { 'SubjectIndex' , 'TaskIndex' , 'TrialIndex', 'TrialOverallIndex' , 'SessionIndex' , 'ROI' } );

verbosity = 0;       
       
%% Put ROI Data in Table
ii = 0;
for j=1:NT
    tasknow = tasks{ j };
    fprintf( 'Preparing the data for task: %s\n' , tasknow );
    for i=1:NS
        xnow = tc{i,j};
        M = size( tc{1,j} , 1 );
           
        if verbosity == 2
            figure( 900 ); clf;
            histogram( xnow(:));
            title( sprintf( 'Task:%s Subject:%d\n' , tasknow , i ));           
            pause;
        end
        
        D.SubjectIndex( ii+1:ii+M ) = subj_idx(i);
        D.TaskIndex( ii+1:ii+M )    = j;
        D.TrialIndex( ii+1:ii+M )   = 1:M;
        D.TrialOverallIndex( ii+1:ii+M )   = (1:M) + toffset( j );
        D.SessionIndex( ii+1:ii+M ) = session_idx(i); % only session 1 data for now
        D.ROI( ii+1:ii+M ,: )       = xnow;
        ii = ii + M;
    end
end

fprintf('Complete\n'); 
