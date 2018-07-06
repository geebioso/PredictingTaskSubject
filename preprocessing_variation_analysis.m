
classifiertype = 0; 
isHPC = 0; 
whsets = [ 1:4 ]; 

baseline = 'whset1'; 

[ AC , best_costs, tab ] = print_result_table_regularization_optimized( classifiertype, whsets, isHPC);


%%  Differences in Task/No Task 


fields = fieldnames( AC.(baseline)); 
F = length(fields); 
D = zeros(F, NPT); 

for f = 1:F 
    fieldnow = fields{f}; 
    
    D(f,:) = AC.(baseline).(fieldnow) - AC.whset27.(fieldnow); 
   
end 

%% Other Differences 

whsets = [ 3 ]; 
fields = fieldnames( AC.(baseline)); 
F = length(fields); 
maxD = zeros(F,NPT); 
for whs = whsets 
   
    whs_field = sprintf('whset%d', whs); 
    for f = 1:F
        
        fieldnow = fields{f};
        difference = AC.(baseline).(fieldnow) - AC.(whs_field).(fieldnow); 
        
        if sum( abs(difference))  > sum( abs(maxD(f,:)))
           maxD(f,:) = difference;  
        end
        
    end
    
end