function [x] = run_simple_theory(ts)
    % run ismple theory
    %   Args:
    %       ts
    %   Returns:
    %       x
    %
    
    
    % cummulative sum of GC's. 
    numWsCumSum = cumsum((ts == 1)  | (ts == 4) );

    %  Find all ligands without any A's or T's, i.e. i and i+4 should
    %  have the same value. We start at the left of the first possible
    %  ligand
    x = [numWsCumSum(5:end) == numWsCumSum(1:end-4); 0; 0; 0; 0];
  
end

