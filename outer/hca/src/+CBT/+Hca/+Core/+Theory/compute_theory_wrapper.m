function [x] = compute_theory_wrapper(ts, sets)
    % theory wrapper.  Can compute multiple theory methods based on the
    % input in the settings files
    
    % Args:
    %    ts : time series structures
    %    sets : settings structure
        
    import CBT.Hca.Core.Theory.run_simple_theory;
    import CBT.Hca.Core.Theory.run_literature_theory;
    
%     if sum(ts > 4) > 0 % if there are nan's, don't compute the theory
%         x = nan(length(ts),1);
%     else
    switch sets.theoryGen.method
        case "v1"
            import CBT.Hca.Core.Theory.run_simple_theory_v1;
            x = run_simple_theory_v1(ts);
        case "test1"
            x = run_simple_theory(ts);
        case "simple"
            x = run_simple_theory(ts);
        case "literature"
             x = run_literature_theory(ts',sets);
        case "custom" % same as literature here
             x = run_literature_theory(ts',sets);
        case "GC"
            x = gc_rate(ts',4);
        case "gcweighted"
            import CBT.Hca.Core.Theory.gcweighted;

            x = gcweighted(ts',4,sets.theoryGen.atPreference);
%              x = run_literature_theory(ts',sets);
        case "enzyme"
             x = zeros(length(ts),1);
             dots = strfind(ts',sets.model.pattern);
             x(dots) = 1;
%         case "TCGA"
%              x = zeros(length(ts),1);
%              dots = strfind(ts',sets.model.pattern);
%              x(dots) = 1;
        otherwise
            error('Incorrect method for theory generation selected' )
    end
%     end


end

