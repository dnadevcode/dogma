function [theory] = run_literature_theory(ts,sets)
    % run literature theory


    import CBT.Hca.Core.Theory.cb_theory;



    % free concentrations of yoyo and netropsin
    cN = sets.theoryGen.concN;
    cY = sets.theoryGen.concY;
    yoyoBindingConstant = sets.model.yoyoBindingConstant;
    values = sets.model.netropsinBindingConstant;



    theory = cb_theory(ts, cN,  cY, yoyoBindingConstant, values, 0);

end

