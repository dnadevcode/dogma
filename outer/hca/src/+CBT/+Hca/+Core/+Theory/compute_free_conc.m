function [ sets ] = compute_free_conc( sets )
    % compute_free_conc
    % Computing free concentrations
    %     Args:
    %         sets (struct): Input settings
    % 
    %     Returns:
    %         sets: Output settings
    %

    import CBT.Hca.Core.Theory.cb_theory;

    if sets.theoryGen.computeFreeConcentrations
        
        %% load lambda sequence, preferably already in the settings file
        lambdaSequence = fastaread(strcat([sets.lambda.fold sets.lambda.name]));

        tic
        disp('Computing free concentrations');
        ntIntSeq = nt2int( lambdaSequence.Sequence, 'ACGTOnly',1);

        probsBinding = @(x) cb_theory(ntIntSeq, x(2),x(1),sets.model.yoyoBindingConstant,sets.model.netropsinBindingConstant, 1000, 2);

        x0 = [sets.theoryGen.concY sets.theoryGen.concN];

        fun = @(x) x0-x-mean(probsBinding(x))*sets.theoryGen.concDNA*0.25;

        % we minimise the sum square
        fun2 = @(x) sum(fun(x).^2);

        % using fminsearch
        [xNew] = fminsearch(fun2,x0);

        sets.theoryGen.concY = xNew(1);
        sets.theoryGen.concN = xNew(2);
        toc
        disp('Finished computing free concentrations');

    end


end

