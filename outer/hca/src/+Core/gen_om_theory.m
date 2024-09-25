function [theory, header, bitmask] = gen_om_theory(name,sets)
    % gen_om_theory / rewritten compute_theory_barcode function

    % if free conc. needs to be generated
%     import CBT.Hca.Core.Theory.compute_free_conc;
%     sets = compute_free_conc(sets);
    
    
    import CBT.Hca.Core.Theory.create_memory_struct;
    [chr1, header] = create_memory_struct(name);
    
      
    if length(chr1.Data) < sets.theoryGen.m
        disp(strcat(['Theory ' name ' skipped because length < theoryGen.m =' num2str(sets.theoryGen.m) ]));
        theory = [];  bitmask = [];
        delete(chr1.Filename);
        clear chr1;
        return;
    end
        
%     sets.theoryGen.k=2^17

%     sets.theoryGen.k = 2^17;%length(chr1.Data)+sets.theoryGen.m;
    % set-up-theory-generation
    if sets.addgc
        import Core.set_up_theory_for_calculation_gc;
    [~, Y, cutPointsL, cutPointsR, extraL, extraR, ...
        ~, ~, k, m, circular, ~, ~, numPx,~, ~,bpPlaceL,bpPlaceR ] = ...
        set_up_theory_for_calculation_gc(sets, chr1);
    else
        import Core.set_up_theory_for_calculation;
        [~, Y, cutPointsL, cutPointsR, extraL, extraR, ...
        ~, ~, k, m, circular, ~, ~, numPx,~, ~,bpPlaceL,bpPlaceR ] = ...
        set_up_theory_for_calculation(sets, chr1);
    end
%     
%     round(bpPlaceL{2}(1)) % todo: make sure that there's exactly 1 bp
% %     difference between last and first element of each batch
%     round(bpPlaceR{1}(end))

%     round(bpPlaceL{3}(1))
%     round(bpPlaceR{2}(end))
% 
%     round(bpPlaceL{4}(1))
%     round(bpPlaceR{3}(end))

    % generate theory & bitmask/ using cutPointsL and cutPointsR
    import Core.generate_theory_cg;
    [theory,bitmask] = generate_theory_cg(Y, chr1, cutPointsL, cutPointsR, extraL, extraR, ...
     k, m, circular, numPx, sets);
    
    
    % clear chr1
    delete(chr1.Filename);
    clear chr1;
    
    
% end
% %% A C G T
% chr2.Data=ones(length(chr1.Data),1);
% chr2.Data((1:5)+5)=3;
% chr2.Data(40000:40004)=3;
% % chr2.Data(end-4:end)=3;
% 
% 
% % chr2.Data(1:end/2)=3;
%  [~, Y, cutPointsL, cutPointsR, extraL, extraR, ...
%         ~, ~, k, m, circular, ~, ~, numPx,~, ~,bpPlaceL,bpPlaceR ] = ...
%         set_up_theory_for_calculation(sets, chr2);
%     sets.theoryGen.method = 'simple';
% circular=0;
%     [theory,bitmask] = generate_theory(Y, chr2, cutPointsL, cutPointsR, extraL, extraR, ...
%      k, m, circular, numPx, sets);
% %  figure,plot(fftshift(theory))
% 
%  circular=1;
%     [theory2,bitmask] = generate_theory(Y, chr2, cutPointsL, cutPointsR, extraL, extraR, ...
%      k, m, circular, numPx, sets);
%  
%   figure,plot(fftshift(theory))
%  hold on
%  plot(fftshift(theory2))
% % plot(fftshift(theory))