function [fig1] = plot_max_coef_vs_filtered( fig1,comparisonStruct,comparisonStructFiltered, numBar, sets, markers )
    % plot_max_coef - pltos three maximum coefficients
    
    maxcoef = cell2mat(cellfun(@(x) x.maxcoef,comparisonStruct,'UniformOutput',false)');
    
    p = plot(maxcoef(1:numBar,:),1:numBar,'ob');
    p(1).Marker = markers(1);
    p(2).Marker = markers(2);
    p(3).Marker = markers(3);
    
    if  sets.genConsensus == 1
        plot([0.1:0.1:1], 0.5+repmat(numBar,10,1));
        p3 = plot(maxcoef(numBar+1,:),numBar+1,'ob');
        p3(1).Marker = markers(1);
        p3(2).Marker = markers(2);
        p3(3).Marker = markers(3);
       % p4(1).Marker = 'none';
    end

    ylabel('Barcode nr.','Interpreter','latex')
    xlabel('Maximum match score','Interpreter','latex')
    xlim([0.5 1])
    ylim([0,numBar+2])
    legend({'$\hat C$','$C_2$','$C_3$','Consensus line'},'Location','sw','Interpreter','latex')
end

