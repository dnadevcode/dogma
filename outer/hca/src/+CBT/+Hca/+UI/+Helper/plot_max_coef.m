function [fig1,maxcoef] = plot_max_coef( fig1,comparisonStruct, numBar, sets, markers )
    % plot_max_coef - pltos three maximum coefficients
    
    maxcoef = cell2mat(cellfun(@(x) x.maxcoef(1:min(3,end)),comparisonStruct,'UniformOutput',false)');
    
    p = plot(fig1,maxcoef,1:size(maxcoef,1),'ob');
    p(1).Marker = markers(1);
    try
        p(2).Marker = markers(2);
    end
    try
        p(3).Marker = markers(3);
    end
        
    
    if  sets.genConsensus == 1
        plot(fig1,0.1:0.1:1, 0.5+repmat(numBar,10,1));
%         p3 = plot(maxcoef(numBar+1,:),numBar+1,'ob');
%         p3(1).Marker = markers(1);
%         p3(2).Marker = markers(2);
%         p3(3).Marker = markers(3);
       % p4(1).Marker = 'none';
    end

    ylabel(fig1,'Barcode nr.','Interpreter','latex')
    xlabel(fig1,'Maximum match score','Interpreter','latex')
	try    
        xlim(fig1,[min(maxcoef(:)) max(maxcoef(:))]);
    catch
         xlim(fig1,[0.5 1]);
    end
    ylim(fig1,[0,size(maxcoef,1)+2])
    if size(maxcoef,2) == 1
        legend(fig1,{'$\hat C$','Theories separator'},'Location','southoutside','Interpreter','latex')
    else
        if size(maxcoef,2)==3
            legend(fig1,{'$\hat C$','$C_2$','$C_3$'},'Location','southoutside','Interpreter','latex')
        else
            legend(fig1,{'$\hat C$','$C_2$','$C_3$','Consensus line'},'Location','southoutside','Interpreter','latex')

        end
    end
    title(fig1,'Max score','Interpreter','latex');
end

