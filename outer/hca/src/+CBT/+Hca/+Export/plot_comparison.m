function [  ] = plot_comparison(ii,dd,len1,pos,orientation,b1, b1Bit,bar,barBit,hcaSessionStruct )
    % plot_comparison
                  
    if orientation(ii,1) == 2
        b1 = fliplr(b1);
        b1Bit = fliplr(b1Bit);
    end
    
    fitPositions = pos(ii,1):pos(ii,1)+length(b1)-1;
    
    %
    
    % temporary fix
    if length(bar) < length(fitPositions)
        bar = [bar zeros(1,length(fitPositions))];
        barBit = [barBit zeros(1,length(fitPositions))];
    end
%     
      has0 = 0;
    if sum(fitPositions<=0) > 1
        has0 = 1;
        indx = find(fitPositions<=0);
        indx = length(bar)-fliplr(indx);
        fitPositions(find(fitPositions<=0)) = indx;     
    end

    
   % barFit = bar(fitPositions);
%    if sum(fitPositions<0) > 0
%        fitPositions = fitPositions+length(bar);
%    end
%    
   barFit = bar(fitPositions);
   
    barBit = barBit(fitPositions);
    
    m1 = mean(barFit(logical(b1Bit)));
    s1= std(barFit(logical(b1Bit)));
    
    m2 = mean(b1(logical(b1Bit)));
    s2= std(b1(logical(b1Bit)));
    
    if has0 == 1;
        [a,b] = find(fitPositions == 1);
        fitPositions = [fitPositions(1:b-1) NaN fitPositions(b:end)];
        b1 = [b1(1:b-1) NaN b1(b:end)];
        barFit = [barFit(1:b-1) NaN barFit(b:end)];
    end

    plot(fitPositions,((b1-m2)/s2) *s1+m1)
    hold on
    plot(fitPositions,barFit)
    %xlim([pos(ii,1) pos(ii,1)+length(b1)-1 ])
    xlim([min(fitPositions) max(fitPositions)])

%     ax = gca;
%     ticks = fitPositions(1):30:fitPositions(end);
%     ticksx = floor(ticks);
%     ax.XTick = [ticks];
%     ax.XTickLabel = [ticksx/1000];
    xlabel('Position (px)','Interpreter','latex')
    ylabel('Rescaled to theoretical intesity','Interpreter','latex')
    if ii <= len1
        name = strcat([hcaSessionStruct.names{ii}]);
    else
        name = 'consensus';
    end
    name = strrep(name,'_',' ');

    name = strrep(name,'kymograph.tif','');

    title(strcat(['Unfiltered barcode ']),'Interpreter','latex');

    legend({strcat(['$\hat C_{\rm M' num2str(ii) '}=$' num2str(dd,'%0.2f')]),hcaSessionStruct.comparisonStructure{ii}.name},'Interpreter','latex')
    %[xcorrs,~,~] =  CBT.Hca.Core.Comparison.get_cc_fft(b1,barFit,b1Bit,barBit);

        

end

