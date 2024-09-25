function [] = plot_theory_vs_experiment( bar,barBit,b1,b1Bit,orientation,pos,maxcoef,name,theoryName,titleT)
    % this plots theory vs experiments, when exp is longer than theory

%     import SignalRegistration.XcorrAlign.get_no_crop_lin_circ_xcorrs;
%     [xcorrs, ~, ~] = get_no_crop_lin_circ_xcorrs(bar,b1,barBit,b1Bit);
% 
%     [xcorrs2,~,~] =  CBT.Hca.Core.Comparison.get_cc_fft(bar,b1,barBit,ones(1,length(b1)));
%     [xcorrs2,~,~] =  CBT.Hca.Core.Comparison.get_cc_fft(bar,b1,barBit,b1Bit);

   % [m,p,o] = CBT.Hca.UI.Helper.get_best_parameters(xcorrs,length(bar) );
  
    if orientation == 2
        pos = pos + length(b1Bit)-length(bar);
        bar = fliplr(bar);
        barBit = fliplr(barBit);
    end
    
% 
% 
	fitPositions = pos:pos+length(bar)-1;
% 
%     if sum(fitPositions<=0) > 1
%         indx = find(fitPositions<=0);
%         indx = length(bar)-fliplr(indx);
%         fitPositions(find(fitPositions<=0)) = indx;
% %                fitPositions(find(fitPositions<=0)) = indx;
%    % fitPositions = fitPositions(1:min(end,length(bar)));
%     end
% a = find(fitPositions>length(bar));

  %  fitPositions(find(fitPositions>length(bar))) = fitPositions(find(fitPositions>length(bar))) -length(bar);

            
    a = find(fitPositions == length(b1)+1);
    if ~isempty(a)
        fitPositions(a:end) = 1:length(fitPositions(a:end));
    else
        a =  find(fitPositions == 0);
        fitPositions(1:a) = (length(b1)-a+1):length(b1);
        a = a+1;
    end

    
    barFit = b1(fitPositions);
    barBit = b1Bit(fitPositions);
%     [xcorrs, ~, ~] = get_no_crop_lin_circ_xcorrs(bar,barFit,barBit,barBit);

    m1 = mean(bar(logical(barBit)));
    s1 = std(bar(logical(barBit)));

    m2 = mean(barFit(logical(barBit)));
    s2 = std(barFit(logical(barBit)));
    
    

    hold on

    plot(((b1-m2)/s2)*s1+m1)
    if ~isempty(a)
        plot([fitPositions(1:a-1) nan fitPositions(a:end)],[bar(1:a-1) nan bar(a:end) ])
    else
        plot([fitPositions(1:end)],[bar(1:end)])
    end
%     hold on
%     plot(fitPositions,barFit)
    %xlim([pos(ii,1) pos(ii,1)+length(b1)-1 ])
    xlabel('pixel nr.')
    ylabel('Rescaled to theoretical intesity')
%     if ii <= len1
%         name = strcat([name]);
%     else
%         name = 'consensus';
%     end
    name = strrep(name,'_',' ');
    name = strrep(name,'kymograph.tif','');
    theoryName =  strrep(theoryName,'_',' ');
   % title(strcat(['Molecule ' name ]));
    title(titleT,'Interpreter','latex')
    legend({strcat(['$\hat C=$' num2str(maxcoef,'%0.2f')]),theoryName(1:min(20,end))},'Interpreter','latex')

end

