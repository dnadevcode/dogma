function [barMat,bars,orBars,reFac,pval] = single_bargroup(bar1,pthresh,PCC_OVERLAP,PCC_MP,len1,len2,lenOverlap,...
        stridx,mpI1,LOCS1,pksUnique1,pksUniquePos1,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,pl)

%   Creates single bargroup for the specified input.

% todo: for this keep only the edges between the considered barcodes. Can
% create overlap graph.

%   Args:
%     PCC_OVERLAP,PCC_MP,len1,len2,lenOverlap

%   Returns:
%       
% index of ROOT barcode
%%
% bar1 = 20;
% figure,plot(PCC_MP(bar1,:))
% figure,plot(PCC_OVERLAP(bar1,:))

pccVals = PCC_OVERLAP(bar1,:);
idxToCheck = find(pccVals>0.3);

pval = nan(1,length(pccVals));
import Zeromodel.beta_ev_pdf;
for j=idxToCheck
        a = 0.13*lenOverlap(bar1,j);
        n = 2*(len1(bar1,j)+len2(bar1,j)-2*lenOverlap(bar1,j));
        pval(j) = 1-Zeromodel.beta_ev_cdf( PCC_OVERLAP(bar1,j), a, 1, n, 1);
end

    numGoodMatch = pval < pthresh;

    peaksToTry = find(numGoodMatch);
    
    % bargroup..
    import Core.plot_bargroup;
    if pl==1
        [barMat,bars,orBars,reFac] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,pl);
    else
        [barMat,bars,orBars,reFac] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);
    end        
    % instead create matrix like in gen reference based assembly
%     figure;hold on;
%     for i=1:length(barMat)
%         plot(  barMat{i}{1},zscore(  barMat{i}{2} )+3*i,'red');
%     end

end

