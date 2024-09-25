function [chi2Score] = chi2_test(histAll, pdf,LU,nQuantiles)
    % chi2_test - goodness of fit test for distribution

    %   Args:
    %       histAll - histogram of values
    %       pdf - probability density function
    %       LU - bin edges, typically [1 length(histAll))]

    %   Returns:
    %       chi2Score - chi^2 for the comparison

    %   Example:
    %       import Zeromodel.chi2_test
    %       [chi2Score] = chi2_test(histAll, pdf,LU,nQuantiles)

    pdf = pdf/sum(pdf);
    cdf = cumsum(pdf);


    % two alternatives: even num of quantiles, or min bin size
    if nargin < 4
        minBinSize = 50;
        binEdges = find(histAll(1:end-1) >= minBinSize);
    else
        % find quantiles
        binEdges = zeros(1,nQuantiles+1);
        for i=1:nQuantiles-1
            % find the value where pvalue=1-cdf > pValThresh
            binEdges(i+1) = find(cdf > i/nQuantiles,1,'first');
        end
        %
        binEdges(1) = LU(1);
        binEdges(end) = LU(2);
    end
   
    
    numCounts = sum(histAll); % number of counts in histogram
    % bin counts
    binCountsFit = zeros(1,length(binEdges)-1);
    histAllGroup = zeros(1,length(binEdges)-1);
    for i=1:length(binEdges)-1
          binCountsFit(i) = numCounts*sum(pdf(binEdges(i):binEdges(i+1)-1));
          histAllGroup(i) = sum(histAll(binEdges(i):binEdges(i+1)-1));
    end

    
%         figure,plot((histAllGroup - binCountsFit).^2./binCountsFit)

    chi2Score = sum ( (histAllGroup - binCountsFit).^2./binCountsFit);
end

