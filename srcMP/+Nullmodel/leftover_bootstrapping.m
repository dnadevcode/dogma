function [pccScoreB] = leftover_bootstrapping(wBoostrapping,pairedPartialBar, partialLength, NN, minL,s)
%{
    Boot-strapp leftover scores to be of length wBoostrapping. Only do this if
    there are at least minL pixels to bootstrap from

%}

if nargin < 6
    s = RandStream('mlfg6331_64'); 
end


pccScoreB = cell(1,length(pairedPartialBar));
for j=1:length(pairedPartialBar)
    if partialLength(j) > minL
        B = pairedPartialBar{j};
    
        pccScore = zeros(1,NN);
        for i=1:NN
            y = datasample(s,1:length(B),wBoostrapping,'Replace',true);
        
            b1 = B(1,y);
            b2 = B(2,y);
            pccScore(i) = zscore(b1(:)',1)*zscore(b2(:),1)/length(b1(:)); % todo: speed up by pre-calculating the residuals and only scale by mean
        
        end
        pccScoreB{j} = pccScore;
    else
        pccScoreB{j} = nan;
    end
end
end

