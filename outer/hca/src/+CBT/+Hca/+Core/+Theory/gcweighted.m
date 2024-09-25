function [x] = gcweighted(seq,blockSize, atPreference)
    if nargin < 3
        atPreference = 1;
    end
    numNT = numel(seq);
    %     blockSize = 500; % how many bp to make the bloc size?
    numBlocks = floor(numNT/blockSize);

    ratio = zeros(numBlocks+1,1);

    A = nt2int('A'); C = nt2int('C'); G = nt2int('G'); T = nt2int('T');

    x = zeros(numNT,1);
    for count = 1:numBlocks
        % calculate the indices for the block
        start = 1 + blockSize*(count-1);
        stop = blockSize*count;
        % extract the block
        block = seq(start:stop);
        % find the GC and AT content
        gc = (sum(block == G | block == C));
        at = (sum(block == A | block == T));
        % calculate the ratio of GC to the total known nucleotides
        ratio(count) = gc/(gc+at*atPreference);
        x(start:stop) =  ratio(count);
    end

    block = seq(stop+1:end);
    gc = (sum(block == G | block == C));
    at = (sum(block == A | block == T));
    ratio(end) = gc/(gc+at*atPreference);
    x(stop+1:end)= ratio(end);
%     for i=1:length(ratio)-1
%         score = repmat(
%     end
    
end

