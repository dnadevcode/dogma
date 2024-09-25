function [pxCutLeft,pxCutRight,px] = px_cut_pos(numWsCumSum, gcSF, pxSize)

% Find px positions, start at 1
px= 1;

pxCutLeft = [];
pxCutRight = [];
pxStart = 1;

numPrev = 0;
sumCur = numWsCumSum(1);
lenCur = 1;
for i = 2:length(numWsCumSum)
    sumCur = numWsCumSum(i)-numPrev;
    lenCur = lenCur+1;
    lenActual = sumCur+(lenCur-sumCur)*gcSF; % CG's are scaled by gcSF

    if (lenActual>=round(pxSize))
        pxCutLeft(px) = pxStart;
        pxCutRight(px) = i;
        pxStart = i+1;

        px = px+1;
        numPrev = numWsCumSum(i);
        lenCur = 0;
    end
end

end

