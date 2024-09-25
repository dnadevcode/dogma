function     posShift = calc_pos_dif(oS,sortedIdsAll,resStruct, lenThry, N)
    %   calc_pos_dif
    %
    %   Args:
    %       oS - overlap structure,sortedIdsAll - sorted indexes from best
    %       to lowest score
    %       resStruct - ground truth positions
    %       N - number of scores to calculate
    %
    %   Returns:
    %       posShift - positional shift between detected and ground truth
    %       position

    % for synthetic, there is two options: compare the position on the theory,
    % or compare the overlap positions
    posShift = zeros(1,N);
    for idxPair=1:N
        
        [curSink,curSource] = ind2sub(size(oS),sortedIdsAll(idxPair));
        pos1 = resStruct{curSink}.pos(1);
        pos2 = resStruct{curSource}.pos(1);
        if resStruct{curSource}.or(1)==1
            posDifTheory = pos1-pos2;
        else
            posDifTheory = -(pos1-pos2+resStruct{curSink}.lengthMatch-resStruct{curSource}.lengthMatch);
        end
        posDifExp = oS(curSink,curSource).pB - oS(curSink,curSource).pA;
    
    
%         posShift(idxPair) = posDifTheory-posDifExp;% min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);

        posShift(idxPair) = min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);
    end

end

