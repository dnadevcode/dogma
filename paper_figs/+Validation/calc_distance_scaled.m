function     [posShift,sinkSource] = calc_distance_scaled(oS,sortedIdsAll,resStruct, lenThry, N)
    %   calc_distance_scaled
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

   import Core.synth_to_table;
   import Core.os_to_table;
    import Core.update_sf_barset;

    % for synthetic, there is two options: compare the position on the theory,
    % or compare the overlap positions
    posShift = zeros(1,N);
    sinkSource = zeros(N,2);
    for idxPair=1:N
        [curSink,curSource] = ind2sub(size(oS),sortedIdsAll(idxPair));
        sinkSource(idxPair,:) = [curSink,curSource];

        % create table for the barcode island:

        % synth to table
        [tableS] = synth_to_table(resStruct([curSource,curSink]));
        % os to table
        [tableR] = os_to_table(oS(curSink,curSource));

        bbS = tableS(1,4);  % best bar stretch (bbS)
        bbO = (tableS(1,3)~=tableR(1,3))+1;


        [tableRscaled] = update_sf_barset(tableR, bbS/tableR(1,4), bbO);

        posDifExp   = tableRscaled(1,1)-tableRscaled(2,1);

        posDifTheory  = tableS(1,1) - tableS(2,1); 

%         pos1 = resStruct{curSink}.pos(1);
%         pos2 = resStruct{curSource}.pos(1);
%         if resStruct{curSource}.or(1)==1
%             posDifTheory = pos1-pos2;
%         else
%             posDifTheory = -(pos1-pos2+resStruct{curSink}.lengthMatch-resStruct{curSource}.lengthMatch);
%         end
%         posDifExp = oS(curSink,curSource).pB - oS(curSink,curSource).pA;
%     
    
%         posShift(idxPair) = posDifTheory-posDifExp;% min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);

        posShift(idxPair) = min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);
    end

end

