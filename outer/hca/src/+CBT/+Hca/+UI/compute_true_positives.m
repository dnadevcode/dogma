function [hcaSessionStruct] = compute_true_positives(hcaSessionStruct, sets )
% compute_true_positives

% ? Add to future features list (or fix directly if easy)? 1. Also add statistics for the filtered barcodes? 2. For both filtered and unfiltered add i) how many barcodes that pass consensus hit correct position ii) how many barcodes that passes BOTH consensus and p-value threshold find the correct position. 3. What is the actual p-value threshold? Can I set it myself = ask in the same window as for position in pixel + error? Would be good if possible :)
    lengths =  cellfun(@length, hcaSessionStruct.theoryGen.theoryBarcodes);
    sets.correctPlace = sets.correctPlace +sum(lengths(1:sets.correctChromosome-1));
    
    
    positions = cell2mat(cellfun(@(x) x.pos(1),hcaSessionStruct.comparisonStructure,'UniformOutput',0));
    left = sets.correctPlace-sets.allowedError;
    right = sets.correctPlace+sets.allowedError;
    allPos =  (left<positions).*(positions<right);
    
       
    % shorthand for if there is consensus
    withC = sets.barcodeConsensusSettings.aborted;
    
    % correctly matching barcodes (without consensus)
    disp('Correctly matching barcodes')
    strLineMol = strcat([num2str(sum(allPos(1:end-1+withC))) ' out of ' num2str(length(allPos)-1+withC)]);
    disp(strLineMol)
        
    % correctly matching consensus
    if withC == 0        
    	disp('Correctly matching consensus')
        strLine = strcat([num2str(allPos(end)) ' out of ' num2str(1)]);
        disp(strLine)

        disp('Correctly matching barcodes that pass the consensus threshold')
        strLineMol = strcat([num2str(sum(allPos(sets.barcodeConsensusSettings.barcodesInConsensus))) ' out of ' num2str(length(sets.barcodeConsensusSettings.barcodesInConsensus))]);
        disp(strLineMol)
    end
    
    
     pValComp = isfield(hcaSessionStruct,'pValueResults');
     if pValComp
         % barcodes passing p-val thresh
        passingPThresh=hcaSessionStruct.pValueResults.pValueMatrix(1:end-1+withC) <sets.pvaluethresh;
        disp('Method 1. Barcodes passing p-value thresh')
        strLine = strcat([num2str(sum(passingPThresh)) ' out of ' num2str(length(passingPThresh))]);
        disp(strLine);
        % barcodes that pass the p-value thresh and match at a correct
        % place
        disp('Method 1. Correcly matching barcodes passing p-value thresh')
        strLine = strcat([num2str(sum(allPos(1:end-1+withC).*passingPThresh')) ' out of ' num2str(sum(passingPThresh))]);
        disp(strLine)
             
         if withC == 0;
             disp('Method 1. Correctly matching barcodes that pass both the p-value thresh and the consensus threshold')
             strLine = strcat([num2str(sum(allPos(sets.barcodeConsensusSettings.barcodesInConsensus).*passingPThresh(sets.barcodeConsensusSettings.barcodesInConsensus)')) ' out of ' num2str(sum(passingPThresh(sets.barcodeConsensusSettings.barcodesInConsensus)))]);

            % strLineMol = strcat([num2str(sum(allPos(sets.barcodeConsensusSettings.barcodesInConsensus))) ' out of ' num2str(length(sets.barcodeConsensusSettings.barcodesInConsensus))]);
             disp(strLine)
         end
     end
     
  
     
     pValComp = isfield(hcaSessionStruct,'pValueResultsOneBarcode');
     if pValComp
         passingPThresh = hcaSessionStruct.pValueResultsOneBarcode.pValueMatrix(1:end-1+withC) <sets.pvaluethresh;
         disp('Method 2. Barcodes passing p-value thresh')
         strLine = strcat([num2str(sum(passingPThresh)) ' out of ' num2str(length(passingPThresh))]);
         disp(strLine);
         disp('Method 2. Correcly matching barcodes passing p-value thresh')
         strLine = strcat([num2str(sum(allPos(1:end-1+withC).*passingPThresh')) ' out of ' num2str(sum(passingPThresh))]);
         disp(strLine)
             
         if withC == 0
             disp('Method 2.  Correctly matching barcodes that pass both the p-value thresh and the consensus threshold')
             strLine = strcat([num2str(sum(allPos(sets.barcodeConsensusSettings.barcodesInConsensus).*passingPThresh(sets.barcodeConsensusSettings.barcodesInConsensus)')) ' out of ' num2str(sum(passingPThresh(sets.barcodeConsensusSettings.barcodesInConsensus)))]);
             disp(strLine)
         end
     end
     
     
    if sets.filterSettings.filter==1
    	positions = cell2mat(cellfun(@(x) x.pos(1),hcaSessionStruct.comparisonStructureFiltered,'UniformOutput',0));
        allPos =  (left<positions).*(positions<right);
        disp('Correctly matching filtered barcodes')
        strLineMol = strcat([num2str(sum(allPos(1:end-1+withC))) ' out of ' num2str(length(allPos)-1)]);
        disp(strLineMol)
        if withC == 0;
            disp('Correctly matching filtered consensus')
            strLine = strcat([num2str(allPos(end)) ' out of ' num2str(1)]);
            disp(strLine)
         
            disp('Correctly matching filtered barcodes that pass the consensus threshold')
            strLineMol = strcat([num2str(sum(allPos(sets.filterSettings.barcodesInConsensus))) ' out of ' num2str(length(sets.filterSettings.barcodesInConsensus))]);
            disp(strLineMol)
        end
    end

  %   pValComp = isfield(hcaSessionStruct,'pValueResults');
%      if pValComp
%         passingPThresh = hcaSessionStruct.pValueResults.pValueMatrixFiltered(1:end-1+withC) < sets.pvaluethresh;
%         disp('Method 1. Filtered barcodes passing p-value thresh')
%         strLine = strcat([num2str(sum(passingPThresh)) ' out of ' num2str(length(passingPThresh))]);
%         disp(strLine);
%         disp('Method 1. Correcly matching filtered barcodes passing p-value thresh')
%         strLine = strcat([num2str(sum(allPos(1:end-1+withC).*passingPThresh')) ' out of ' num2str(sum(passingPThresh))]);
%         disp(strLine)
% 
%         if withC == 0;
%             disp('Method 1.  Correctly matching filtered barcodes that pass both the p-value thresh and the consensus threshold')
%             strLine = strcat([num2str(sum(allPos(sets.filterSettings.barcodesInConsensus).*passingPThresh(sets.filterSettings.barcodesInConsensus)')) ' out of ' num2str(sum(passingPThresh(sets.filterSettings.barcodesInConsensus)))]);
%             disp(strLine)
%         end
%      end
    
%      pValComp = isfield(hcaSessionStruct,'pValueResultsOneBarcode');
%      if pValComp
%         passingPThresh = hcaSessionStruct.pValueResultsOneBarcode.pValueMatrixFiltered(1:end-1+withC) < sets.pvaluethresh;
%         disp('Method 2. Filtered barcodes passing p-value thresh')
%         strLine = strcat([num2str(sum(passingPThresh)) ' out of ' num2str(length(passingPThresh))]);
%         disp(strLine);
%         disp('Method 2. Correcly matching filtered barcodes passing p-value thresh')
%         strLine = strcat([num2str(sum(allPos(1:end-1+withC).*passingPThresh')) ' out of ' num2str(sum(passingPThresh))]);
%         disp(strLine)
%              
%          if withC== 0;
%              disp('Method 2.  Correctly matching filtered barcodes that pass both the p-value thresh and the consensus threshold')
%              strLine = strcat([num2str(sum(allPos(sets.filterSettings.barcodesInConsensus).*passingPThresh(sets.filterSettings.barcodesInConsensus)')) ' out of ' num2str(sum(passingPThresh(sets.filterSettings.barcodesInConsensus)))]);
%              disp(strLine)
%          end
%      end

end

