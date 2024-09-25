function [outConsensus] = barcodes_to_scaled(barcodeGen, comparisonStruct, theoryStruct)
    % convert original barcodes to correctly scaled along the theory
    % TODO: do this from the output table

%     if nargin < 4 || isempty(islandTableRepresentation)
        % This scaling is based on placement along the theory


%     if ~iscell(comparisonStruct)
%         import Plot.islandsPosStruct;
%         posStruct = islandsPosStruct(barcodeIslandsData,barIslands);
%     end


        % convert the reference coefs to p-vals, etc.
        orientation = cell2mat(cellfun(@(x) x.or,comparisonStruct,'UniformOutput',0)');
        pos = cell2mat(cellfun(@(x) x.pos,comparisonStruct,'UniformOutput',0)');

%         maxcoef = cell2mat(cellfun(@(x) x.maxcoef,comparisonStruct,'UniformOutput',0)');
        strF = cell2mat(cellfun(@(x) x.bestBarStretch,comparisonStruct,'UniformOutput',0)');
        lengthsbar = cellfun(@(x) length(x.rawBarcode),barcodeGen);

        posCur = pos;

        if nargin >=3 &&~isempty(theoryStruct)
            % load theory file. just single one
            theorBar = theoryStruct{comparisonStruct{1}.idx}.rawBarcode;
          
            %% plot multi theories with some nan's in between
            outConsensus = nan(length(barcodeGen)+1,length(theorBar));
            outConsensus(1,:) = theorBar;
            lenT = length(theorBar);
        else % need to estimate the length
            posCur = posCur-min(pos)+1;
            %% plot multi theories with some nan's in between
            lenT = max(posCur+round(lengthsbar'.*strF));

            outConsensus = nan(length(barcodeGen)+1,lenT);

        end

        for k=1:length(barcodeGen)
            ii = k;
            expBar = barcodeGen{ii}.rawBarcode;
            expBit = barcodeGen{ii}.rawBitmask;
    
            expLen = length(expBar);
    
            % interpolate to the length which gave best CC value
            expBar = interp1(expBar, linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch));
            expBit = expBit(round(linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch)));
            expBar(~expBit)= nan;
            
            if orientation(ii,1)~=1
                expBar = fliplr(expBar);
                expBit = fliplr(expBit);
            end
        
            if posCur(ii) < 1
                posCur(ii) = posCur(ii)+lenT;
            end
            posEnd = posCur(ii)+length(expBar)-1;
            numEltsOver = posEnd - lenT;
            numFirst = lenT - posCur(ii)+1;
            if posEnd > lenT
    %             seqToplot =  (expBar-nanmean(expBar))/nanstd(expBar,1)+5;
                outConsensus(k+1,[posCur(ii):lenT 1:numEltsOver]) = [expBar(1:numFirst) expBar(numFirst+1:end)] ;
            else
                outConsensus(k+1,posCur(ii):posCur(ii)+length(expBar)-1) = expBar;
            end
                
        end

%     else % independent on the theory
%     end

end

