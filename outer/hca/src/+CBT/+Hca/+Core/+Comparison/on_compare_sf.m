function [ rezMaxM] = on_compare_sf(barGen,theoryStruct,comparisonMethod,w,numPixelsAroundBestTheoryMask)
    % on_compare_sf 
    % Compares experiments with a stretch factor to single theory
    %     Args:
    %         sets: settings structure
    %         barcodeGen: barcode structure
    %         theoryStruct: theory structure
    % 
    %     Returns:
    %         rezMaxM: rezult structure

    import mp.mp_profile_stomp_nan_dna;
    switch comparisonMethod
        case 'mass_pcc'
            % choose k just higher than the length of small sequence for
            % best precision. (larger k though could increase speed)
            comparisonFun = @(x,y,z,w,u) unmasked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))),theoryStruct.isLinearTF,numPixelsAroundBestTheoryMask,1);
        case 'mpnan'
            % todo: this function itself does not include possibility of
            % being circular, so inputed X1, X2 should include extra
            % values..
            comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_profile_stomp_nan_dna(X1',X2', bitX1', bitX2', w,2^(4+nextpow2(length(X1))),numPixelsAroundBestTheoryMask);
        otherwise
            error('undefined method');
    end
            
    theorBar = theoryStruct.rawBarcode;
    if ~isempty(theoryStruct.rawBitmask)
        theorBit = theoryStruct.rawBitmask;
    else
        theorBit = ones(1,length(theorBar));
    end
        
%     theorBar = filter_barcode(theorBar,sets);
    
    rezMaxM = cell(1,5);
    rezMaxM{1} = zeros(length(barGen),length(barGen{1}.rescaled)); % this saves double
    rezMaxM{2} = int16(zeros(length(barGen),length(barGen{1}.rescaled))); % int
    rezMaxM{3} = int8(zeros(length(barGen),length(barGen{1}.rescaled)));
    rezMaxM{4} = int16(zeros(length(barGen),length(barGen{1}.rescaled)));
    rezMaxM{5} = int16(zeros(length(barGen),length(barGen{1}.rescaled)));
    rezMaxM{6}.info = {'1)PCC, 2) Pos, 3) Or, 4) SecondPos 5) Len '};

    ltheory = numel(theorBit);
%     cc = zeros(length(barGen{1}.rescaled),length(barGen));
    % for all the barcodes run
    % parfor
    for idx=1:length(barGen)
        for idy=1:length(barGen{idx}.rescaled) % loop over stretch factors
            if (isequal(comparisonMethod,'mpnan') && ((sum(barGen{idx}.rescaled{idy}.rawBitmask) < w)) || ltheory < w)% in case barcode stretch outside, or theory smaller
                rezMaxM{1}(idx,idy) = nan;
                rezMaxM{2}(idx,idy) = 1;
                rezMaxM{3}(idx,idy) = 1;
                rezMaxM{4}(idx,idy) = 1;
                rezMaxM{5}(idx,idy) = 1;
            else
                [rezMaxM{1}(idx,idy),rezMaxM{2}(idx,idy),rezMaxM{3}(idx,idy),rezMaxM{4}(idx,idy),rezMaxM{5}(idx,idy),~] =...
                    comparisonFun(barGen{idx}.rescaled{idy}.rawBarcode, theorBar, barGen{idx}.rescaled{idy}.rawBitmask,theorBit,w);
            end
%             cc(idy,idx) = rezMaxM(idx,idy).maxcoef(1);

%             [rezMaxM{idx}{idy}.maxcoef,rezMaxM{idx}{idy}.pos,rezMaxM{idx}{idy}.or,rezMaxM{idx}{idy}.secondPos,rezMaxM{idx}{idy}.lengthMatch,~] =...
%                 comparisonFun(barGen{idx}.rescaled{idy}.rawBarcode, theorBar, barGen{idx}.rescaled{idy}.rawBitmask,theorBit,w);
%             cc(idy,idx) = rezMaxM{idx}{idy}.maxcoef(1);
        end
    end

%     rezMax{1} = rezMaxM; % make single cell
%     rezMax{2} = cc;
end