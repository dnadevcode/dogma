function [hcaSessionStruct] = compare_figure_of_merit(hcaSessionStruct, sets )
    % the main function to compare fragments of human chromosome vs. theory
    
    % using figure-of-merit function.
    % see http://www.pnas.org/content/110/13/4893.full.pdf?with-ds=yes
    % article
    
    % input hcaSessionStruct, sets 
    
    % output comparisonStructure? (move it outside hcaSessionStruct?)
        
    disp('Starting comparing exp to theory...')
    tic
    
    rawBarcodes = hcaSessionStruct.rawBarcodes;
    rawBitmasks = hcaSessionStruct.rawBitmasks;
    if sets.barcodeConsensusSettings.aborted==0
        rawBarcodes = [rawBarcodes; hcaSessionStruct.consensusStruct.barcode];
        rawBitmasks = [rawBitmasks hcaSessionStruct.consensusStruct.bitmask];
    end
    
    
    hcaSessionStruct.comparedStructure = cell(1,length(hcaSessionStruct.theoryGen.theoryBarcodes));
    %%%%%%%%%%%%%%%%%   
    % unfiltered comparison
    for barNr = 1:length(hcaSessionStruct.theoryGen.theoryBarcodes)
        theorBar = hcaSessionStruct.theoryGen.theoryBarcodes{barNr};
        theorBit = hcaSessionStruct.theoryGen.bitmask{barNr};
        % what should be the best filter size?
        filterSize = sets.filterSettings.filterSize;

        % create the rezult structure
        comparisonStructure = cell(length(rawBarcodes),1);

        % stretch factors
        stretchFactors = sets.barcodeConsensusSettings.stretchFactors;

        % make available to all workers.
        theorBar = parallel.pool.Constant(theorBar);
        theorBit = parallel.pool.Constant(theorBit);
         
  
        % we use the standard xcorralign function (check if can change this to
        % improve speed)
        import SignalRegistration.XcorrAlign.get_no_crop_lin_circ_xcorrs;
        import CBT.Hca.Core.Comparison.get_cc_fft;
        % for all the barcodes run
        for i=1:length(rawBarcodes)
            % xcorrMax stores the  maximum coefficients
            xcorrMax = zeros(1,length(stretchFactors));

            % rezMaz stores the results for one barcode
            rezMax = cell(1,length(stretchFactors));


            % barTested barcode to be tested
            barTested = rawBarcodes{i};

            % length of this barcode
            lenBarTested = length(barTested);

            % barBitmask - bitmask of this barcode
            barBitmask = rawBitmasks{i};

            % make them available to workers
            barTested = parallel.pool.Constant(barTested);
            barBitmask = parallel.pool.Constant(barBitmask);
            tic
            % run the parloop for the stretch factors
            parfor j=1:length(stretchFactors)
                % here interpolate both barcode and bitmask 
                barC = interp1(barTested.Value, linspace(1,lenBarTested,lenBarTested*stretchFactors(j)));

                % if we need to filter the barcode,
%                 if str==1
%                     barC = imgaussfilt(barC, filterSize);
%                 end

                barB = barBitmask.Value(round(linspace(1,lenBarTested,lenBarTested*stretchFactors(j))));
                
                
                % compute the scores.
                % here separate which barcode is the longer and which is the
                % shorter (todo: check how this works when we compare fragment
                % to fragment theory)
                if length(barC) > length(theorBar.Value)
                   % [xcorrs,~,~] =  CBT.Hca.Core.Comparison.get_cc_fft(theorBar.Value,barC, theorBit.Value,barB);
                    [xcorrs, ~, ~] = SignalRegistration.XcorrAlign.get_no_crop_lin_circ_xcorrs(theorBar.Value,barC, theorBit.Value,barB);
                else
                    tempb = [theorBar.Value theorBar.Value(1:length(barC(barB)))];                 
                    barC = zscore(barC(barB));
                    fliped = fliplr(barC);
                    for jj=1:length(theorBar.Value)
                       xcorrs(1,jj) = 1/(2*length(barC))*sum((barC-zscore(tempb(jj:jj+length(barC)-1))).^2);
                       xcorrs(2,jj) = 1/(2*length(barC))*sum((fliped-zscore(tempb(jj:jj+length(barC)-1))).^2);
                    end                  %  [xcorrs, ~, ~] = CBT.Hca.Core.Comparison.get_cc_fft(barC, theorBar.Value, barB,theorBit.Value);
                   % [xcorrs, ~, ~] = get_no_crop_lin_circ_xcorrs(barC, theorBar.Value, barB,theorBit.Value);
                end 

                % now find the maximum score for this stretching parameter
                xcorrMax(j) = min(xcorrs(:));
                [rezMax{j}.maxcoef,rezMax{j}.pos,rezMax{j}.or] = CBT.Hca.UI.Helper.get_best_parameters(xcorrs,length(barC) );
            end       
            toc
            % find which stretching parameter had the best score
            [~,b] = max(xcorrMax);

            % select the results for this best stretching parameter and output
            % them.
            comparisonStructure{i} = rezMax{b};
            comparisonStructure{i}.bestBarStretch = stretchFactors(b);
%             if str==1
%                     barC = imgaussfilt(interp1(rawBarcodes{i}, linspace(1,length(rawBarcodes{i}),length(rawBarcodes{i})*comparisonStructure{i}.bestBarStretch)), filterSize);
%                     comparisonStructure{i}.bestStretchedBar = barC;
%             else
               comparisonStructure{i}.bestStretchedBar = interp1(rawBarcodes{i}, linspace(1,length(rawBarcodes{i}),length(rawBarcodes{i})*comparisonStructure{i}.bestBarStretch));
%             end
            comparisonStructure{i}.bestStretchedBitmask = rawBitmasks{i}(round(linspace(1,length(rawBitmasks{i}),length(rawBarcodes{i})*comparisonStructure{i}.bestBarStretch)));
        end 
        hcaSessionStruct.comparedStructure{barNr} = comparisonStructure;
   
    end
    
    timePassed = toc;
    disp(strcat(['Experiments were compared to theory in ' num2str(timePassed) ' seconds']));
end

