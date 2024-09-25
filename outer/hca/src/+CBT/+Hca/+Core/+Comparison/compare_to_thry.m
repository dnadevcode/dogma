function [ maxVals ] = compare_to_thry( theoryStruct, barnr, lenBarcodes, sets )

    if nargin < 4
        import CBT.Hca.Import.set_default_settings_code;
        sets = set_default_settings_code();
    end
    
    if nargin < 3
        lenBarcodes = 300;
    end
    
    if nargin < 2
        barnr = 24;
    end

    % number of barcodes. Might be some pixels left at the end of the
    % theory
    numBars = floor((length(theoryStruct.theoryGen.theoryBarcodes{barnr}))/lenBarcodes);
    
    maxVals = zeros(1,numBars);
    
    for ii=1:numBars
        disp(strcat(['Computing comparison for bar nr ' num2str(ii) ' out of ' num2str(numBars) ]));
        % which pixels to take
        xx = (1+(ii-1)*lenBarcodes):(ii*lenBarcodes);
        % we cut out a barcode
        bar1 = theoryStruct.theoryGen.theoryBarcodes{barnr}(xx);
        
        % create new hcastructure to save this barcode
        hcaSessionStructNew = [];
        hcaSessionStructNew.rawBarcodes{1} = bar1;
        % TODO: do we add some randomization for this barcode? can add
        % as an option randomization to get something like 0.9 cc value at
        % the best place
        
        hcaSessionStructNew.rawBitmasks{1} = ones(1,length(xx));
        hcaSessionStructNew.theoryGen = theoryStruct.theoryGen;
        % remove the exact pixels of this barcode from theory
        hcaSessionStructNew.theoryGen.theoryBarcodes{barnr}(xx) = [];
        % no consensus

        % compare to theory
        import CBT.Hca.UI.compare_theory_to_exp;
        hcaSessionStruct = compare_theory_to_exp(hcaSessionStructNew, sets);

        % combine results
        import CBT.Hca.UI.combine_chromosome_results;
        hcaSessionStruct = combine_chromosome_results(hcaSessionStruct,sets);
        
        % put the results here
        maxVals(ii) = hcaSessionStruct.comparisonStructure{1}.maxcoef(1);
    end



end

