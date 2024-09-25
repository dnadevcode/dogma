function [  ] = plot_comparison_exp_vs_exp(selectedIndices,comparisonStruct,theoryStruct,barcodeGen  )
    % plot_comparison_exp_vs_exp plot one barcode vs the other
    % we don't compare them anew, but just plot them on the same plot based
    % on their positions on the theory barcode. 

 import CBT.Hca.UI.Helper.load_theory_and_stretch_ex;

 if length(selectedIndices) >= 2
    % first barcode is used as a reference
    nr1 = selectedIndices(1);
    % load both barcode and corresponding best theory
    [ barT,~, bar1, bit1] = load_theory_and_stretch_ex(nr1, theoryStruct, comparisonStruct,barcodeGen );
    pos1 =  comparisonStruct{nr1}.pos(1);
    max1 = comparisonStruct{nr1}.maxcoef(1);
    or1 = comparisonStruct{nr1}.or(1);
    %flip barcodes in case not the same orientation
    if or1 == 2
        bar1 = fliplr(bar1);
        bit1 = fliplr(bit1);
    end
    % assign nan's where bitmask is 0
    bar1(~bit1) = nan;
         % plot barcode, here we zscore
	figure,
    plot((barT-nanmean(barT))./nanstd(barT),'black')
    hold on            
	plot(pos1:pos1+length(bar1)-1,(bar1-nanmean(bar1))./nanstd(bar1))
    legendInfo = {};
	legendInfo{1} = strrep(theoryStruct{comparisonStruct{nr1}.idx}.name,'_','\_');
    legendInfo{2} = strcat(['$C_{' strrep(barcodeGen{nr1}.name,'_','\_') '}$ = ' num2str(max1,'%.4f ')]);

     for ii=2:length(selectedIndices)   
         % if the best theories are not matching, not much sense in
         % placing them on the theory together
        [ ~,~, bar2, bit2] = load_theory_and_stretch_ex(selectedIndices(ii), theoryStruct, comparisonStruct,barcodeGen );


        bar2(~bit2) = nan;

        % start positions
        pos2 =  comparisonStruct{selectedIndices(ii)}.pos(1);
        % max coefficients
        max2 = comparisonStruct{selectedIndices(ii)}.maxcoef(1);
        % orientations
        or2 = comparisonStruct{selectedIndices(ii)}.or(1);
        
        if or2 == 2
            bar2 = fliplr(bar2);
        end

        plot(pos2:pos2+length(bar2)-1,(bar2-nanmean(bar2))./nanstd(bar2))
        if ii< 5
            legendInfo{ii+1} =  strcat(['$C_{' strrep(barcodeGen{selectedIndices(ii)}.name,'_','\_') '}$ = ' num2str(max2,'%.4f ')]);
        end
    end
	legend(legendInfo,'Interpreter','latex','location', 's');

 end
        

end

