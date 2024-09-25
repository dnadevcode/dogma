function [  ] = export_exp_vs_theory(comparisonStructAll,theoryStruct,barcodeGenC,sets,barNr,theoryNr  )
    % input  comparisonStructAll,theoryStruct,barcodeGenC,sets

%     barNr = 3;
%     theoryNr = 2;

    for i=1:length(comparisonStructAll)
        for j=1:length(comparisonStructAll{1})
            try
                comparisonStructAll{i}{j}.idx = i;
                comparisonStructAll{i}{j}.name = theoryStruct{i}.name;
            catch
            end  
        end
    end
    CBT.Hca.Export.plot_comparison_vs_theory(comparisonStructAll{theoryNr},theoryStruct,barcodeGenC,barNr,sets.export.savetxt);




end

