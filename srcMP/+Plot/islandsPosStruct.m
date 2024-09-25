
function comparisonStruct = islandsPosStruct(barcodeIslandsData,barcodeIslands)
    % islandsPosStruct
    %
    %   Args:
    %       barcodeIslandsData 
    %       barcodeIslands
    %   Returns:
    %       comparisonStruct - comparison structure for a specific
    %       consensus cluster


    comparisonStruct = cell(1,size(cell2mat(barcodeIslandsData'),1));
    t= 1;
    for i=1:length(barcodeIslandsData)
        for j=1:size(barcodeIslandsData{i},1)
            comparisonStruct{t}.maxcoef = 0.99; % for now fix this. TODO: For experimental data, pass as an argument/from data. Can calculate coef wrt consensus later
            comparisonStruct{t}.pos = barcodeIslandsData{i}(j,1);
            comparisonStruct{t}.or = barcodeIslandsData{i}(j,3);
            comparisonStruct{t}.bestBarStretch = barcodeIslandsData{i}(j,4);

            comparisonStruct{t}.lengthMatch = barcodeIslandsData{i}(j,2)-barcodeIslandsData{i}(j,1)+1;
            comparisonStruct{t}.idx = i;
            comparisonStruct{t}.barid = barcodeIslands{i}(j);
            t = t+1;
        end
    end


    
end