function [consensusBar] = create_consensus(dataTable, bars)
    % create consensus from aligned data table

    % sort barcodes based on starting position, so the next barcode would
    % have an overlap with current barcode

    % stretch barcodes to correct lengths
    barsC = bars;
    % convert bars
    for i=1:length(barsC)
        lenBarTested = length( barsC{i}.rawBarcode );
        barsC{i}.rawBarcode = interp1(barsC{i}.rawBarcode, linspace(1,lenBarTested,lenBarTested*dataTable(i,4)));
        barsC{i}.rawBitmask =  barsC{i}.rawBitmask (round(linspace(1,lenBarTested,lenBarTested*dataTable(i,4))));

%         barsC{i}.rawBarcode = barsC{i}.rawBarcode( barsC{i}.rawBitmask);

         if dataTable(i,3)==2
             barsC{i}.rawBarcode = fliplr(barsC{i}.rawBarcode);
             barsC{i}.rawBitmask = fliplr(barsC{i}.rawBitmask);
         end         
    end

    % sorted values, for easier consensus calculation
    [val,id] = sort(dataTable(:,1));


    barZscored = cell(1,length(bars));

    ii=1;
    m = mean(barsC{id(ii)}.rawBarcode(barsC{id(ii)}.rawBitmask));
    s = std(barsC{id(ii)}.rawBarcode(barsC{id(ii)}.rawBitmask),1);
    barZscored{1}.rawBarcode = (barsC{id(ii)}.rawBarcode-m)./s;
    barZscored{1}.rawBitmask = barsC{id(ii)}.rawBitmask;

    barZscored{1}.rawBarcode((barZscored{1}.rawBitmask~=1)) = nan;
    for ii = 2:length(bars)
%     dataTable(id(1),1:2)

        matchStart = dataTable(id(ii),1);
        matchStop = min(dataTable(id(ii-1),2),dataTable(id(ii),2));
    
        b1 =  barZscored{ii-1}.rawBarcode;
        b1(~barZscored{ii-1}.rawBitmask) = nan;
    
        b2 = barsC{id(ii)}.rawBarcode;
        b2(~barsC{id(ii)}.rawBitmask) = nan;
    
        valsb2 = b2(1:matchStop-matchStart+1);
        valsb1 = b1(matchStart-dataTable(id(ii-1),1):matchStart-dataTable(id(ii-1),1)+length(valsb2)-1);
    
        mask = ~isnan(valsb2.*valsb1);
    
        mean1 = mean(valsb1(mask));
        std1 = std(valsb1(mask),1);
    
        mean2 = mean(valsb2(mask));
        std2 = std(valsb2(mask),1); 

        barZscored{ii}.rawBarcode = (b2-mean2)./std2*std1+mean1;
        barZscored{ii}.rawBitmask = barsC{id(ii)}.rawBitmask;
    end


    
    % now just need to average the z-scored barcodes
    maxV = val+cellfun(@(x) length(x.rawBarcode),barZscored)';
    minV = min(maxV);
    len = max(maxV-minV+1);

    consensusBar = zeros(1,len);
%     for i=1:length(val);
%         con
% 
%     end

    figure,plot(val(1):val(1)+length(barZscored{1}.rawBarcode)-1,barZscored{1}.rawBarcode)
    hold on
%     figure
    ix1 = 2;
    plot(val(ix1):val(ix1)+length(barZscored{ix1}.rawBarcode)-1,barZscored{ix1}.rawBarcode)
        ix1 = 3;
    plot(val(ix1):val(ix1)+length(barZscored{ix1}.rawBarcode)-1,barZscored{ix1}.rawBarcode)

        ix1 = 4;
    plot(val(ix1):val(ix1)+length(barZscored{ix1}.rawBarcode)-1,barZscored{ix1}.rawBarcode)

        ix1 = 5;
    plot(val(ix1):val(ix1)+length(barZscored{ix1}.rawBarcode)-1,barZscored{ix1}.rawBarcode)

% 
%    hold on
%     val = 3;
%     plot(id(val):id(val)+length(barZscored{val}.rawBarcode)-1,barZscored{val}.rawBarcode)
%    hold on
%     val = 4;
%     plot(id(val):id(val)+length(barZscored{val}.rawBarcode)-1,barZscored{val}.rawBarcode)
%    hold on
%     val = 5;
%     plot(id(val):id(val)+length(barZscored{val}.rawBarcode)-1,barZscored{val}.rawBarcode)
% 






end

