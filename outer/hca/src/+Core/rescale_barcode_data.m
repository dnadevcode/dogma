function [barGen] = rescale_barcode_data(barGen,stretchFactors,initialStretch, maxRescaleFactor)

%   Args:
%       barGen - barcode cell structure
%       stretchFactors - stretch factors
%       initialStretch - initial length re-scaling factor (if i.e. barcodes
%       come from different datasets)
%   Returns:
%       barGen - with included rescaled struct

if nargin < 3
    initialStretch = ones(1,length(barGen));
end

% if nargin < 4
%     maxRescaleFactor = 0;
% end

for i=1:length(barGen)
    lenBarTested = length(barGen{i}.rawBarcode);

    for j=1:length(stretchFactors)
        barGen{i}.rescaled{j}.rawBarcode =  interp1(barGen{i}.rawBarcode, linspace(1,lenBarTested,lenBarTested*initialStretch(i)*stretchFactors(j)));
        barGen{i}.rescaled{j}.rawBitmask  = barGen{i}.rawBitmask(round(linspace(1,lenBarTested,lenBarTested*initialStretch(i)*stretchFactors(j))));
    end

end

end

