function [val,sorti,valRes] = scaled_pdif_vs_theory(bT,theoryStruct,barcodeGen,barIslands,barcodeIslandsData,sF,sets)

% scaled_pdif_vs_theory - place barcode from the barcode island against the
% theory
%   Args:
%       bT - theory
%       barcodeGen - exp barcodes
%       barIslands - island elts
%       barcodeIslandsData - table
%       sF - stretch factors for comparing against thry

%   Return:
%       valRes


% Here explain the re-scaling strategy that is used.

% This is a validation method.
% Here we're comparing theory against exeriments to get PDIF. There are a
% number of methods that one can use to compare experiments vs. theory.


import Core.update_sf_barset;
import Core.synth_to_table;

%% Validation. For all islands
theoryStruct.rawBitmask = bT{1}.rawBitmask;
% tS.length = length(bT{1}.rawBitmask);
pDif = cell(1,length(barIslands));
sorti =  cell(1,length(barIslands));
val =  cell(1,length(barIslands));

valRes = cell(1,length(barIslands));
for idx = 1:length(barIslands)

    %%1) Run comparison (all barcodes)
    bars = barcodeGen(barIslands{idx}); 
    barx = barIslands{idx};

    % Compare individual barcodes: todo: change to hca barcode
    % comparison
    [compI, rezI, ~] = compare_to_t(bars,theoryStruct,sF,sets); % maybe calculate Stouffer score here also?

    % get the barcode with minimum p-value
    [a] = cellfun(@(x) x.pval,compI);
    [minV, selBarId] = min(a);
    %
    
    % In practice, we might want to choose another barcode to anchor from
    % (i.e. if we're looking at specific place along the genome)
    import Validation.scaled_specific_barcode_vs_theory;
    [valRes{idx},val{idx},sorti{idx}] = scaled_specific_barcode_vs_theory(compI,selBarId,barcodeIslandsData,theoryStruct,idx);
    valRes{idx}.compI = compI;
    valRes{idx}.bars = bars;



end


end

