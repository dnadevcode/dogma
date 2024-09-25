function [ kymoStructs ] = extract_kymos_structs(rawKymos, rawKymoFilepaths, pixelsWidths_bps,alignedKymos,alignedKymosStretchFactors,shiftAlignedKymos)
    if nargin < 3
        passesFilters = ones(length(rawKymos),1);
        alignedKymos = cell(length(rawKymos),1);
        alignedKymosStretchFactors = cell(length(rawKymos),1);
        shiftAlignedKymos = cell(length(rawKymos),1);
        pixelsWidths_bps = nan(length(rawKymoFilepaths),1);

    end
                pixelsWidths_bps = -ones(length(rawKymoFilepaths),1);

    dataType = {'Kymo'};
    import Fancy.Utils.data_hash;
    dataHashes = cellfun(@data_hash, rawKymos, 'UniformOutput', false);

    pixelsWidths_bps =  arrayfun(@(tf) tf,pixelsWidths_bps, 'UniformOutput', false);
    passesFilters = arrayfun(@(tf) tf,passesFilters, 'UniformOutput', false);
    
    rawKymoFileMoleculeIdxs =transpose(1:length(rawKymos));
   %     rawKymoFileMoleculeIdxs = arrayfun(@(tf) tf,rawKymoFileMoleculeIdxs, 'UniformOutput', false);
    rawKymoFilepathsShort = cellfun(@(tf) tf(max(1,1):end), rawKymoFilepaths,'UniformOutput', false)
    displayName = strcat(arrayfun(@(fileMoleculeIdx) ['Mol # ', num2str(fileMoleculeIdx) ',' ' '], rawKymoFileMoleculeIdxs, 'UniformOutput', false),rawKymoFilepathsShort);

    kymoStructs = struct(...
        'passesFilters', passesFilters,...
        'filePath', rawKymoFilepaths,...
        'unalignedKymo', rawKymos,...
        'dataHash', dataHashes,...
        'alignedKymos', alignedKymos,...
        'displayName',displayName,...
        'alignedKymosStretchFactors', alignedKymosStretchFactors,...
        'shiftAlignedKymos', shiftAlignedKymos,...
        'bpsPerPx', pixelsWidths_bps,...
        'type', dataType...
    );
    kymoStructs = arrayfun(@(kymoStruct) kymoStruct, kymoStructs, 'UniformOutput', false);
end

