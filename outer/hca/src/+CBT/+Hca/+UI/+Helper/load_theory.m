function [ theoryStruct ] = load_theory( sets )
    % load_theory
    % Loads theory names from path
    %     Args:
    %         sets: Input settings for loading theory
    % 
    %     Returns:
    %         theoryStruct: theory structure
    % 
    %     Example:
    %          
    
    % todo: make it recognize that the file list is txt files,
    % and then just save these names to txt struct
    
    % In case theory is in a single file, read of it's contents into a
    % folder and then save the names of the theories
    
    theoryStruct = {};

    import CBT.Hca.UI.Helper.load_theory_mat;
    import CBT.Hca.UI.Helper.load_theory_into_txts;

    % loop through theories
    for idx=1:length(sets.theoryFile)
        [~,part,st] = fileparts(sets.theoryFile{idx});
        if isequal(st,'.txt')
            % file already exists, just load the info from header
            names = fullfile(sets.theoryFileFold{idx},sets.theoryFile{idx});
            parts = strsplit(names,'_');
            theoryStructN.filename = names;
            theoryStructN.name = part;
            theoryStructN.length = round(str2num(parts{end-5}));
            theoryStructN.meanBpExt_nm  = str2num(parts{end-4});
            theoryStructN.pixelWidth_nm =  str2num(parts{end-3});
            theoryStructN.psfSigmaWidth_nm = str2num(parts{end-2});
            theoryStructN.isLinearTF = str2num(parts{end-1});
            theoryStruct = [theoryStruct theoryStructN];
        else
            % all theories
           [bars, bits, names, meanbpnm, pixelWidth_nm, psfSigmaWidth_nm, isLinearTF] = ...
            load_theory_mat(sets,idx);

            if ~sets.theory.theoryDontSaveTxts
                [theoryStruct] = load_theory_into_txts(theoryStruct, sets.theoryFileFold{idx}, bars,bits, names, meanbpnm, pixelWidth_nm, psfSigmaWidth_nm, isLinearTF,sets.theory.precision);
            else
                % add info about barcodes as struct
                theoryStruct = cell2struct([bars;...
                bits;arrayfun(@(x) isLinearTF,1:length(bars),'un',false);...
                names;...
                cellfun(@(x) length(x),bars,'un',false);...
                arrayfun(@(x) meanbpnm,1:length(bars),'un',false);...
                arrayfun(@(x) pixelWidth_nm,1:length(bars),'un',false);...
                arrayfun(@(x) psfSigmaWidth_nm,1:length(bars),'un',false)]',...
                {'rawBarcode','rawBitmask', 'isLinearTF','name','length','meanBpExt_nm','pixelWidth_nm','psfSigmaWidth_nm'},2);
            end
        end
    end
end

