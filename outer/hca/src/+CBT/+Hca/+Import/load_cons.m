function [data] = load_cons(uniqueNames,uniqueFold,outfold)
    % Args:
    %   fold,sample1
    % Returns:
    %   A, B
    %
    % all consensuses
    % fold = '/media/albyback/My Passport/DATA/SV/mbio_sv/New labels_220kbp all barcodes/*.mat';
%     fold = '/home/albyback/git/sv/data/mbio_sv/New labels_220kbp all barcodes/*.mat';
%     files = dir(fold);
    if nargin < 3
        outfold = 'outdata';
    end

    % take the first sample
    % sample1 = {'P11K0','P11K22'}; 
%     sample1 = {'P6K0','P6K25'}; 
% 
%     settings = 'sample1_settings.txt';

    % camera res - unknown (. nm/bp stretch for each sample - unknown. Hence
    % convertion px to bp - unknown. Henceforth need stretching to get the best!

    % load data
    data = cell(1,length(uniqueNames));
    for i=1:length(uniqueNames)
        data{i} = load(fullfile(uniqueFold{i},uniqueNames{i}));
        
           % load first
        fname = fullfile(outfold,strcat([ uniqueNames{i} '.txt']));
        fileID = fopen(fname,'w');
        fprintf(fileID,'%2.16f ',data{i}.clusterConsensusData.barcode);
        fclose(fileID);
        data{i}.fname = fname;
        data{i}.name = fname;
        data{i}.rawBarcode = data{i}.clusterConsensusData.barcode;
        data{i}.rawBitmask = data{i}.clusterConsensusData.bitmask;
        
        % now make sure that if there is bitmask, all zeros are shifted to
        % be at the end of the barcode (alt equally spaced..)
        zerosLoc = find(data{i}.rawBitmask==0,1,'last');
        if ~isempty(zerosLoc)
            data{i}.rawBarcode = circshift(data{i}.rawBarcode,[0,-zerosLoc]);
            data{i}.rawBitmask = circshift(data{i}.rawBitmask,[0,-zerosLoc]);
        end
     
        %TODO : add check to see if there are a few reasons where bitmask
        %is 0
        
         data{i}.circ = isempty(find(data{i}.clusterConsensusData.bitmask==0));
%         iscirculard = isempty(find(data{end}.clusterConsensusData.bitmask==0));

        % set if this comparison is circular
%         sets.circ = iscircularq;

        % should we resize to larger size so that rescaling would have smaller affect?
        %         data{i}.barcode = data{i}.clusterConsensusData.barcode;
        %         data{i}.bitmask = data{i}.clusterConsensusData.bitmask;
    end

end

