function [maxPCC,totPlacementAttempts] = gen_scores_real(...
    bg,bg2,rlLeft,rlRight,minOverlapLen,...
    NN,sF,out,gap,nTries, SCAMP_LINE)

    %   Args:
    %       bg - barcodes from first group
    %       bg2 - barcodes from second group
    %       rlLeft - left random length
    %       rlRight - right random length
    %       minOverlapLen - min overlap length,...
    %       NN,sF,out,SCAMP_LINE

    %   Returns:
    %
    if nargin < 9
        gap = 100;
        nTries = 1;%length(namesBar);
    end

    if nargin < 11
        SCAMP_LINE = 'SCAMP';
    end

    
    mkdir(out);
    mp1 = cell(1,NN);
    maxPCC=cell(1,NN);
    totPlacementAttempts = zeros(1,NN);
    
    barLensLeft = cellfun(@(x) sum(x.rawBitmask),bg);
    barLensRight = cellfun(@(x) sum(x.rawBitmask),bg2);
    
    goodLensLeft = barLensLeft>rlLeft+(NN-1)*gap; % barcodes longer than length required

    goodLensRight = barLensRight>rlRight; % barcodes longer than length required
    
    nTries =min(nTries,sum(goodLensRight));
for k=1:NN
    curLen = rlLeft+(k-1)*gap;

    maxPCC{k} = cell(1,nTries);
    
    % create bar structure with random cut-out from the original barcode
    barcodes = cellfun(@(x) x.rawBarcode(x.rawBitmask),bg(goodLensLeft),'un',false);
    barcodes = cellfun(@(x) x(randi(length(x)-curLen)+(1:curLen)),barcodes,'un',false); % random cut-out of length curLen
    bitmasks = cellfun(@(x) true(size(x)), barcodes, 'un', 0);
    barSynth = cell2struct([barcodes; bitmasks]',{'rawBarcode','rawBitmask'},2);
    
    barcodes2 = cellfun(@(x) x.rawBarcode(x.rawBitmask),bg2(goodLensRight),'un',false);
    barcodes2 = cellfun(@(x) x(randi(length(x)-rlRight)+(1:rlRight)),barcodes2,'un',false);
    bitmasks2 = cellfun(@(x) true(size(x)), barcodes2, 'un', 0);
    barSynth2 = cell2struct([barcodes2; bitmasks2]',{'rawBarcode','rawBitmask'},2);

    % have to save separately..
    foldSynth= 'synthtest';
    % have to save separately..
    [namesBar, stridx] = Core.save_bars_rescaled_txt(barSynth2,sF,foldSynth);

    % now save the same barcodes as long vectors since we will compare against
    % these.\todo: possibiility of using one long barcode?
    [names2, baridx2] = Core.save_long_bar_txt_skip(barSynth,inf,foldSynth);

    for ii=1:nTries
        command = strcat([SCAMP_LINE ' --window=' num2str(minOverlapLen+(ii-1)*gap) ' --input_a_file_name='...
           names2{1} ' --input_b_file_name=' namesBar{ii} ' --num_cpu_workers=28 --no_gpu --output_pearson --print_debug_info'...
           ' --output_a_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
           ' --output_a_index_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
            [status,cmdout] = system(command);
        % now extract from the output best score for each comparison
        fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
        raw2 = textscan(fid, '%s ');
        fclose(fid);
        nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
        mp1{k}{ii} = nan(length(nonanValues),1);
        mp1{k}{ii}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');


%         mp1{k}{ii} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
        %mpI1{k}{ii} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));  
        [PKS,LOCS] = sort(mp1{k}{ii},'descend','MissingPlacement','last');
        LOCS = LOCS(~isnan(PKS));
        PKS = PKS(~isnan(PKS));
        barNum = baridx2{1}(LOCS);

        [LOCSUIQUE,LOCSUNIQUEPos] = unique(barNum,'stable');
        maxPCC{k}{ii}(LOCSUIQUE) = PKS(LOCSUNIQUEPos);
% 
%         PKS(LOCSUNIQUEPos(LOCSUIQUE==2))
%         max(mp1{k}{ii}( baridx2{1}==2))

    end
    
    % total number of placement attempts (including sF and reverse barcode)
    totPlacementAttempts(k) = length(stridx{1})-(2*length(sF)-1)*minOverlapLen;


%     toc

%     % find maxpcc for each barcode
%     [PKS,LOCS] = sort(mp1{k},'descend');
%     LOCS = LOCS(~isnan(PKS));
%     PKS = PKS(~isnan(PKS));
%     barNum = baridx2{1}(LOCS);
%     [LOCSUIQUE,LOCSUNIQUEPos] = unique(barNum,'stable');
%     maxPCC{k}(LOCSUIQUE) = PKS(LOCSUNIQUEPos);
%     totLen(k) = length(stridx{1})-length(sF)*minOverlapLen;

end


end

