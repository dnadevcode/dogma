function [maxPCC,totLen,barcodes2] = gen_scores(...
    NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,RAND_LENGTH_2,MIN_OVERLAP_PIXELS,...
    NN,sF,out,SCAMP_LINE)


% NN=100;
mkdir(out);
% parpool('local',28)
mp1 = cell(1,NN);
mp2 = cell(1,NN);
PKS1=cell(1,NN);
LOCS1=cell(1,NN);
pksUnique1 = cell(1,NN);
pksUniquePos1=cell(1,NN);
barcodePair1=cell(1,NN);
maxPCC=cell(1,NN);
orPair1=cell(1,NN);
maxMP = zeros(1,NN);
totLen = zeros(1,NN);

import Nullmodel.gen_random;

for k=1:NN
    % run loop
%     tic
    k
    maxPCC{k} = nan(1,NUM_RAND_FRAGMENTS);
    

    [barcodes] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN+(k-1)*10,0);
    %
    bitmasks = cellfun(@(x) true(size(x)), [barcodes], 'un', 0);


    % create bar structure 
    barSynth = cell2struct([barcodes; bitmasks]',{'rawBarcode','rawBitmask'},2);


    [barcodes2] = gen_random(1,PSF_WIDTH_PIXELS,RAND_LENGTH_2,0);
    bitmasks2 = cellfun(@(x) true(size(x)), [barcodes2], 'un', 0);

    barSynth2 = cell2struct([barcodes2; bitmasks2]',{'rawBarcode','rawBitmask'},2);


    % barSynth = [cell2struct(barcodes,{'rawBarcode'},1);cell2struct(bitmasks,{'bitmasks'},1)];


    % have to save separately..
    foldSynth= 'synthtest';
    % have to save separately..
    [namesBar, stridx] = Core.save_bars_rescaled_txt(barSynth2,sF,foldSynth);

    % now save the same barcodes as long vectors since we will compare against
    % these.\todo: possibiility of using one long barcode?


    [names2, baridx2] = Core.save_long_bar_txt_skip(barSynth,inf,foldSynth);



%     tic
    command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name='...
       names2{1} ' --input_b_file_name=' namesBar{1} ' --num_cpu_workers=30 --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);

    [status,cmdout] = system(command);
%     toc
%     mp1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
    mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));
    fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
    raw2 = textscan(fid, '%s ');
    fclose(fid);
    nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
    mp1{k} = nan(length(nonanValues),1);
    mp1{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
%     toc

    % find maxpcc for each barcode
    [PKS,LOCS] = sort(mp1{k},'descend');
    LOCS = LOCS(~isnan(PKS));
    PKS = PKS(~isnan(PKS));
    barNum = baridx2{1}(LOCS);
    [LOCSUIQUE,LOCSUNIQUEPos] = unique(barNum,'stable');
    maxPCC{k}(LOCSUIQUE) = PKS(LOCSUNIQUEPos);
    totLen(k) = length(stridx{1})-length(sF)*MIN_OVERLAP_PIXELS;

% 
%  
%     toc
end


end

