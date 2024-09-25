function [valsPassingThresh,pccValsMax,mp1,mpI1] = find_periodic_barcodes(barStruct,timestamp,numWorkers)

    if nargin < 3
        numWorkers = 5;
        thresh = 0.9;
        MIN_OVERLAP_PIXELS = 50;

    end
    
    sF = 1; % will check self similarity
    foldSynth = timestamp;
    mkdir(foldSynth);
    [namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth,0);

%     [namesBar, stridx, baridx] = Core.save_bars_rescaled_txt_stack(barStruct,sF,foldSynth);
%     stridx = stridx(1:end-MIN_OVERLAP_PIXELS+1);
%     baridx = baridx(1:end-MIN_OVERLAP_PIXELS+1);
    %%
    mpI1  = [];
    mp1 = [];
    %     numWorkers = 4;
    for k=1:length(namesBar)
        k

        if ispc
            command = strcat(['scamp --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
            namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
            ' --output_a_file_name=.\' fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
            ' --output_a_index_file_name=.\' fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
        end
        
        [~,message] = system(command);
        
        mpI1{k} = importdata(fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_index' ])));
    
    %also get PCC scores.
%     tic
    % a bit slower import to make sure imports correctly on win
    fid = fopen(fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_mp' ])));
    raw2 = textscan(fid, '%s ');
    fclose(fid);
    nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
    mp1{k} = nan(length(nonanValues),1);
    mp1{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
%     toc


    end
    %%
    
    % maximum values,, 2 max, since two locations should be
    % similar
    
    sortedMP = cellfun(@(x) sort(x,'desc','MissingPlacement','last'),mp1,'un',false);
    
    pccValsMax = cellfun(@(x) mean(x(1:2)), sortedMP);
    
    valsPassingThresh = (pccValsMax<thresh);
    
    
end

% 
% 
% function [mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE,MIN_OVERLAP_PIXELS,numWorkers,names2,namesBar)
% % parpool(c)
% NN=length(namesBar);
% import Zeromodel.beta_ev_cdf;
% 
% out=strcat('output',timestamp);
% % mkdir(out);
% % parpool('local',28)
% mp1 = cell(1,NN);
% % mp2 = cell(1,NN);
% % PKS1=cell(1,NN);
% % LOCS1=cell(1,NN);
% % pksUnique1 = cell(1,NN);
% % pksUniquePos1=cell(1,NN);
% % barcodePair1=cell(1,NN);
% % rescaleFactorPair1=cell(1,NN);
% % orPair1=cell(1,NN);
% maxMP = zeros(1,NN);
% % pvals =zeros(1,NN);
% % tic
% for k=1:NN
% %     k
%     % win
%     if ispc
%     command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
%        names2{k} ' --input_b_file_name=.\' namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
%        ' --output_a_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
%        ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
%     else
%        command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
%        names2{k} ' --input_b_file_name=\' namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
%        ' --output_a_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
%        ' --output_a_index_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
%     end
%     [~,message] = system(command);
% %     tic
% %     mp1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
% %     toc
%     mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));
%     
%     %also get PCC scores.
% %     tic
%     % a bit slower import to make sure imports correctly on win
%     fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
%     raw2 = textscan(fid, '%s ');
%     fclose(fid);
%     nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
%     mp1{k} = nan(length(nonanValues),1);
%     mp1{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
% %     toc
%    
%     
%     maxMP(k) = max(mp1{k});
%     
% %     [PKS1{k},LOCS1{k},pksUnique1{k},pksUniquePos1{k},barcodePair1{k},rescaleFactorPair1{k},orPair1{k}] = ...
% %     Core.unique_peak_barcodes(mp1{k},mpI1{k}, 0.3,baridx2,stridx,k);
% % 
% % %     % will also compute p-value
% %     a = 0.122*MIN_OVERLAP_PIXELS;
% %     len2 = lengths(barcodePair1{k}(1));
% %     len1 = lengths(k);
% %     n = 2*(len2-MIN_OVERLAP_PIXELS)*(len1-MIN_OVERLAP_PIXELS)/100;
% % % 
% %     pvals(k) = 1-Zeromodel.beta_ev_cdf(maxMP(k), a, 1, n, 1);
% % %
% 
%     
% end
% % delete(gcp('nocreate'))
% 
% 
% end
% 
