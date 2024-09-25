function [mpI1,mp1,maxMP,mpI1B,mp1B] = calc_overlap_both(barStruct,timestamp,SCAMP_LINE,MIN_OVERLAP_PIXELS,numWorkers,names2,namesBar)
% parpool(c)
NN=length(namesBar);
import Zeromodel.beta_ev_cdf;

out=strcat('output',timestamp);
% mkdir(out);
% parpool('local',28)
mp1 = cell(1,NN);
% mp2 = cell(1,NN);
% PKS1=cell(1,NN);
% LOCS1=cell(1,NN);
% pksUnique1 = cell(1,NN);
% pksUniquePos1=cell(1,NN);
% barcodePair1=cell(1,NN);
% rescaleFactorPair1=cell(1,NN);
% orPair1=cell(1,NN);
maxMP = zeros(1,NN);
% pvals =zeros(1,NN);
% tic
for k=1:NN
%     k
    % win
    if ispc
    command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
       names2{k} ' --input_b_file_name=.\' namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
    else
       command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
       names2{k} ' --keep_rows --input_b_file_name=\' namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ...
        ' --output_b_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mpB' ])) ...
       ' --output_b_index_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_indexB' ]))]);
    end
    [~,message] = system(command);
%     tic
%     mp1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
%     toc
    mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));
    
    mpI1B{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_indexB' ])));

    %also get PCC scores.
%     tic
    % a bit slower import to make sure imports correctly on win
    fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
    raw2 = textscan(fid, '%s ');
    fclose(fid);
    nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
    mp1{k} = nan(length(nonanValues),1);
    mp1{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
%     toc
   
    
    maxMP(k) = max(mp1{k});
    
%     [PKS1{k},LOCS1{k},pksUnique1{k},pksUniquePos1{k},barcodePair1{k},rescaleFactorPair1{k},orPair1{k}] = ...
%     Core.unique_peak_barcodes(mp1{k},mpI1{k}, 0.3,baridx2,stridx,k);
% 
% %     % will also compute p-value
%     a = 0.122*MIN_OVERLAP_PIXELS;
%     len2 = lengths(barcodePair1{k}(1));
%     len1 = lengths(k);
%     n = 2*(len2-MIN_OVERLAP_PIXELS)*(len1-MIN_OVERLAP_PIXELS)/100;
% % 
%     pvals(k) = 1-Zeromodel.beta_ev_cdf(maxMP(k), a, 1, n, 1);
% %
    fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mpB' ])));
    raw2 = textscan(fid, '%s ');
    fclose(fid);
    nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
    mp1B{k} = nan(length(nonanValues),1);
    mp1B{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
%
    
end
% delete(gcp('nocreate'))


end

