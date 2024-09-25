function [mpI1,mp1,maxMP] = calc_overlap_first(barStruct,timestamp,SCAMP_LINE,MIN_OVERLAP_PIXELS,numWorkers,names2)
    % parpool(c)rcalc_overlap_neighbors
    %%
    % calc_overlap_neighbors - calculates MP scores for neighboring barcodes
    NN=length(barStruct)-1;
    import Zeromodel.beta_ev_cdf;

    out=strcat('output',timestamp);
    % mkdir(out);
    % parpool('local',28)
    mp1 = cell(1,NN);

    maxMP = zeros(1,NN);
    tic
    for k=1:NN
    %     k
        % win
        if ispc
        command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
           names2{1} ' --input_b_file_name=.\' names2{k+1} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
           ' --output_a_file_name=.\' fullfile(out,strcat([num2str(1) '_' num2str(k+1) '_mp' ])) ...
           ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(1) '_' num2str(k+1) '_index' ])) ]);
        else
           command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
           names2{1} ' --input_b_file_name=\' names2{k+1} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
           ' --output_a_file_name=\' fullfile(out,strcat([num2str(1) '_' num2str(k+1) '_mp' ])) ...
           ' --output_a_index_file_name=\' fullfile(out,strcat([num2str(1) '_' num2str(k+1) '_index' ])) ]);
        end
        [~,message] = system(command);

        mpI1{k} = importdata(fullfile(out,strcat([num2str(1) '_' num2str(k+1) '_index' ])));

        %also get PCC scores.
    %     tic
        % a bit slower import to make sure imports correctly on win
        fid = fopen(fullfile(out,strcat([num2str(1) '_' num2str(k+1) '_mp' ])));
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


    end
    % delete(gcp('nocreate'))


end

