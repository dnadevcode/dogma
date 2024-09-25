function [mpI1,mp1,maxMP] = calc_overlap(barStruct, foldSynth, SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers, names2, namesBar)
    % calc_overlap
    
    %   Args:
    %       MIN_OVERLAP_PIXELS
    %       names2
    %       namesBar

    %   Returns:
    %       mpI1 - matrix profile index
    %       mp1 - matrix profile
    %       maxMP - maximum of matrix profile

    if isempty(SCAMP_LINE)
        SCAMP_LINE = 'SCAMP';
    end

    NN=length(namesBar);
    
    out = foldSynth;
    mp1 = cell(1,NN);
    mpI1 = cell(1,NN);
    maxMP = zeros(1,NN);
    command = cell(1,NN);
    for k=1:NN
        % win
        if ispc
        command{k} = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
           names2{k} ' --input_b_file_name=.\' namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
           ' --output_a_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
           ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
        else
           command{k} = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
           names2{k} ' --input_b_file_name=\' namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
           ' --output_a_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
           ' --output_a_index_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
        end
    end

%     if ispc
        for k=1:NN
            [~, message] = system(command{k}); % single command so we don't need to go in and out of matlab
        end
        
%     else
%         [~,message] = system(strjoin(command, ';')); % single command so we don't need to go in and out of matlab
%     end


    for k=1:NN
        % todo: if scamp didn't work, we need to run matlab's version

        
        %     tic
        %     mp1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
        %     toc
        mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));
        
        % a bit slower import to make sure imports correctly on win
        fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
        raw2 = textscan(fid, '%s ');
        fclose(fid);
        nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
        mp1{k} = nan(length(nonanValues),1);
        mp1{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
      
        maxMP(k) = max(mp1{k});  
    end

end