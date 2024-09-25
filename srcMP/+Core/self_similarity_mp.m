function [MPIthr1, MPthr] = self_similarity_mp(w, theoryStructRev, numWorkers)
    SCAMP_LINE = 'SCAMP';

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    out = ['output_',timestamp];
    [status,msg] = mkdir(out);

    %
    MPIthr1 = cell(1,length(theoryStructRev));
    MPthr = cell(1,length(theoryStructRev));
    for idxTheory = 1:length(theoryStructRev)
        k=1;
        if ispc
            command = strcat([SCAMP_LINE ' --window=' num2str(w) ' --input_a_file_name=.\'...
                theoryStructRev{idxTheory}.filename ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
                ' --output_a_file_name=.\' fullfile(out,strcat([num2str(idxTheory) '_' num2str(k) '_mp' ])) ...
                ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(idxTheory) '_' num2str(k) '_index' ])) ]);
        else

            command = strcat([SCAMP_LINE ' --window=' num2str(w) ' --input_a_file_name=\'...
                theoryStructRev{idxTheory}.filename   ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
                ' --output_a_file_name=\' fullfile(out,strcat([num2str(idxTheory) '_' num2str(k) '_mp' ])) ...
                ' --output_a_index_file_name=\' fullfile(out,strcat([num2str(idxTheory) '_' num2str(k) '_index' ])) ]);
        end
        [~,message] = system(command); %--output_pearson

        MPIthr1{idxTheory} = importdata(fullfile(out,strcat([num2str(idxTheory) '_' num2str(k) '_index' ])));

        %also get PCC scores.
        fid = fopen(fullfile(out,strcat([num2str(idxTheory) '_' num2str(k) '_mp' ])));
        raw2 = textscan(fid, '%s ');
        fclose(fid);
        nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
        MPthr{idxTheory} = nan(length(nonanValues),1);
        MPthr{idxTheory}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
%         thr{idxTheory} = importdata(   theoryStructRev{idxTheory}.filename);
    end

        [~,~] =  rmdir(out,'s'); % remove temporary folder
end