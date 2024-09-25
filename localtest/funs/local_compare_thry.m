function [MPthr, MPIthr1] = local_compare_thry(w,fastaFile, nmbp, nmPerPx, psffac,numWorkers,theoryStructRev)
    % Compare theories locally to find regions that are matched worse
    
    %   Args:
    %       fastaFile - fasta files
    %       MIN_OVERLAP_PIXELS - minimum overlap pixels
    %       nmbp - nanometer per basepar
    %       nmPerPx - nanometer per pixel (camera)
    %       psffac - point spread function scaling factor

    %   Returns:
    %       MPthr
    %       MPIthr1

    if nargin < 1
        w = 150;
    end

    if nargin < 2
        fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta','/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta','sequence1.fasta','sequence2.fasta','sequence3.fasta'};
        nmPerPx = 261;
        nmbp = 0.225;
        psffac = 1;
        numWorkers = 4;
    end

    if nargin < 7
        % generate theoretical
        theoryStructRev = Thry.gen_save_theoretical(nmPerPx, psffac, fastaFile, nmbp);
    end

    % calculate self-similarity using window width w
    [MPIthr1, MPthr] = Core.self_similarity_mp(w, theoryStructRev, numWorkers);

    thr{1} = importdata(   theoryStructRev{1}.filename);


%     mpAll = MPthr{1}(1:end/2);
    
    %
    % barStruct = [];
    % barStruct(1).rawBarcode = thr{1};
    % barStruct(1).rawBitmask = ones(1,length(barStruct(1).rawBarcode ));
    %
    % barStruct(2).rawBarcode = thr{1};
    % barStruct(2).rawBitmask = ones(1,length(barStruct(2).rawBarcode ));
    %
    % stridx = {[ones(1,theoryStruct{1}.length) nan -ones(1,theoryStruct{1}.length)]};
    % import Core.mp_res_to_struct;
    % [overlapStruct] = mp_res_to_struct({MPthr{1},[]},{MPIthr1{1},[]},{2*ones(1,length(MPthr{1}))},stridx,MIN_OVERLAP_PIXELS,1,barStruct)


    plotex = 0;
    if plotex
        %% plot
        figure
        plot(MPthr{1})

        figure
        plot(MPIthr1{1})
        % tiledlayout(3,2)
        % for i=1:length(MPthr)
        %     nexttile; hold on
        %     plot(MPthr{i}(1:end/2))
        %     ylim([0.5 1])
        %
        % end
        % lgd=legend
        % grid on
        %%
        %     toc
        % figure,plot(mp1{1})
        % figure,plot(MPthr{1})
        % figure,plot(MPthr{2})

        [a,b] = max(MPthr{1});
        loc = MPIthr1{1}(b)+1;

        bar1 = thr{1}(b:b+w-1);
        bar2 = thr{1}(loc:loc+w-1);

        figure,plot(bar1)
        hold on
        plot(bar2)
        zscore(bar1,1)*zscore(bar2,1)'/length(bar1)

        b1 = zscore(bar1,1);
        b2 =zscore(bar2,1);

        % 1-sum((b1-b2).^2)/(2*length(bar1))
        % 1-sum((b1-b2).^2)/(2*length(bar1))

        %%
        [a,b] = findpeaks(MPthr{1}(1:end/2),'MinPeakDistance',200,'SortStr','descend');

        ix = 5;
        pos2 = MPIthr1{1}(b(ix)+1);

        bar1 = thr{1}(b(ix):b(ix)+w-1);
        bar2 =  thr{1}(pos2:pos2+w-1);

        lastPx = find(isnan(thr{1}));

        figure,plot(bar1)
        hold on
        plot(bar2)
        legend({['Starts ' num2str(mod(b(ix),lastPx))],['Starts ' num2str(mod(pos2,lastPx))]},'location','southoutside')
        title(['PCC = ' num2str(a(ix) )])

        pcc = @(x,y) zscore(x,1)*zscore(y,1)'/length(x)

        pcc(bar1,bar2)


    end
    %% both
    %
    % theoryStruct2{1}.filename = 'flip.txt';
    % writematrix(fliplr(importdata( theoryStruct{1}.filename )),theoryStruct2{1}.filename );
    %
    % % idxTheory = 2;
    % k=1;
    % MIN_OVERLAP_PIXELS = 150;
    % numWorkers = 30;
    % out = timestamp;
    % if ispc
    %     command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
    %        theoryStruct2{1}.filename  ' --input_a_file_name=.\'...
    %        theoryStruct{2}.filename ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
    %        ' --output_a_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
    %        ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
    % else
    %
    %        command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
    %         theoryStruct2{1}.filename ' --input_a_file_name=\'...
    %         theoryStruct{2}.filename ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
    %        ' --output_a_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
    %        ' --output_a_index_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
    % end
    % [~,message] = system(command);
    %
    %    mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));
    %
    %     %also get PCC scores.
    % %     tic
    %     % a bit slower import to make sure imports correctly on win
    % fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
    % raw2 = textscan(fid, '%s ');
    % fclose(fid);
    % nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
    % mp1{k} = nan(length(nonanValues),1);
    % mp1{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
    % %
    %
    % figure,plot(mp1{1})

end

