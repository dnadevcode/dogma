function [uniqueNames,uniqueFolds,data] = load_data_wrapper(type,outfold,sets,fold)
    % load any kind of data for HCA!
    %
    %
%     import CBT.Hca.Import.load_data_wrapper;
    
    if nargin < 4
        % default data inputs
        switch type
            case 'consensus'
                sets.fold = '/home/albyback/git/sv/data/mbio_sv/New labels_220kbp all barcodes/';
            case 'kymographs'
                sets.fold = '/home/albyback/git/sv/data/runDataMBIO/';
                addpath(genpath(sets.fold));
            case 'database'
                sets.fold = '/home/albyback/git/sv/data/runDataDatabase/';
            case 'theory' 
                sets.fold = '/home/albyback/git/sv/data/runData/';
                theories = '/home/albyback/git/sv/data/runData/*.fasta';

%                 [tifsName,lambdasName] = create_database_files(source,setsF,resultFold);
            otherwise
        end
    else
        sets.fold = fold;
    end
    
    sets.outfold = outfold;
%           sets.outfold = outfold;

    mkdir(outfold);
    
    import CBT.Hca.Import.load_mbio;
    import CBT.Hca.Import.load_kymograph;
    import CBT.Hca.Import.load_cons; % alternative can calculate consensus here..
    % open file to save txts, could have a different name with a timestamp
    % i.e.
    fd = fopen(fullfile(outfold,'consensus_mat_files.txt'),'w');
    
    switch type
        case 'consensus'
%             sets.fold = '/home/albyback/git/sv/data/mbio_sv/New labels_220kbp all barcodes/';
            
            % read data / would rather load data from txt file..
            listing = dir(strcat(sets.fold,'*.mat') );
            uniqueNames = arrayfun(@(x) x.name, listing,'UniformOutput',false);
            uniqueFolds = arrayfun(@(x) x.folder, listing,'UniformOutput',false);
            for i=1:length(uniqueNames)
                fprintf(fd,'%s\n', fullfile(uniqueFolds{i},uniqueNames{i}));
            end
            [data] = load_cons(uniqueNames,uniqueFolds,outfold);

%             [sets,uniqueNames,uniqueFolds] = load_mbio(sets);
        case 'kymographs'
            % the folders with kymographs have to have a specific structure
            % for us to be able to extract all of them
            [uniqueNames,lambdasName] = create_database_files(sets.fold,sets,outfold);
            uniqueFolds = sets.fold;
            for i=1:length(uniqueNames)
                fprintf(fd,'%s\n', uniqueNames{i});
            end
            [data] = load_kymograph(uniqueNames,lambdasName,sets);
        case 'database'
            [uniqueNames,lambdasName] = create_database_files(sets.fold,sets,outfold);
            uniqueFolds = sets.fold;
            for i=1:length(uniqueNames)
                fprintf(fd,'%s\n', uniqueNames{i});
            end
            [data] = load_kymograph(uniqueNames,lambdasName,sets);
        case 'theory'
            [uniqueNames,lambdasName] = create_database_files(sets.fold,sets,outfold);
            uniqueFolds = sets.fold;
            for i=1:length(uniqueNames)
                fprintf(fd,'%s\n', uniqueNames{i});
            end
            [data] = load_kymograph(uniqueNames,lambdasName,sets);
%             uniqueNames= [uniqueNames 
            % import settings
%             import CBT.Hca.Import.import_hca_settings;
%             [sets2] = import_hca_settings(sets.txt);
            st = length(data)+1;
            [t,theoriesTxt,theoryStruct, sets,theoryGen] = HCA_theory_parallel(theories,data{1}.bpnm,sets.txt2);
           data = [data theoryStruct];
            for i=1:length(theoryStruct)
                data{st}.barcode = theoryGen.theoryBarcodes{i};
                data{st}.fname = theoryStruct{i}.filename;
                data{st}.circ = ~sets.theoryGen.isLinearTF;
                if  data{st}.circ
                    data{st}.bitmask = ones(1,length(data{st}.barcode ));
                else % this case add some zeros
                    data{st}.bitmask = ones(1,length(data{st}.barcode ));
                    data{st}.bitmask(1:5) = 0;
                    data{st}.bitmask(end-4:end) = 0;
                    data{st}.bpnm = data{1}.bpnm;
                end
                st = st+1;
            end
            % sets.resultsDir = 'out/';
            % mkdir(sets.resultsDir);

            % double check puuh_theory.txt
            % generate theory
%             sets.theoryGen.meanBpExt_nm = bpnmTheory;
%             sets.theoryGen.pixelWidth_nm = 159.2;
%             [t,theoriesTxt,theoryStruct, sets,theoryGen] = HCA_theory_parallel(theories,bpnm,'puuh_theory.txt');


        case 'synthetic'
            % here we create synthetic data
            % Options: linear-linear, linear-circular, circular-linear, and
            % circular-circular
            % if we want linear
            
            %generates a listoflinear barcodes
            svTypes = sets.svType*ones(1,sets.numSamples);
            len1 = sets.length1;
            len2 = sets.length2;
            circ = [sets.circ sets.circ];
            sigma = sets.kernelSigma; % same as first theory
            noise = 1-sets.pccScore;
            
            
            import CBT.Hca.Import.gen_synthetic_sv;
            [bar1,bar2,matchTable] = gen_synthetic_sv(len1,len2, svTypes, circ, sigma,noise);
            uniqueNames =[];
            uniqueFolds = [];
            
            data = cell(1,2*length(bar1));
            for i=1:length(bar1);
                %                 % save DATA
                rS = fullfile(outfold,strcat([ num2str(i) '_query_seq.txt']));
                fid = fopen(rS,'w');

                fprintf(fid, '%5.4f ', bar1{i}  );
                fprintf(fid, '\n');
                fclose(fid);
                data{2*i-1}.filename = rS;
                data{2*i-1}.name = rS;

                % save QUERY (this includes the variation)
                rS = fullfile(outfold, strcat([ num2str(i) '_data_seq.txt']));
                fid = fopen(rS,'w');
                if sets.circ % save txt as circular..
                     fprintf(fid, '%5.4f ',bar2{i}   ); % circularity used later..
%                     fprintf(fid, '%5.4f ',[bar2{i} bar2{i}(1:sets.c-1)]   );
                else
                    fprintf(fid, '%5.4f ',bar2{i}   );
                end
                fprintf(fid, '\n');
                fclose(fid);
                data{2*i}.filename = rS;
                data{2*i}.name = rS;


                data{2*i-1}.rawBarcode = bar1{i};
                 data{2*i-1}.length = length( bar1{i});

                data{2*i}.rawBarcode = bar2{i};
                 data{2*i}.length = length( bar2{i});
                data{2*i-1}.matchTable =matchTable{i};
%                 data{2*i}.fname 
                data{2*i-1}.circ = circ(1);
                data{2*i}.circ = circ(2);
                data{2*i-1}.isLinearTF = ~circ(1);
                data{2*i}.isLinearTF = ~circ(2);
                data{2*i}.rawBitmask = ones(1,length(data{2*i}.rawBarcode ));
                data{2*i-1}.rawBitmask = ones(1,length(data{2*i-1}.rawBarcode ));
                data{2*i-1}.bpnm = 0.2; 
                data{2*i}.bpnm = 0.2;

            end
% %             
% %             import Rand.generate_linear_sv;
% %             [bar1,bar2,matchTable,lengths]  = arrayfun(@(x) generate_linear_sv(len1, lenVar, x,sets),sets.svList,'UniformOutput',false);
%             
% %             idx=5;
% %             res.bar1 = bar1{idx};
% %             res.bar2 = bar2{idx};
% %             res.matchTable = matchTable{idx};
% %             res.pass = ones(1,size(matchTable{idx},1));
% % 
% %             sets.fold = outfold;
% % %             plot full result
% %             import Plot.plot_sv_full;
% %             cellfun(@(x) plot_sv_full( x,sets,'test111.eps',2,0,0.5),{res},'UniformOutput',false);
% 
%             %
%             % test - all pccs should be equal one
%             
%             % these can be saved as .tifs
%             
%             % or
% %             import Rand.generate_one_type_sv;
% %             [randStruct] = generate_one_type_sv(sets,vals(i));
        case "random"
        %generates a listoflinear barcodes
            len1 = sets.length1;

%             len2 = sets.length2;
%             svTypes = zeros(1,sets.numSamples/2);

%             import CBT.Hca.Core.Pvalue.gensv;
%             [bar1,bar2,~,~]  = arrayfun(@(x) gensv(len1,len2, x,sets.kernelsigma,sets.islinear),svList,'UniformOutput',false);

            for i=1:sets.pvalue.numRnd;
                bar1 = normrnd(0,1,1,len1);

                if ~sets.circ
                    bar1 =  imgaussfilt(bar1,sets.kernelSigma);
                else
                    bar1T = [bar1 bar1(1:10)];
                    bar1T  =  imgaussfilt(bar1T,sets.kernelSigma);
                    bar1 = bar1T(1:length(bar1));
                end
                
                data{i}.rawBarcode = bar1;
                data{i}.rawBitmask = ones(1,length(data{i}.rawBarcode ));
                data{i}.isLinearTF = ~sets.circ;

    
            end
                uniqueNames=[];uniqueFolds=[];
%             matchTable = [];
%            lengths  = [len1 lenVar];
%        
            

        otherwise
            
    end
    
    fclose(fd);
    
    
    
% [file,fold] = uigetfile('*.mat','select query file','MultiSelect','on', sets.fold);

% listing = dir(strcat(sets.fold,'*.mat') );
% names = arrayfun(@(x) x.name, listing,'UniformOutput',false);
% folds = listing(1).folder;
% 
% sets.fold = '/home/albyback/git/sv/data/mbio_sv/New labels_220kbp all barcodes/';
% % [file,fold] = uigetfile('*.mat','select query file','MultiSelect','on', sets.fold);
% 
% listing = dir(strcat(sets.fold,'*.mat') );
% names = arrayfun(@(x) x.name, listing,'UniformOutput',false);
% folds = listing(1).folder;
% 
% 
% 
%     for i=1:length(fileNames)
%         
%         
%     end

end

