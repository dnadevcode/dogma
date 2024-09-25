function [bgAll,bG,kymoStructs,synthStr, synthStr2,theoryStruct] = load_chrom_data(fold, method, depth, timeframes)
    %   load_chrom_data
    %
    %   Args:
    %       fold - folder with data
    %       method - kymo extraction method
    %       depth - folder depth
    %       timeframes - number of timeframes
    %
    %   Returns:
    %       bgAll - cell with all barcodes
    %       bG - cell with cell of barcoces from each data
    %       kymoStructs - kymograph data

    % Example:
%     userDir = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/ecoli_2_new/12_20_2/';
%     [bgAll,bG,kymoStructs] = load_chrom_data(userDir, 'tiff', 1, 20);

    if nargin < 2 || isempty(method)
        method = 1;
    end

    if nargin < 3 || isempty(depth)
        depth = 0;
    end

    if nargin < 4 % number of time-frames
        timeframes = 20;
    end  

    import Core.gen_masks_simple;
    import Core.get_all_kymos_from_dir;
    import Core.kymo_sig_std;
    import Plot.plot_kymo_masks;

    import Core.Default.read_default_sets;
    import Core.shrink_finder_fun; % needs HCA

    synthStr = [];
    synthStr2 = [];
    theoryStruct = [];

    switch method
        case 'synth'
            import Core.load_synth_data;
            [bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(50,2.72,850,100,7,0.05,1,10000,150);
            kymoStructs = [];
            bgAll = horzcat(bG{:})';

        case 'kymo'
            % step1: get all kymos
            [kymoStructs] = get_all_kymos_from_dir(fold, timeframes, depth);
            
            if ~isfield(kymoStructs{1}{1},'unalignedBitmask')
                [kymoStructs] = gen_masks_simple(kymoStructs); % gen masks (using gaussian mixture)
            end
%               plot_kymo_masks(kymoStructs,1)

            [kymoStructs] = kymo_sig_std(kymoStructs); % gen sigstd

            % Step2: shrink sorter
            [hcaSets] = read_default_sets('shrinksortersets.txt');
            kymoStructsUpdated = cell(1,length(kymoStructs));
            kymoKeep = cell(1,length(kymoStructs));

            for i=1:length(kymoStructs)
                [kymoStructsUpdated{i},kymoKeep{i}] = shrink_finder_fun(hcaSets, kymoStructs{i},0);
%                 goodKymos = find(sum(isnan(kymoKeep{i})')~=2);
            end
%             badKymos = find(sum(isnan(kymoKeep)')==2);
%               plot_kymo_masks(kymoStructsUpdated,1)
%             import Featuretrack.path_track;
%             path_track(kymoStructsUpdated{1}{2}.shiftalignedKymo,kymoStructsUpdated{1}{2}.feature_paths)

            %
            % align (all timeframes
            [kymoStructs] = gen_alingned_data(kymoStructsUpdated);
        
            % extract barcodes
            [bG] = extract_barcodes(kymoStructs);
    
            bgAll = horzcat(bG{:});

            % kymoParams.nmpx = 
% kymoParams.nmbp =

        case 'tiff'
            
            % settings file 
            import DBM4.UI.find_default_settings_path;
            defaultSettingsFilepath = find_default_settings_path('DBMnew.ini');
            import Fancy.IO.ini2struct;
            dbmOSW.DBMSettingsstruct = ini2struct(defaultSettingsFilepath);
            dbmOSW.DBMSettingsstruct.minLen = 150;
            dbmOSW.DBMSettingsstruct.minOverlap = 150;
            
            if depth==1
                
                d = dir(fold);
                
                dfolders = d([d(:).isdir]);
                dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
            else
                dfolders(1).folder = fold;
                dfolders(1).name = '';
            end
            
            allKymos = cell(1,length(dfolders));
            bG = cell(1,length(dfolders));
            kymoStructs = cell(1,length(dfolders));
            
            for idFold = 1:length(dfolders)    
                foldT = fullfile(dfolders(idFold).folder,dfolders(idFold).name);
            
                import DBM4.GenomAs.run_genome_assembly_pipeline;
                [~,bG{idFold},kymoStructs{idFold}] = run_genome_assembly_pipeline(foldT, dbmOSW);
            % bG = [];
            end
            bgAll = horzcat(bG(:));
            
            
            kymoStructs = cell(1,length(length(dfolders)));
            for idFold = 1:length(dfolders)  
                try
            
                spltName = strsplit(dfolders(idFold).name ,'_');
                spltName2 = strsplit(spltName{end},'nm');
                spltName3 = strsplit(spltName{end-1},'nm');
                nmBpidFold= str2double(spltName2{1});
                nmpxidFold = str2double(spltName3{1});
            
                    kymoStructs{idFold}{1}.idFold = idFold;
                    kymoStructs{idFold}{1}.nmBpidFold = nmBpidFold;
                    kymoStructs{idFold}{1}.nmpxidFold = nmpxidFold;
            
                catch
                end
            end

        otherwise 

    end
 
% figure,tiledlayout(ceil(sqrt(length(kymoStructs))),ceil(length(kymoStructs)/sqrt(length(kymoStructs))),'TileSpacing','none','Padding','none')
% for i=1:length(kymoStructs)
% %         nexttile;        imagesc(imresize(kymoStructs{ii}{i}.unalignedKymo,[200 500]));    title(num2str(i));
%         nexttile;        imshowpair(imresize(kymoStructs{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  );    title(num2str(i));
% end


end

