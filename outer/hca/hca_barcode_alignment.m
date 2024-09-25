function [] = hca_barcode_alignment(useGUI, hcaSets)
    %   Replicates HCA_Gui, nicer graphical user interface
    %   Written by Albertas Dvirnas

    % GUI globals
    duplicates = []; % to store duplicate settings between windows
    shrinksorter = []; % to store shrinksorter
    hcaaligner = [];
    hcatheory = [];

    % Other globals:
    versionHCA = [];
    iL = [];
    t1 =[]; t2 = []; t3 = []; t4 = [];
    tsAlignmentVisual =[];

    if nargin < 2
        import CBT.Hca.Import.import_hca_settings;
        [hcaSets] = import_hca_settings('hca_settings.txt');
    end
        
    if nargin >=1
        hcaSets.useGUI = useGUI;  
    end

    if hcaSets.useGUI
        % generate menu
        % https://se.mathworks.com/help/matlab/ref/uimenu.html
%         [hcaSets,tsHCC,textList,textListT,itemsList,...
%                 hHomeScreen,hPanelRawKymos,hPanelAlignedKymos,hPanelTimeAverages,hAdditional,tshAdd]= generate_gui();
        [hFig,hPanel,h] = generate_gui();     
    else

    end
 
    function [hFig,hPanel,h] = generate_gui()

            mFilePath = mfilename('fullpath');
            mfolders = split(mFilePath, {'\', '/'});
            versionHCA = importdata(fullfile(mfolders{1:end-1},'VERSION'));

            hFig = figure('Name', ['HCA Barcode Alignment v' versionHCA{1}], ...
                'Units', 'normalized', ...
                'InnerPosition', [0.05 0.1 0.8 0.8], ...
                'NumberTitle', 'off',...
                'HandleVisibility', 'on', ...
                'MenuBar', 'none',...
                'ToolBar', 'none' ...
            );

            % Fancier UI
            %             hFig = uifigure('Name', ['HCA Barcode Alignment v' versionHCA{1}], ...
            %                 'Units', 'normalized', ...
            %                 'InnerPosition', [0.05 0.1 0.8 0.8], ...
            %                 'NumberTitle', 'off',...
            %                 'HandleVisibility', 'on' ...
            %             );
            m = uimenu(hFig,'Text','HCA');
            cells1 = {'HCA v4.7','Theory','Alignment'};
            mSub = cellfun(@(x) uimenu(m,'Text',x),cells1,'un',false);
            mSub{1}.MenuSelectedFcn = @SelectedOldHCA;

            cellsTheory = {'Generate theory','Concatenate theories'};

            mSubTheory  = cellfun(@(x) uimenu(mSub{2},'Text',x),cellsTheory,'un',false);

            mSubTheory{1}.MenuSelectedFcn = @SelectedHCAtheory;
            mSubTheory{2}.MenuSelectedFcn = @SelectedConcatenateTheories;

            cellsAlignment = {'Run Alignment','Shrink sorter','Duplicates sorter','Save Alignment Result','Load session'};
            mSubALignment  = cellfun(@(x) uimenu(mSub{3},'Text',x),cellsAlignment,'un',false);
            mSubALignment{1}.MenuSelectedFcn = @SelectedHCAAlignment;
            mSubALignment{2}.MenuSelectedFcn = @SelectedShrinkSorter;
            mSubALignment{3}.MenuSelectedFcn = @SelectedDuplicatesSorter;
            mSubALignment{4}.MenuSelectedFcn = @SelectedSaveAlignmentSession;
            mSubALignment{5}.MenuSelectedFcn = @SelectedLoadAlignmentSession;


            hPanel = uipanel('Parent', hFig,  'Units', 'normalized','Position', [0 0 1 1]);
            h = uitabgroup('Parent',hPanel, 'Units', 'normalized', 'Position',   [0 0 1 1]);
    end

    function [structFiles] = get_files_function(tsSet,structFiles,run_handle)
        % get_files_function, create a basic UI element with inport,
        % settings, and run button
        % v5.2.0

        structSets = structFiles.sets;
        structNames = structFiles.names;

        fnames = fieldnames(structSets.default);
        tISets = ones(1,length(fnames));
        if isfield(structSets,'clIdx')
            tISets(structSets.clIdx) = 0;
        end

        if isfield(structSets,'testName')
            testName = structSets.testName;
        else
            testName = 'kymo_example.tif';
        end
    
        if isfield(structSets,'fileext')
            fileext = structSets.fileext;
        else
            fileext = '*.tif';
        end

        structFiles = CreateSelectWindow(structFiles,tsSet, testName,fileext,0.95);
        structFiles.runButton = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{'Run'},'Callback',run_handle,'Units', 'normal', 'Position', [0.7 0.1 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]

        if isfield(structSets, 'showtheoryimport')
            structFiles.theory = CreateSelectWindow(structFiles,tsSet, '','.mat',0.85, 'Theories to compare against (.mat/.txt)');
        end

        checkListIdx = find(~tISets);
        itemListIdx = find(tISets);

        iL = cell(1,length(checkListIdx));
        for i = 1:length(checkListIdx)
            iL{i} = uicontrol('Parent', tsSet, 'Style', 'checkbox','Value', structSets.default.(fnames{checkListIdx(i)}),'String',structNames{checkListIdx(i)},'Units', 'normal', 'Position', [0.45 .78-0.05*i 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        end
                
        positionsText = cell(1,length(itemListIdx));
        positionsBox = cell(1,length(itemListIdx));

        for i=1:length(itemListIdx) % these will be in two columns
            positionsText{i} =   [0.2-0.2*mod(i,2) .83-0.1*ceil(i/2) 0.2 0.03];
            positionsBox{i} =   [0.2-0.2*mod(i,2) .78-0.1*ceil(i/2) 0.15 0.05];
        end

        tL = cell(1,length(itemListIdx));

        for i=1:length(itemListIdx)
            tL{i} = uicontrol('Parent', tsSet, 'Style', 'text','String',structNames{itemListIdx(i)},'Units', 'normal', 'Position', positionsText{i},'HorizontalAlignment','Left');
            tL{i} = uicontrol('Parent', tsSet, 'Style', 'edit','String',{num2str(structSets.default.(fnames{itemListIdx(i)}))},'Units', 'normal', 'Position', positionsBox{i});
        end

        % save in duplicates
        structFiles.Item = iL;
        structFiles.Text = tL;
    end

    function structFiles = CreateSelectWindow(structFiles,tsSet, testName, fileext, posY, txt)
        if nargin < 6
           % KYMO
            txt = 'Kymos to import (.tif/.mat)';
        end
        structFiles.dotText = uicontrol('Parent', tsSet, 'Style', 'text','String',{txt},'Units', 'normal', 'Position', [0.2 posY 0.2 0.03],'HorizontalAlignment','Left');
        if ~isempty(testName)
            testName = fullfile(fileparts(mfilename('fullpath')),'files',testName);
        end
        structFiles.dotImport = uicontrol('Parent', tsSet, 'Style', 'edit','String',{testName},'Units', 'normal', 'Position', [0 posY-0.05 0.5 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        set(structFiles.dotImport, 'Min', 0, 'Max', 25);% limit to 10 files via gui;
        structFiles.dotButton = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{'Browse folder'},'Callback',@(src, event) selection_folder(structFiles.dotImport,event,fileext),'Units', 'normal', 'Position', [0.6 posY-0.05 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        structFiles.dotButtonFile = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{'Browse file'},'Callback',@(src, event) selection_file(structFiles.dotImport,event),'Units', 'normal', 'Position', [0.7 posY-0.05 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        structFiles.dotButtonUI = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{'uigetfiles'},'Callback',@(src, event) selection_uipickfiles(structFiles.dotImport,event),'Units', 'normal', 'Position', [0.8 posY-0.05 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]

    end

    %% Selected functions
    function SelectedOldHCA(~,~)
            disp('Running HCA_Gui v. < 5.5.0')
            delete(hFig);
            HCA_Gui;
    end

    function SelectedDuplicatesSorter(~,~) % duplicates sorter
        disp('Running Duplicates sorting')
        if hcaSets.useGUI
            if isempty(t4)
                t4 = uitab(h, 'title', 'Duplicates sorter');
            end
            h.SelectedTab = t4; 
        end
        tsDuplicates = uitabgroup('Parent',t4, 'Units', 'normalized','Position', [0.01 0.01 0.99 0.99]);
        tsDuplicatesSettings = uitab(tsDuplicates, 'title', 'Duplicate sorter settings');
        [duplicates.sets,duplicates.names] = Core.Default.read_default_sets('duplicatessets.txt');

        [duplicates] = get_files_function(tsDuplicatesSettings,duplicates, @run_duplicatessorter);
    end


    function SelectedShrinkSorter(~,~) % shrink sorter
        disp('Running shrink sorting')
        if hcaSets.useGUI
            if isempty(t3)
                t3 = uitab(h, 'title', 'Shrink sorter');
            end
            h.SelectedTab = t3; 
        end
        tsShrink = uitabgroup('Parent',t3, 'Units', 'normalized','Position', [0.01 0.01 0.99 0.99]);
        tsShrinkSettings = uitab(tsShrink, 'title', 'Shrink sorter settings');
        % get default settings
        [shrinksorter.sets,shrinksorter.names] = Core.Default.read_default_sets('shrinksortersets.txt');
        % use get_files_funcion
        [shrinksorter] = get_files_function(tsShrinkSettings, shrinksorter, @run_shrinksorter);
    end

    function SelectedLoadAlignmentSession(~,~)
       [sessionfile,path] = uigetfile('*.mat','Select session result file to open');
       load(fullfile(path,sessionfile),'theoryStruct','rezMax','hcaSets','barcodeGenC');

        if isempty(t1)
            t1 = uitab(h, 'title', 'HCA alignment');
        end
        h.SelectedTab = t1; 
        
        tsAlignment = uitabgroup('Parent',t1, 'Units', 'normalized','Position', [0.01 0.01 0.99 0.99]);
        tsAlignmentVisual = uitab(tsAlignment, 'title', 'Visual results');
        

        import CBT.Hca.Core.Comparison.combine_theory_results;
        [comparisonStruct] = combine_theory_results(theoryStruct, rezMax);
   

        import Core.run_visual_fun;
        if ~isempty(comparisonStruct)
            run_visual_fun(barcodeGenC,[], comparisonStruct, theoryStruct, hcaSets, tsAlignmentVisual);  
        end    

    end
    function SelectedSaveAlignmentSession(~,~)
        disp('Saving HCA_alignment session')

        timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

        saveFold = fileparts(hcaSets.kymofolder{1});


        saveName = [saveFold,'session_',timestamp,'.mat'];
        if exist('rezMax','var')
            save(saveName,'barcodeGenC','hcaSets','rezMax','theoryStruct');
            disp(['Session saved at ', saveName])
        else
            disp('Results are not saved. They need to be in the workspace in order to be saved')
        end
    end

    function SelectedHCAAlignment(~,~)
            disp('Running HCA_alignment')

            if hcaSets.useGUI
                if isempty(t1)
                    t1 = uitab(h, 'title', 'HCA alignment');
                end
                h.SelectedTab = t1; 

                tsAlignment = uitabgroup('Parent',t1, 'Units', 'normalized','Position', [0.01 0.01 0.99 0.99]);
                tsAlignmentSettings = uitab(tsAlignment, 'title', 'Alignment settings');
                tsAlignmentVisual = uitab(tsAlignment, 'title', 'Visual results');

                [hcaaligner.sets,hcaaligner.names] = Core.Default.read_default_sets('hcaalignmentsets.txt');
                [hcaaligner] = get_files_function(tsAlignmentSettings,hcaaligner, @run_alignment);
            end

    end

    function SelectedHCAtheory(~,~)
            disp('Running HCA_theory')

            if hcaSets.useGUI
                if isempty(t2)
                    t2 = uitab(h, 'title', 'HCA theory');
                end
                h.SelectedTab = t2; 

                tsTheory = uitabgroup('Parent',t2, 'Units', 'normalized','Position', [0.01 0.01 0.99 0.99]);
                tsTheorySettings = uitab(tsTheory, 'title', 'Theory settings');

                [hcatheory.sets,hcatheory.names] = Core.Default.read_default_sets('hcasets.txt');
                [hcatheory] = get_files_function(tsTheorySettings,hcatheory, @run_theory);
            end
    end

    function SelectedConcatenateTheories(~,~)
        disp('Running Concatenate theories')
        [FILENAME] = uipickfiles;
        data1 = cell(1,length(FILENAME));
        for i=1:length(FILENAME)
            data1{i} = load(FILENAME{i});
        end

        conc = 0;

        theoryGen = data1{1}.theoryGen;
        for i=2:length(FILENAME)
            if isequal(data1{1}.theoryGen.sets,data1{i}.theoryGen.sets)
                theoryGen.theoryBarcodes = [theoryGen.theoryBarcodes data1{i}.theoryGen.theoryBarcodes];
                theoryGen.theoryBitmasks = [theoryGen.theoryBitmasks data1{i}.theoryGen.theoryBitmasks];
                theoryGen.theoryNames = [theoryGen.theoryNames data1{i}.theoryGen.theoryNames];
                theoryGen.theoryIdx = [theoryGen.theoryIdx data1{i}.theoryGen.theoryIdx];
                theoryGen.bpNm = [theoryGen.bpNm data1{i}.theoryGen.bpNm];               
                conc = 1;       
            else
                disp(['Settings are not the same in the two files, cant merge ', num2str(i), ' dataset']);
            end
        end
        if conc
        [fileN,fileL] = uiputfile(strrep(FILENAME{1},'.mat','_concatenated.mat'));
        save(fullfile(fileL,fileN),'theoryGen');
        disp(['Concatenated to ', fullfile(fileL,fileN)]);
        end
    end

%% run functions
    function run_duplicatessorter(~,~)
            display(['Started shrink sorter hca_alignment v',versionHCA{1}]);

            duplicatesSets = save_settings_hca(duplicates.sets, duplicates.Item,duplicates.Text);

            duplicatesSets.kymofolder = duplicates.dotImport.String;

            import Core.duplicatessorter;
            [duplicateInfo,oS] = duplicatessorter(duplicatesSets);

            % save duplicate results info

    end


    function run_shrinksorter(~,~)
            display(['Started shrink sorter hca_alignment v',versionHCA{1}]);

            sorterSets = save_settings_hca(shrinksorter.sets, shrinksorter.Item,shrinksorter.Text);

            sorterSets.kymofolder = shrinksorter.dotImport.String;

            import Core.shrink_finder_fun;
            shrink_finder_fun(sorterSets);

    end

    function run_alignment(~,~)
            
            display(['Started analysis hca_alignment v',versionHCA{1}])
            
            hcaSetsDefault.default = hcaSets;

            hcaSetsDefault = save_settings_hca(hcaaligner.sets, hcaaligner.Item,hcaaligner.Text,hcaSetsDefault);

                       % incorporate directly to list!
            if  hcaSetsDefault.default.comparisonMethod==1
                 hcaSetsDefault.default.comparisonMethod  = 'mass_pcc';
            else
                 hcaSetsDefault.default.comparisonMethod  = 'mpnan';
            end
            hcaSets = hcaSetsDefault.default;
            hcaSets.kymofolder = hcaaligner.dotImport.String;

            hcaSets.theoryfolder = hcaaligner.theory.dotImport.String;

            % run alignment            
            import Core.run_hca_alignment;
            [barcodeGenC,consensusStruct, comparisonStruct, theoryStruct, hcaSets] = run_hca_alignment(hcaSets);

            import Core.run_visual_fun;
            if ~isempty(comparisonStruct)

                run_visual_fun(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct, hcaSets, tsAlignmentVisual);

            end           
    end

        function run_theory(~,~)
            
            display(['Started analysis hca_theory v',versionHCA{1}])

            hcatheorySets = save_settings_hca(hcatheory.sets, hcatheory.Item,hcatheory.Text);
            hcatheorySets = hcatheorySets.default;

            hcatheorySets.folder = hcatheory.dotImport.String;

            disp(['N = ',num2str(length(  hcatheorySets.folder )), ' sequences to run'])

            import Core.run_hca_theory;
            run_hca_theory(hcatheorySets);


        end

    % select tifs in folder
    function selection_folder(src, ~, fileext)
        [rawNames] = uigetdir(pwd,strcat('Select folder with file(s) to process'));
        rawFiles = [dir(fullfile(rawNames,fileext))];
        rawNames = arrayfun(@(x) fullfile(rawFiles(x).folder,rawFiles(x).name),1:length(rawFiles),'UniformOutput',false);
        set(src, 'String', rawNames);
    end   

    function selection_file(src, ~) % {'*.tif';'*.mat';}
        [FILENAME, PATHNAME] = uigetfile(fullfile(pwd,'*.*'),strcat(['Select file(s) to process']),'MultiSelect','on');
        name = fullfile(PATHNAME,FILENAME);
        if ~iscell(name)
            name  = {name};
        end
        set(src, 'String', name);
    end  

    function selection_uipickfiles(src, ~)
        [FILENAME] = uipickfiles;
        name = FILENAME';
        if ~iscell(name)
            name  = {name};
        end
        set(src, 'String', name);

    end

    % save settings for particular items/checklist, works for any settings
    function hcaSetsDefault = save_settings_hca(hcaSets,iL,tL,hcaSetsDefault)

        if nargin <4 
            hcaSetsDefault = struct();
        end

        fnames = fieldnames(hcaSets.default);
        tISets = ones(1,length(fnames));
        if isfield(hcaSets,'clIdx')
            tISets(hcaSets.clIdx) = 0;
        end
    
        checkListIdx = find(~tISets);
        itemListIdx = find(tISets);

        for i = 1:length(checkListIdx)
            hcaSetsDefault.default.(fnames{checkListIdx(i)}) = iL{i}.Value;
        end

        for i = 1:length(itemListIdx)
            if ~isnan(str2double(tL{i}.String))
                hcaSetsDefault.default.(fnames{itemListIdx(i)}) = str2double(tL{i}.String);
            else
                hcaSetsDefault.default.(fnames{itemListIdx(i)}) = tL{i}.String{1};
            end
        end   

    end
    
end