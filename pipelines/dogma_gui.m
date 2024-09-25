function [] = dogma_gui(useGUI, bargiSets)
    %   Replicates bargi_run, nicer graphical user interface
    %   Written by Albertas Dvirnas

    % globals
     bargiStruct = [];

    % load default settings
    import Core.Default.read_default_bargi_sets; % loads default bargi_settings from file
    if nargin < 2
        [bargiSets,bargiNames] = Core.Default.read_default_bargi_sets('bargi_settings.txt');
    end
        
    if nargin >=1
        bargiSets.useGUI = useGUI;  
    end

    if bargiSets.useGUI
        % generate menu
        % https://se.mathworks.com/help/matlab/ref/uimenu.html
        [hFig,hPanel,h,ts,tsV] = generate_gui();     

    else
% 
% 
    end
%  

    function [hFig,hPanel,h,ts,tsV] = generate_gui()
            % generate graphical user interface (GUI)

            mFilePath = fileparts(fileparts(mfilename('fullpath')));
            versionBargi = importdata(fullfile(mFilePath,'VERSION'));
            bargiSets.version =  versionBargi{1};

            % todo: change to uifigure
            hFig = figure('Name', ['Bardenas ' versionBargi{1}], ...
                'Units', 'normalized', ...
                'InnerPosition', [0.05 0.1 0.8 0.8], ...
                'NumberTitle', 'off',...
                'HandleVisibility', 'on',...
                'MenuBar', 'none',...
                'ToolBar', 'none');

            m = uimenu(hFig,'Text','Bardenas');
            cells1 = {'Run Bardenas','Load Results'};

            mSub = cellfun(@(x) uimenu(m,'Text',x),cells1,'un',false);
            mSub{1}.MenuSelectedFcn = @bargi_create_gui;
            mSub{2}.MenuSelectedFcn = @bargi_load_results;

            hPanel = uipanel('Parent', hFig,  'Units', 'normalized','Position', [0 0 1 1]);
            h = uitabgroup('Parent',hPanel,  'Units', 'normalized','Position', [0 0 1 1]);

            ts = uitab(h, 'title', 'Bardenas settings');
            tsVisual = uitab(h, 'title', 'Visual results');
            tsVs = uitabgroup('Parent',tsVisual, 'Units', 'normalized','Position', [0.01 0.01 0.99 0.99]);
            tsV.single = uitab(tsVs, 'title', 'Single');
            tsV.islands = uitab(tsVs, 'title', 'Islands');
            tsV.graph = uitab(tsVs, 'title', 'Graph');
            tsV.block = uitab(tsVs, 'title', 'Block');


    end

    function bargi_load_results(~,~)
        sessionFileLoc = uigetfile(pwd);


    end



    function bargi_create_gui(~,~)

        disp('Running  bargi')
        if bargiSets.useGUI
            if isempty(ts)
                ts = uitab(h, 'title', 'Bargi');
            end
            h.SelectedTab = ts; 
        end
%         tsBargi = uitabgroup('Parent',ts);
        tempsets.sets = bargiSets;
        tempsets.names = bargiNames;

%         tsBargiSettings = uitab(tsBargi, 'title', 'Bargi settings');
        [bargiStruct] = get_files_function(ts,tempsets, @run_bargi);
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

        structFiles.dotImport = uicontrol('Parent', tsSet, 'Style', 'edit','String',{fullfile(fileparts(mfilename('fullpath')),'files',testName)},'Units', 'normal', 'Position', [0 0.9 0.5 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        set(structFiles.dotImport, 'Min', 0, 'Max', 25);% limit to 10 files via gui;
        structFiles.dotButton = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{'Browse folder'},'Callback',@(src, event) selection_folder(structFiles.dotImport,event,fileext),'Units', 'normal', 'Position', [0.6 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        structFiles.dotButtonFile = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{'Browse file'},'Callback',@(src, event) selection_file(structFiles.dotImport,event),'Units', 'normal', 'Position', [0.7 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        structFiles.dotButtonUI = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{'uigetfiles'},'Callback',@(src, event) selection_uipickfiles(structFiles.dotImport,event),'Units', 'normal', 'Position', [0.8 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        structFiles.runButton = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{'Run'},'Callback',run_handle,'Units', 'normal', 'Position', [0.7 0.2 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
%         structFiles.runButtonIslands = uicontrol('Parent', tsSet, 'Style', 'pushbutton','String',{''},'Callback',run_handle,'Units', 'normal', 'Position', [0.7 0.2 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]

        checkListIdx = find(~tISets);
        itemListIdx = find(tISets);

        iL = cell(1,length(checkListIdx));
        for i = 1:length(checkListIdx)
            iL{i} = uicontrol('Parent', tsSet, 'Style', 'checkbox','Value', structSets.default.(fnames{checkListIdx(i)}),'String',structNames{checkListIdx(i)},'Units', 'normal', 'Position', [0.45 .83-0.05*i 0.3 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        end
                
        positionsText = cell(1,length(itemListIdx));
        positionsBox = cell(1,length(itemListIdx));

        for i=1:length(itemListIdx) % these will be in two columns
            positionsText{i} =   [0.2-0.2*mod(i,2) .88-0.1*ceil(i/2) 0.2 0.03];
            positionsBox{i} =   [0.2-0.2*mod(i,2) .83-0.1*ceil(i/2) 0.15 0.05];
        end

        tL = cell(1,length(itemListIdx));

        for i=1:length(itemListIdx)
            tL{i} = uicontrol('Parent', tsSet, 'Style', 'text','String',structNames{itemListIdx(i)},'Units', 'normal', 'Position', positionsText{i},'HorizontalAlignment','Left');
            tL{i} = uicontrol('Parent', tsSet, 'Style', 'edit','String',{num2str(structSets.default.(fnames{itemListIdx(i)}))},'Units', 'normal', 'Position', positionsBox{i});
        end

        % save
        structFiles.Item = iL;
        structFiles.Text = tL;
    end


    function run_bargi(src, event)

%             display(['Started analysis bargi v',versionBargi])
%             bargiSets.kymofolder = dotImport.String;


            bargiStructSets = save_settings_bargi(bargiStruct.sets, bargiStruct.Item,bargiStruct.Text);
            bargiStructSets = bargiStructSets.default;

            bargiStructSets.kymofolder = bargiStruct.dotImport.String;

            disp(['N = ',num2str(length(  bargiStructSets.kymofolder )), ' sequences to run'])

            dogma_run(bargiStructSets,[],[],tsV);
    end


    % save settings for particular items/checklist, works for any settings
    function hcaSetsDefault = save_settings_bargi(hcaSets,iL,tL,hcaSetsDefault)

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
        disp(['Selected ', num2str(length(name)) , ' files']);
    end  

    function selection_uipickfiles(src, ~)
        [FILENAME] = uipickfiles;
        name = FILENAME';
        if ~iscell(name)
            name  = {name};
        end
        set(src, 'String', name);

    end




end