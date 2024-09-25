function [kymoStructs] = get_all_kymos_from_dir(fold, timeframes, depth)

if nargin < 2
    timeframes = inf;
end

if nargin < 3
    depth = 0;
end


if ~iscell(fold)
    d = dir(fold);
    dfolders = d([d(:).isdir]);
    dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));

    allKymos = cell(1,length(dfolders));
    allMat = cell(1,length(dfolders));
    
    
    
    
    % foldName = cell(1,length(dfolders));
    for idFold = 1:length(dfolders)   
    
        if depth == 1
            subDir = dir(fullfile(dfolders(idFold).folder,dfolders(idFold).name));
            subDir(ismember( {subDir.name}, {'.', '..'})) = [];  %remove . and ..
            subDir = subDir(1); % first. todo: make sure this is folder...
        else
            subDir = dfolders(idFold);
        end
    
    
        files = dir(fullfile(subDir.folder,subDir.name,'*.tif'));
        info.foldName{idFold} = subDir.name;
        filesC = arrayfun(@(x) fullfile(files(x).folder,files(x).name),1:length(files),'un',false);
        allKymos{idFold} = filesC;
    
        if isempty(files)
            files = dir(fullfile(subDir.folder,subDir.name,'*.mat'));
            filesC = arrayfun(@(x) fullfile(files(x).folder,files(x).name),1:length(files),'un',false);
            allMat{idFold} = filesC;
        end
    end

else

    [~,~,ed] = fileparts(fold{1});
    if isequal(ed,'.mat')
        allMat{1} = fold;
        allKymos{1}  = [];
    else

        allKymos{1} = fold;
        allMat{1} = [];
    end

end
%% now run through each

kymoStructs = cell(1,length(allKymos));
for idFold = 1:length(allKymos)  
    try
    spltName = strsplit(  info.foldName{idFold}  ,'_');
    spltName2 = strsplit(spltName{end},'nm');
    spltName3 = strsplit(spltName{end-1},'nm');
    nmBpidFold= str2double(spltName2{1});
    nmpxidFold = str2double(spltName3{1});
    catch
        nmBpidFold = nan;
        nmpxidFold = nan;
    end
    for idKym=1:length(allKymos{idFold})  % todo: should get directly from dbmstruct
         kymoStructs{idFold}{idKym}.name = allKymos{idFold}{idKym};
        % if old structure
        if length(imfinfo(allKymos{idFold}{idKym}))==1
            kymoStructs{idFold}{idKym}.unalignedKymo = imread( allKymos{idFold}{idKym},1);

        else
            kymoStructs{idFold}{idKym}.unalignedKymo = imread( allKymos{idFold}{idKym},2);
            kymoStructs{idFold}{idKym}.unalignedBitmask = imread( allKymos{idFold}{idKym},3);
            kymoStructs{idFold}{idKym}.unalignedBitmask = kymoStructs{idFold}{idKym}.unalignedBitmask(1:min(end,timeframes),:);

        end
        kymoStructs{idFold}{idKym}.unalignedKymo = kymoStructs{idFold}{idKym}.unalignedKymo(1:min(end,timeframes),:);
        kymoStructs{idFold}{idKym}.idFold = idFold;
        kymoStructs{idFold}{idKym}.nmBpidFold = nmBpidFold;
        kymoStructs{idFold}{idKym}.nmpxidFold = nmpxidFold;
    end

    if ~isempty(allMat{idFold})
        data = load(allMat{idFold}{1}); % currently load just first one
        kymoStructs{idFold} = data.kymoStructs;
        for k=1:length( kymoStructs{idFold})
            kymoStructs{idFold}{k}.nmBpidFold = nmBpidFold;
            kymoStructs{idFold}{k}.nmpxidFold = nmpxidFold;
        end
    end
end
        

end

