function [sets,uniqueNames,uniqueFolds] = load_mbio(sets)


    % read data / would rather load data from txt file..
    listing = dir(strcat(sets.fold,'*.mat') );
    sets.names = arrayfun(@(x) x.name, listing,'UniformOutput',false);
    sets.folds = arrayfun(@(x) x.folder, listing,'UniformOutput',false);

    % split info different experiments by name
    classes = cell(1,length(sets.names));
    indexes = zeros(1,length(sets.names));

    Index1 = cellfun(@(x) strfind(x, 'P'), sets.names ); 
    Index2 = cellfun(@(x) strfind(x, 'K'), sets.names ); 
    Index3 = cellfun(@(x) strfind(x, '.'), sets.names ); 
    % 
    for i=1:length(sets.names);
        classes{i} =sets.names{i}(Index1(i)+1:Index2(i)-1);
        indexes(i) =  str2num(sets.names{i}(Index2(i)+1:Index3(i)-1));
    end

    % extract unique classes;
     [uc, firstPos, idc] = unique( classes ) ;
     counts = accumarray( idc, ones(size(idc)) ) ;

    uniqueClasses = cell(1,length(uc));
    uniqueIdx = cell(1,length(uc));
    uniqueNames =cell(1,length(uc));
    uniqueFolds = cell(1,length(uc));
    idx = 1;

    sets.morethanone = [];
    for i=1:length(uc)
        uniqueClasses{i} = find(idc==i);
        uniqueIdx{i} = indexes(find(idc==i));
        uniqueNames{i} = (sets.names(find(idc==i)));
        uniqueNames{i} = natsortfiles(uniqueNames{i});
        uniqueFolds{i} =  (sets.folds(find(idc==i)));
    end

    sets.morethanone = find(  cellfun(@(x) length(x) ,uniqueClasses)>=2);

    
end

