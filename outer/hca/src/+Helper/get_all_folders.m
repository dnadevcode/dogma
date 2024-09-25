function [barN, twoList,foldname] = get_all_folders(dirName)


if nargin < 1
    dirName = uigetdir();
end
% '/proj/snic2022-5-384/users/x_albdv/data/discriminative/pyo/';

dr = dir(dirName);
folderNames = {dr([dr.isdir]).name};
dr = dr(~ismember(folderNames ,{'.','..'}));

numDirs = length(dr);
subdirs = zeros(1,numDirs);
barN = cell(1,numDirs);
for i=1:numDirs
    dr2 = dir(fullfile(dr(i).folder,dr(i).name));
    dr2 = dr2([dr2.isdir]);
    folderNames = {dr2([dr2.isdir]).name};
    dr2 = dr2(~ismember(folderNames ,{'.','..'}));
    subdirs(i) = numel(dr2);
    for j=1:numel(dr2)
        barN{i}(j) = numel(dir(fullfile(dr2(j).folder,dr2(j).name,'kymos','*.tif')));
        foldname{i}{j} = dr2(j).name;
    end
end


twoList = zeros(sum(subdirs),2);
id=1;
for i=1:length(subdirs)
    for j=1:subdirs(i)
        twoList(id,:) = [i j];
        id = id+1;
    end
end


end

