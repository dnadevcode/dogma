function [mp1,mpI1,indexesT,indexesT2,telapsed] = compare_rescaled_theory(bar, theoryStruct,sF, wmin,numWorkers)
% compare_rescaled_theory Compare barcodes in barto re-scaled theory

if nargin < 5
    numWorkers = 30;
end

tmpNameB = [tempname(pwd),'.txt'];
tmpNameA = [tempname(pwd),'.txt'];

tmpNameMP = [tempname(pwd),'.txt'];
tmpNameMPI = [tempname(pwd),'.txt'];

theoryStruct.rawBitmask = zeros(1,length(theoryStruct.rawBarcode));
import Core.rescale_barcode_data;
[tsRescaled] = rescale_barcode_data({theoryStruct},sF);


tempCell = cell(1, 4*numel(tsRescaled{1}.rescaled)); % 1 is forward, 3 is reverse, 2 and 4 is nan
tempCell(1:4:end) = cellfun(@(x) x.rawBarcode,tsRescaled{1}.rescaled,'un',false);
tempCell(3:4:end) = cellfun(@(x) fliplr(x),tempCell(1:4:end) ,'un',false);
tempCell(2:2:end) = {NaN};
% Concatenate the cell array into a single vector
vecConcat = cat(2, tempCell{:});

% convert vec to integers for simpler calculation
newvecB = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);

writematrix(newvecB,tmpNameB,'Delimiter',' ');

%%
% blend forward and reverse
tempCell2 = cell(1, 2*numel(bar)); % 1 is forward, 3 is reverse, 2 and 4 is nan
tempCell2(1:2:end) =  cellfun(@(x) x.rawBarcode,bar,'un',false);
tempCell2(2:2:end) = {NaN};
% Concatenate the cell array into a single vector
vecConcat2 = cat(2, tempCell2{:});

% convert vec to integers for simpler calculation
newvecA2 = round((vecConcat2-min(vecConcat2))/(max(vecConcat2)-min(vecConcat2))*256);

writematrix(newvecA2,tmpNameA,'Delimiter',' ');

%     strjoin([theoryStruct(1:2).rawBarcode])
com = strcat(['SCAMP --window=' num2str(wmin) ' --input_a_file_name='...
    tmpNameA ' --input_b_file_name=' ...
    tmpNameB ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
    ' --output_a_file_name=' tmpNameMP ...
    ' --output_a_index_file_name=' tmpNameMPI]);

% tic
t0 = tic;

[a, val ] = system(com);
% toc
telapsed = toc(t0);



fid = fopen(tmpNameMP);
raw2 = textscan(fid, '%s ');
fclose(fid);
nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
mp1 = nan(length(nonanValues),1);
mp1(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');

if nargout >=2
    mpI1 = importdata(tmpNameMPI);
end

% if we want to delete temporary files
delete(tmpNameB);
delete(tmpNameA);
delete(tmpNameMP);
delete(tmpNameMPI);


% newvecB = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);

allLens = cellfun(@(x) length(x.rawBarcode),bar);
% nanIndices = find(isnan(newvecA2));

% first barcode positions
indexesT = cell(1,numel(allLens));
% Split the vector based on NaN delimiters
for i = 1:numel(allLens) % every two, since doesn't matter if we look at forward or reverse
    if i == 1 
        if i==numel(allLens)
            indexesT{i} = [1 allLens(1)-wmin+1];
        else
        indexesT{i} = [1 allLens(1)];
        end
    else
        if i==numel(allLens)
            indexesT{i}  =  [indexesT{i-1}(2)+2 indexesT{i-1}(2)+2+allLens(i)-1-wmin];
        else
            indexesT{i}  =  [indexesT{i-1}(2)+2 indexesT{i-1}(2)+1+allLens(i)];
        end
    end
end

% second barcode positions
indexesT2 = cell(1,numel(tsRescaled{1}.rescaled));

% % Split the vector based on NaN delimiters
% for i = 1:numel(allLens) % every two, since doesn't matter if we look at forward or reverse
%     if i == 1
%         indexesT{i} = [1 allLens(1)];
%     else
%         if i==numel(allLens)
%             indexesT{i}  =  [indexesT{i-1}(2)+2 indexesT{i-1}(2)+2+allLens(i)-1-wmin];
%         else
%             indexesT{i}  =  [indexesT{i-1}(2)+2 indexesT{i-1}(2)+1+allLens(i)];
%         end
%     end
% end



%
% %% For rest of the theories
% restOfTheories = find(idSpecies~=idSeq);
% thryToCalc = theoryStruct(idSpecies~=idSeq);
%
% % thryToCalc = theoryStruct(1:100);
% tempCell = cell(1, 4*numel(thryToCalc)); % 1 is forward, 3 is reverse, 2 and 4 is nan
% tempCell(1:4:end) = {thryToCalc(:).rawBarcode};
% tempCell(3:4:end) = cellfun(@(x) fliplr(x),{thryToCalc(:).rawBarcode},'un',false);
% tempCell(2:2:end) = {NaN};
% % Concatenate the cell array into a single vector
% vecConcat = cat(2, tempCell{:});
%
% % convert vec to integers for simpler calculation
% newvecB = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);
%
% nanIndices = find(isnan(newvecB));
%
% % Preallocate memory
% % indexesT = cell(1,numel(nanIndices)/2);
% % % Split the vector based on NaN delimiters
% % for i = 1:numel(nanIndices) % every two, since doesn't matter if we look at forward or reverse
% %     if i == 1
% %         indexesT{i} = [1 nanIndices(i)-1];
% %     else
% %         indexesT{i}  = [nanIndices(i-1)+1 nanIndices(i)-1];
% %     end
% % end
%
% %     barIndex = zeros(1,length(newvec));
%     indexesT = cell(1,numel(nanIndices)/2);
%     % Split the vector based on NaN delimiters
%     for i = 1:2:numel(nanIndices) % every two, since doesn't matter if we look at forward or reverse
%         if i == 1
% %             barIndex(1:nanIndices(i+1)-1) = i*ones(1, nanIndices(i+1)-1);
%             indexesT{i} = [1 nanIndices(i+1)-1];
%         else
% %             barIndex(nanIndices(i-1)+1:nanIndices(i+1)-1) = ceil(i/2)*ones(1, nanIndices(i+1)-1-(nanIndices(i-1)+1)+1);
%             indexesT{ceil(i/2)}  = [nanIndices(i-1)+1 nanIndices(i+1)-1];
%         end
%     end
%
%
%
% writematrix(newvecB','barB.txt','Delimiter',' ');
%
% % numWorkers = 4;
%
% com= strcat(['SCAMP --window=' num2str(wmin) ' --input_a_file_name='...
%     fullfile(pwd,'barA.txt') ' --input_b_file_name=' ...
%     fullfile(pwd,'barB.txt') ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
%     ' --output_a_file_name=' 'bar_mp' ...
%     ' --output_a_index_file_name=' 'bar_index']);
%
% tic
% [a,val ] = system(com);
% toc
end

