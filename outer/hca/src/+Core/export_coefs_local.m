function [matFilepath] = export_coefs_local(thryNames,maxCoef,maxOr,maxPos,maxlen, bestSF, barcodeNames,matDirpath,timestamp)

% export_coefs_local
if nargin < 9
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
end

if nargin < 10
    writemode = 'overwrite';
end

T = table(thryNames');

for i=1:length(barcodeNames)
    [d,name,ext] = fileparts(barcodeNames{i});

    N = matlab.lang.makeValidName(name);

    T2 = table(maxCoef{i},maxlen{i}', maxPos{i}',bestSF{i}' ,'VariableNames',{N,strcat(['len_'  num2str(i)]),strcat(['pos_'  num2str(i)]),strcat(['stretch_'  num2str(i)])});
    T = [T T2];
end
disp('Saving ccvals table');
matFilepath = exptxt(T, matDirpath,timestamp, writemode);

    function [matFilepath,timestamp] = exptxt(T, matDirpath, timestamp, writemode)
        % Exports table as txt

        if nargin < 3
            timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
            writemode = 'overwrite';
        end
        matFilename = strcat([ 'table_' timestamp  '.txt']);
        %     matFilename2 = strcat([ 'table_' timestamp  '.dat']);


        if isequal(matDirpath, 0)
            return;
        end
        matFilepath = strcat([matDirpath, matFilename]);
        writetable(T,matFilepath,'Delimiter','\t','WriteMode',writemode)
        %     matFilepath = strcat([matDirpath, matFilename2]);

        %     writetable(T,matFilepath,'WriteRowNames',true,'WriteMode',writemode)
        fprintf('Saved table to ''%s''\n', matFilepath);


    end


end

