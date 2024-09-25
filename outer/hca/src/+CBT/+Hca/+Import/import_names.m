function [names] = import_names(files)

   try 
        fid = fopen(files); 
        names = textscan(fid,'%s','delimiter','\n'); fclose(fid);
%         for i=1:length(theories{1})
%             [sets.theoryFileFold{i}, name, ext] = fileparts(theories{1}{i});
%             sets.theoryFile{i} = strcat([name ext]);
%         end
        names = names{1};
   catch
        error('No valid theories provided, please check the provided file');
   end
   
end

