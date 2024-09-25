function [nmbp, nmpx] = extract_extension_params(dirName)
% extract_extension_params
% Every folder name has nmpx and nmPerbp parameters in its name, we extract
% them using this function 
% smart extract nmPerbp from folder name

spltName = strsplit(dirName ,'_');
spltName2 = strsplit(spltName{end},'nm');
nmbp = str2double(spltName2{1}); % extract nm per bp


nmpxOptions = [110 130 159 254];
% extract nm/px, sometimes this is saved in bpperpx
spltName3 = strsplit(spltName{end-1},'nm');
if numel(spltName3)==1
    spltName3 = strsplit(spltName{end-1},'bp');
    [minDistance, indexOfMin] = min(abs(nmpxOptions-nmbp*str2double(spltName3{1})));
    nmpx = nmpxOptions(indexOfMin);
else
    nmpx = str2double(spltName3{1});

end



end

