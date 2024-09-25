function [is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies)

    %   Args:
    %
    %   Returns:
    %       is_distinct - whether match is distinct
    %       numMatchingSpecies - number of matching sequences
    %       uniqueMatchSequences - unique matchin sequences


    %Returns:   
    %   truePositives - true positives based on allspecies
    %   discSpecies - how many strains is it discriminative to
    %   discAll - binary vector of all strains its discriminative to
    %   allNums - all discriminative locations
    %   refNums - refnums of discriminative locations
    %   signMatch - numver of unique discriminative locations

% idxCor = zeros(1,length(speciesName)); % can be multiple species/ strains we are interested in
% for i = 1:length(speciesName)
%     idxCor(i) = find(cellfun(@(x) isequal(speciesName{i},x),uniqueSpeciesNames));
% end

numMatchingSpecies = zeros(1,length(refNums));
uniqueMatchSequences = cell(1,length(refNums));
is_distinct = zeros(1,length(refNums));
for i = 1:length(refNums)
    idxSpecies = idSpecies(refNums{i});
    uniqueMatchSequences{i} = unique(idxSpecies);

    numMatchingSpecies(i) = length(uniqueMatchSequences{i});
    if numMatchingSpecies(i)==1
        is_distinct(i) = 1;
    end

end


% allSpecies = find(speciesLevel);
% allSpecUnique = unique(idc(find(speciesLevel)));
% try
% idc(idc==allSpecUnique(2)) = allSpecUnique(1);
% catch
% end
% discAll = cellfun(@(x) ismember(x,allSpecies),refNums,'UniformOutput',false);
% 
% discSpecies = cellfun(@(x) sum(ismember(x,allSpecies)==0),refNums,'UniformOutput',true);
% 
% truePositives = sum(discSpecies==0); % also count false positives?;
% 
% discAny = cellfun(@(x) unique(idc(x)),refNums,'UniformOutput',false);
% 
% positives = cellfun(@(x) length(x),discAny);
% fp = sum(positives==1)-truePositives;

end

