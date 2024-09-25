function [uniqueSpeciesNames,idSpecies] = extract_species_name(thryNames)
    %   Args:
    %       thryNames - theory struct
    %   Returns:
    %       uc - names of unique species
    %       idc - index of where these species appear
    
    % also return unique identifiers for each species
    species = arrayfun(@(x) strsplit(thryNames{x},' '),1:length(thryNames),'un',false);
    species = cellfun(@(x) [x{2},' ',x{3}],species,'un',false);
    [uniqueSpeciesNames, ~, idSpecies] = unique( species ) ; % todo: if some species are considered the same, include that here

end

