function [is_distinct,numMatchingSpecies,uniqueMatchSequences,refNums,refNumBad] = identify_discriminative(hcaSets,barcodeGenC, rezMax,thryNames)

% speciesName = {'Streptococcus pyogenes'};

% rezMax,bestBarStretch,bestLength

% unique species names
import Core.Discriminative.extract_species_name;
[uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name(thryNames);

% discriminative locations
import Core.Discriminative.disc_locations;
[refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations(rezMax, 0.05);

import Core.Discriminative.disc_true;
[is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);

import Core.Discriminative.save_disc;
if ~isempty(barcodeGenC)
    save_disc(barcodeGenC,thryNames,refNums,rezMax,bestCoefs,is_distinct, refNumBad,bestCoefsBad,hcaSets)
end

end

