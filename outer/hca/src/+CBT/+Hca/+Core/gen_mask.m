function [kymoStructs,goodMols] = gen_mask(kymoStructs)


        val = cellfun(@(x) strsplit(x.name,'int-'),kymoStructs,'un',false);
%         val = cellfun(@(x) strsplit(x.name,'filter-'),kymoStructs,'un',false);

        val2 = cellfun(@(x) strsplit(x{2},'_'),val,'un',false);
        intVal = cellfun(@(x) str2num(x{1}),val2);
        [uniqueInt,b] = unique(intVal); % different intensity values

        import OptMap.MoleculeDetection.EdgeDetection.median_filt_alt;
         
% todo: extract int from info file to put to name
        means = cell(1,length(uniqueInt));
        goodMols = ones(1,length(kymoStructs));
        for ii=1:length(uniqueInt)
            toRun = find(intVal==uniqueInt(ii));
            [bitmask, posY,mat,meanD,varD,badMol] = median_filt_alt(cellfun(@(x) x.unalignedKymo,kymoStructs(toRun),'un',false), [5 15]);
            goodMols(toRun(find(badMol))) = 0;
            for i=1:length(toRun)  
                kymoStructs{toRun(i)}.unalignedBitmask = bitmask{i};
            end
        end
        
end

