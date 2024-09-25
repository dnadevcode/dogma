function [ kymoStructs ] = align_kymos( sets, kymoStructs )
    % align_kymos
    % Runs alignment of kymographs. Currently two choices of
    % methods - ssdalign and nralign. Add possibility of more methods in  
    %
    %     Args:
    %         sets: settings structure
    %         unalignedKymos: unaligned kymographs
    % 
    %     Returns:
    %         alignedKymo: aligned kymographs
    %         leftEdgeIdxs: left edge indices of the molecule
    %         rightEdgeIdxs: left edge indices of the molecule
   
    disp('Starting kymo alignment...')
    
    % the two methods that could be used
    import OptMap.KymoAlignment.SSDAlign.ssd_algn;
    import OptMap.KymoAlignment.NRAlign.nralign;
    import OptMap.KymoAlignment.SPAlign.spalign;

    tic %
    
    if ~isfield(sets,'edgeDetectionSettings') % default edge detection
        sets.edgeDetectionSettings.method = 'Otsu';
    end

    edgeDetectionSettings = sets.edgeDetectionSettings;

    if ~isfield(sets,'alignMethod')
        sets.alignMethod = 1; %nralign default
    end

    switch sets.alignMethod
        case 2
            ssdCoef = cell(1,length(kymoStructs));
            for i=1:length(kymoStructs)
                [kymoStructs{i}.alignedKymo,ssdCoef{i}] = ssd_algn(double(kymoStructs{i}.unalignedKymo),sets);
            end   
                %         case 3
                %           gen_dtw_mean not inluced in hca
                %             for i=1:length(kymoStructs)
                %                 import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
                %                 [ leftEdgeIdxs,rightEdgeIdxs,~] = approx_main_kymo_molecule_edges(double(kymoStructs{i}.unalignedKymo), sets.edgeDetectionSettings);
                %      
                %                 [mean_mat,kymoStructs{i}.alignedKymo,f_ssg, mean_mm, f_mm,X] =  gen_dtw_mean(double(kymoStructs{i}.unalignedKymo),leftEdgeIdxs,rightEdgeIdxs);
                %                 kymoStructs{i}.leftEdgeIdxs = 1;
                %                 kymoStructs{i}.rightEdgeIdxs = length(kymoStructs{i}.alignedKymo);
                %             end
        case 1
            for i=1:length(kymoStructs)
                if isfield(kymoStructs{i},'unalignedBitmask')
               % Luis: Now kymostructs include featuresIdxs and
               % shiftalignedKymo
               [kymoStructs{i}.alignedKymo,~,kymoStructs{i}.shiftalignedKymo,kymoStructs{i}.alignedMask,~,~,~,kymoStructs{i}.featuresIdxs] = nralign(double(kymoStructs{i}.unalignedKymo), false,  kymoStructs{i}.unalignedBitmask(1:size(kymoStructs{i}.unalignedKymo,1),:) );                                                                         
                else
                    %Maybe delete aligned Mask from here
                    [kymoStructs{i}.alignedKymo,~,kymoStructs{i}.shiftalignedKymo,~,~,~,~,kymoStructs{i}.featuresIdxs] = nralign(double(kymoStructs{i}.unalignedKymo));
                end
  
            end  
        case 0
            edgeDetectionSettings = sets.edgeDetectionSettings;
            for i=1:length(kymoStructs)
                kymoStructs{i}.alignedKymo = double(kymoStructs{i}.unalignedKymo);
            end
            
        case 4
            % this case is based on Ostenato code - where we find the main
            % region to align
            
            % need to do a little of work here
        case 5
            % consensus based alignment
        case 3
            if ~isfield(kymoStructs{1},'unalignedBitmask')
                % based on int setting
                val = cellfun(@(x) strsplit(x.name,'int-'),kymoStructs,'un',false);
                val2 = cellfun(@(x) strsplit(x{2},'_'),val,'un',false);
                intVal = cellfun(@(x) str2num(x{1}),val2);
                [a,b] = unique(intVal);
                for ii=a
                    import OptMap.MoleculeDetection.EdgeDetection.median_filt;
                    toRun = find(intVal==ii);
                    [bitmask, posY,mat] = median_filt(cellfun(@(x) x.unalignedKymo,kymoStructs(toRun),'un',false), [5 15]);
                    for i=1:length(toRun)
              
                        kymoStructs{toRun(i)}.unalignedBitmask = bitmask{i};
                    end
                end
                
                % simplify:

            end
           for i=1:length(kymoStructs)
%                 if isfield(kymoStructs{i},'unalignedBitmask')
                    [kymoStructs{i}.alignedKymo,kymoStructs{i}.alignedMask,~,~] = ...
                        spalign(double(kymoStructs{i}.unalignedKymo),kymoStructs{i}.unalignedBitmask,sets.minOverlap,sets.maxShift,sets.skipPreAlign, sets.detPeaks);
           end  
%            figure;
%  imshowpair(imresize(kymoStructs{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  )


%                     valsAlign = sum(kymoStructs{i}.unalignedBitmask)<size(kymoStructs{i}.unalignedBitmask,1);
%                     st = find(valsAlign==0,1,'First');
%                     stop = find(valsAlign==0,1,'Last');
%                     % chop left and right
%                     kymoStructs{i}.alignedMask(:,1:st-1) = 0;
%                     kymoStructs{i}.alignedMask(:,stop+1:end) = 0;


%                     [kymoStructs{i}.alignedKymo,~,~,kymoStructs{i}.alignedMask] = nralign(double(kymoStructs{i}.unalignedKymo), [],  kymoStructs{i}.unalignedBitmask(1:size(kymoStructs{i}.unalignedKymo,1),:) );
%                 else
%                       kymoStructs{i}.unalignedBitmask = bitmask{1};
% %                     kymoStructs{i}.alignedKymo = nralign(double(kymoStructs{i}.unalignedKymo));
%                     [kymoStructs{i}.alignedKymo,kymoStructs{i}.alignedMask] = ...
%                         spalign(double(kymoStructs{i}.unalignedKymo),kymoStructs{i}.unalignedBitmask);
%                 end

        otherwise
    end

    % edge detection // could be skipped if only one row.
    import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
    for i=1:length(kymoStructs)
         if isfield(kymoStructs{i},'alignedMask')
                kymoStructs{i}.leftEdgeIdxs = arrayfun(@(frameNum) find(kymoStructs{i}.alignedMask(frameNum, :), 1, 'first'), 1:size(kymoStructs{i}.alignedMask,1));
                kymoStructs{i}.rightEdgeIdxs = arrayfun(@(frameNum) find(kymoStructs{i}.alignedMask(frameNum, :), 1, 'last'), 1:size(kymoStructs{i}.alignedMask,1));
         else
            [ kymoStructs{i}.leftEdgeIdxs,kymoStructs{i}.rightEdgeIdxs,kymoStructs{i}.alignedMask] = approx_main_kymo_molecule_edges(kymoStructs{i}.alignedKymo, edgeDetectionSettings);       
         end
     end

    timePassed = toc;
    disp(strcat(['All kymos were aligned in ' num2str(timePassed) ' seconds']));

	%assignin('base','hcaSessionStruct',hcaSessionStruct)

end

