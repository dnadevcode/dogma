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
%     import OptMap.KymoAlignment.SSDAlign.ssd_algn;
%     import OptMap.KymoAlignment.NRAlign.nralign;
    tic %
    
    edgeDetectionSettings = sets.edgeDetectionSettings;

    switch sets.alignMethod
        case 2
            ssdCoef = cell(1,length(kymoStructs));
            for i=1:length(kymoStructs)
                [kymoStructs{i}.alignedKymo,ssdCoef{i}] = ssd_algn(double(kymoStructs{i}.unalignedKymo),sets);
            end   
        case 3
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
                kymoStructs{i}.alignedKymo = nralign(double(kymoStructs{i}.unalignedKymo));
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
        otherwise
    end

    % edge detection // could be skipped if only one row.
    import TempStuff.approx_main_kymo_molecule_edges;% OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
    for i=1:length(kymoStructs)
        [ kymoStructs{i}.leftEdgeIdxs,kymoStructs{i}.rightEdgeIdxs,~] = approx_main_kymo_molecule_edges(kymoStructs{i}.alignedKymo, edgeDetectionSettings);       
    end

    timePassed = toc;
    disp(strcat(['All kymos were aligned in ' num2str(timePassed) ' seconds']));

	%assignin('base','hcaSessionStruct',hcaSessionStruct)

end

