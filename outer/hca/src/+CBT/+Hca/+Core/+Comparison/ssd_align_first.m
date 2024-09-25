function [ alignedKym,unAlignedKymoMoleculeMask,alignedKymoMoleculeMask,backgroundKym, ssdCoef] = ssd_align_first(kymoToAlign,barcodeConsensusSettings,edgeDetectionSettings,sets)

    % input
    % kymoToAlign,barcodeConsensusSettings,edgeDetectionSettings,sets
    
    % output 
    % alignedKym,unAlignedKymoMoleculeMask,alignedKymoMoleculeMask,backgroundKym, ssdCoef
    
    if nargin < 4
        stretchPar = 1;
%         alowedShift = 10*stretchPar;
%         pxMax = 500*stretchPar;
%         shiftInd = 40*stretchPar;
    else
        stretchPar = sets.stretchPar;
%         alowedShift = sets.alowedShift;
%         pxMax = sets.pxMax;
%         shiftInd = sets.shiftInd; 
    end

    ssdCoef.left= [];
    ssdCoef.cor = [];
    ssdCoef.tot = [];
    ssdCoef.shift = [];
    import CBT.Hca.Core.Comparison.SSD_fft_all;
    
    filterSize = barcodeConsensusSettings.psfSigmaWidth_nm/barcodeConsensusSettings.prestretchPixelWidth_nm;
   
    % generate first barcode        
   % kym1 = imgaussfilt(kymoToAlign(1,:), filterSize);
    
    edgeDetectionSettings.otsuApproxSettings.globalThreshTF = false;
    % instead, approx left and right edges for each row separately
    import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
   % [leftEdgeIdxs, rightEdgeIdxs, alignedKymoMoleculeMask] = approx_main_kymo_molecule_edges(kymoToAlign, edgeDetectionSettings);
   
   % Do edge detection row by row!
   
%    
   leftEdgeIdxs = zeros(1,size(kymoToAlign,1));
   rightEdgeIdxs = zeros(1,size(kymoToAlign,1));
  % unAlignedKymoMoleculeMask = nan(size(kymoToAlign));
  % kym = cell(1,size(kymoToAlign,1));
%    for i=1:size(kymoToAlign,1)
%        kym{i} = imgaussfilt(kymoToAlign(i,:),filterSize);
%       [leftEdgeIdxs(i), rightEdgeIdxs(i), unAlignedKymoMoleculeMask(i,:)] = approx_main_kymo_molecule_edges(kym{i}, edgeDetectionSettings);
%    end
%    
% %    
   unAlignedKymoMoleculeMask = zeros(size(kymoToAlign));
   kym = cell(1,size(kymoToAlign,1));
   for i=1:size(kymoToAlign,1)
       
        kym{i} = imgaussfilt(kymoToAlign(i,:),stretchPar*filterSize);
        [idx1,~] = kmeans(kym{i}',2);
        [~,~,idx1] = unique(idx1,'stable');
        idx1 = idx1-1;
        
        %
        leftEdgeIdxs(i) = find(idx1,1,'first');
        rightEdgeIdxs(i) = find(idx1,1,'last');
        %
        unAlignedKymoMoleculeMask(i,leftEdgeIdxs(i):rightEdgeIdxs(i))=ones(1,length(rightEdgeIdxs(i))-length(leftEdgeIdxs(i))+1);
                

        %unAlignedKymoMoleculeMask(i,:)=idx1;
       %rez = kmeans(kym{i}',2);

     % [leftEdgeIdxs2(i), rightEdgeIdxs2(i), alignedKymoMoleculeMask2(i,:)] = approx_main_kymo_molecule_edges(kym{i}, edgeDetectionSettings);
   end
% %    
    alignedKymoMoleculeMask = unAlignedKymoMoleculeMask;

%    	[~,leftEdgeIdxs]=max(alignedKymoMoleculeMask,[],2);
%     [~,rightEdgeIdxs]=min(alignedKymoMoleculeMask(:,max(leftEdgeIdxs):end),[],2);
% 	rightEdgeIdxs = rightEdgeIdxs+max(leftEdgeIdxs)-2;

%    for i=1:20
%     figure,plot(kym{i})
%     hold on
%     plot(leftEdgeIdxs(i):rightEdgeIdxs(i),kym{i}(leftEdgeIdxs(i):rightEdgeIdxs(i)))
%    end

%     kym = nan(size(kymoToAlign,1),stretchPar*pxMax);
%     kymoToAlign
   % figure,plot
   
%     backgroundKym = nan(size(kymoToAlign));
%     backgroundKym(~alignedKymoMoleculeMask) = kymoToAlign(~alignedKymoMoleculeMask);
%     
%    %[leftEdgeIdxs, rightEdgeIdxs, alignedKymoMoleculeMask] = approx_main_kymo_molecule_edges(kym1, edgeDetectionSettings);
   

    
    kym1mol = kym{1}(leftEdgeIdxs(1):rightEdgeIdxs(1)); 
    edgePixels = round(barcodeConsensusSettings.prestretchUntrustedEdgeLenUnrounded_pixels);  
    kym1bit = ones(1,length(kym1mol));
    kym1bit(1:edgePixels) = zeros(1,length(edgePixels));
    kym1bit(end-edgePixels+1:end) = zeros(1,length(edgePixels));
    
%     figure,plot(kym1)
%     hold on
%     plot(leftEdgeIdxs:rightEdgeIdxs,kym1(leftEdgeIdxs:rightEdgeIdxs))
%     
%     kym1int = interp1(kym1mol, linspace(1,length(kym1mol),stretchPar*length(kym1mol)));
%     kym1bitInt = kym1bit(round(linspace(1,length(kym1bit),length(kym1bit)*stretchPar)));

    %kym1bitInt = interp1(kym1bit,linspace(1,length(kym1),stretchPar*length(kym1)));
    

    alignedKym = nan(size(kymoToAlign));
    %firstRow = kymoToAlign(1,leftEdgeIdxs(1):rightEdgeIdxs(1));
    alignedKym(1,:) = kymoToAlign(1,:);
    
    %totalShift = 0;
            tic
    %for i=2:2
    for i=2:size(kymoToAlign,1)
         
         kym2mol = kym{i}(leftEdgeIdxs(i):rightEdgeIdxs(i)); 
         kym2bit = ones(1,length(kym2mol));
         kym2bit(1:edgePixels) = zeros(1,length(edgePixels));
         kym2bit(end-edgePixels+1:end) = zeros(1,length(edgePixels));
        % kym2int = interp1(kym2mol, linspace(1,length(kym2mol),stretchPar*length(kym2mol)));
         %kym2bitInt = kym2bit(round(linspace(1,length(kym2bit),length(kym2bit)*stretchPar)));

        [ssdV,ssdB, indices] = CBT.Hca.Core.Comparison.SSD_fft_all(kym1mol/mean(kym1mol),kym2mol/mean(kym2mol), kym1bit, kym2bit, 2*edgePixels);

         %[ssdV,ssdB, indices] = CBT.Hca.Core.Comparison.SSD_fft_all(kym1int,kym2int, kym1bitInt, kym2bitInt, alowedShift);
         
         ssdCoef.left = [ ssdCoef.left; ssdV];
         [a,b] = min(ssdV);
         relShift = indices(b);

         ssdCoef.cor = [ ssdCoef.cor b];
         %ssdCoef.right{i} = ssdB;
         
      %   nextRow =[nan(1,shiftInd) kymoToAlign(i,leftEdgeIdxs(i):rightEdgeIdxs(i)) nan(1,shiftInd)];
         ssdCoef.shift  =[ssdCoef.shift relShift+leftEdgeIdxs(1)-leftEdgeIdxs(i)];
         alignedKym(i,:) = circshift(kymoToAlign(i,:),[0,relShift+leftEdgeIdxs(1)-leftEdgeIdxs(i)]);
         alignedKymoMoleculeMask(i,:) = circshift(alignedKymoMoleculeMask(i,:),[0,relShift+leftEdgeIdxs(1)-leftEdgeIdxs(i)]);
       %  alignedKym(i,1:length(nextRow)*stretchPar) = circshift( interp1(nextRow, linspace(1,length(nextRow),stretchPar*length(nextRow))),[0,totalShift]);
%          kym1int = kym2int;
%          kym1bitInt = kym2bitInt;

    end
    backgroundKym = nan(size(alignedKym));
    backgroundKym(~alignedKymoMoleculeMask) = alignedKym(~alignedKymoMoleculeMask);
    
    %alignedKymBitMask = ~isnan(alignedKym);
   % figure,imshow(alignedKym,[])

end

