function [bitmaskIndices , cutOutStartPos] = generate_random_cutouts(cutoutSize , ....
           noOfCutouts , bitmaskCellArray  )
%
% Provides random cut-outs of specific size 
% from a list of input bitmasks
% 
% Input:
% cutoutSize = size of cut-out region (pixels)
% noOfCutouts = number of cut-outs
% bitmaskCellArray = Input set of bitmasks (cell array)
% 
%
% Output:
% barcodeIndices = random barcode indices
% cutOutStartPos = random start postions 
%
% Dependencies: None.
%

% Figure out lengths, effective lengths and number of possible
% cut-outs for each barcode
N = length(bitmaskCellArray);  % The number of bitmasks we have
startPosBitmask=zeros(1,N);
endPosBitmask=zeros(1,N);
lengths=zeros(1,N);
effectiveLengths=zeros(1,N);
noOfPossibleCutouts=zeros(1,N); % for a given experimental barcode, 
                                % counts how many possible cut-outs there are
for bitmaskIdx=1:N 
    %
    lengths(bitmaskIdx)=length(bitmaskCellArray{bitmaskIdx});
    %
    bitmaskIdxTemp=find(bitmaskCellArray{bitmaskIdx}==1);
    startPosBitmask(bitmaskIdx)=min(bitmaskIdxTemp);
    endPosBitmask(bitmaskIdx)=max(bitmaskIdxTemp);
    effectiveLengths(bitmaskIdx)=endPosBitmask(bitmaskIdx)-startPosBitmask(bitmaskIdx) + 1;
    %
    noOfPossibleCutouts(bitmaskIdx)=effectiveLengths(bitmaskIdx) - cutoutSize + 1;
    %
end
if max(noOfPossibleCutouts) < 1
   disp('Warning: the cut-out size is larger than the effective lengths of all barcodes')
end



% Select randomly barcodes with probability proportional 
% to the possible number of cut-outs from that barcode. 
% Provide probability = 0 to a barcode if the possible 
% number of cut-outs is smaller than 1.
weightVec=zeros(1,N);
for bitmaskIdx=1:N
    if noOfPossibleCutouts(bitmaskIdx) < 1
        weightVec(bitmaskIdx)=0;
    else
        weightVec(bitmaskIdx) = noOfPossibleCutouts(bitmaskIdx);
        %weightVec(barcodeIdx)=effectiveLength(barcodeIdx);
        %weightVec(barcodeIdx)=lengths(barcodeIdx);
    end
end
% Turn weights into probabilities and generate 
% random numbers i, with probabilities p_i (with replacement)
probVec=weightVec/sum(weightVec);
iVec=1:N;   
bitmaskIdxSelected = randsample(iVec,noOfCutouts,true,probVec);


% Find a random region of size cutOutSize in the selected barcodes
for counter=1:noOfCutouts
    
    bitmaskIdx=bitmaskIdxSelected(counter);
    
    firstAllowedStartPos=startPosBitmask(bitmaskIdx);
    lastAllowedStartPos=endPosBitmask(bitmaskIdx)-cutoutSize+1;
    randStartPos(counter) = randi([firstAllowedStartPos lastAllowedStartPos],1);
    
end

% Return output
bitmaskIndices = bitmaskIdxSelected;
cutOutStartPos = randStartPos;

end

