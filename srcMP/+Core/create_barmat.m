function [barMat,shiftVec,barM,barMzscored] = create_barmat(barS,strFac,mpInbr,LOCSnbr,barIdxst)

% here we plot a bargroup. Later this will be put into bargroup class with
% operations (such as bargroup addition/averaging)

% ix = 4
% In short: A - long barcode, i.e. all possible
% B - short barcode (with all possible re-scale factors)

barMat = cell(length(barS),1);

ix = 1;
b2 = barS(ix).rawBarcode(barS(ix).rawBitmask); % first barcode

barMat{1}{1} = 1:length(b2);
barMat{1}{2} = b2;
barMat{1}{3} = ix;

% orBars = zeros(1,length(pksUniquePos));
% reFac = zeros(1,length(pksUniquePos));

% if nargin>=12
%     f=figure;plot(zscore(b2,1),'black'); hold on
%     title('Bargroup')
% end

import Core.get_two_bars_alignment_params;

initShift = 0;

shiftVec = zeros(1,length(barS));
stopPx  = zeros(1,length(barS));
stopPx(1) = length(b2);
for ix=1:length(mpInbr)
 
    
    [~, pAnbr, pBnbr, ~, ~] = get_two_bars_alignment_params(strFac{ix},mpInbr{ix}, LOCSnbr{ix}(1), barIdxst{ix},1); 

    
    % maybe consider more reasonable names
    nextBar = barS(ix+1).rawBarcode(barS(ix+1).rawBitmask);

    
    barMat{ix+1}{1} = initShift+(-pAnbr+1+pBnbr:-pAnbr+pBnbr+length(nextBar));
    barMat{ix+1}{2} = nextBar;
    barMat{ix+1}{3} = ix;
    
    initShift = initShift -pAnbr+1+pBnbr;

    shiftVec(ix+1) = initShift;
    stopPx(ix+1)= initShift-pAnbr+pBnbr+length(nextBar);

end

shiftMin = min(shiftVec);
shifMax = max(stopPx);

barM = nan(length(barS),shifMax-shiftMin);
barMzscored = nan(length(barS),shifMax-shiftMin);

for ix=1:length(barS)
    barM(ix,barMat{ix}{1}-shiftMin+1) =  zscore(barMat{ix}{2})+5*ix;
        barMzscored(ix,barMat{ix}{1}-shiftMin+1) =  zscore(barMat{ix}{2});

end
