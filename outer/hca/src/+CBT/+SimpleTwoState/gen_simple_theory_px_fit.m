function [barLong] = gen_simple_theory_px_fit(pxStruct,gcSF,pxSize,nmpx,isC,sigma,kN,nm,cY,cN,kY,ligandLength,yoyoBindingProb,idsElt)
    % create a simple theory with Gaussian convolution at the px level
    % gen_simple_theory_bp will create at the bp level

    if nargin < 9 % 
        nm = 300;
    end

    if nargin < 10
        cY = 0.02; % concentration yoyo-1
        cN = 6; % concentration netropsin
        kY = 26; % overall kY constant
        kN = 30; % overall (max) kN constant
        ligandLength = 4;
    end
     
     
    pxPsf = nm/nmpx;
    
    % left right cut positions
%     import CBT.SimpleTwoState.px_cut_pos;
%     [pxCutLeft, pxCutRight, px] = px_cut_pos( ntSeq, gcSF, pxSize);
    pxCutLeft = pxStruct.pxCutLeft;
    pxCutRight = pxStruct.pxCutRight;
    px = pxStruct.px;

    
    k=1;
    
    % constants 
    % k=1;
    % constFun = arrayfun(@(x) netrConst*exp(-x/sigma(k)),ones(1,length(sortedSubseq)).*(4-sortv));

%     numGCs = zeros(1,px-1);
    numGCsScaled  = zeros(1,px-1);
    ll = pxCutRight(1)-pxCutLeft(1)-ligandLength+1;
%     numGCsScaled = arrayfun(@(x) sum(yoyoBindingProb(idsElt(pxCutLeft(x):pxCutRight(x)-ligandLength))/(pxCutRight(1)-pxCutLeft(1)-ligandLength+1)),1:px-1);
    % create px based map
    for i=1:px-1
        probsIndividual = yoyoBindingProb(idsElt(pxCutLeft(i):pxCutRight(i)-ligandLength));
        numGCsScaled(i) = sum(probsIndividual/ll);
%         if i==1
%            numGCs(i) =  pxCutRight(i)-pxCutLeft(i)+1 - numWsCumSum(pxCutRight(i));
%         else
%            numGCs(i) = pxCutRight(i)-pxCutLeft(i)+1 -( numWsCumSum(pxCutRight(i))-numWsCumSum(pxCutLeft(i)-1));    
%         end
%        all4mers = numWsCumSum(pxCutLeft(i)+ligandLength:pxCutRight(i))-numWsCumSum(pxCutLeft(i):pxCutRight(i)-ligandLength);
%        numATs(i,:) = arrayfun(@(x) sum(all4mers==x)/length(all4mers),ligandLength:-1:0);
    
    end

%     constFun = arrayfun(@(x) kN*exp(-x/sigma(k)),(0:ligandLength));

%     probYoyoBasedAT = cY*kY./(1+cY*kY+cN.*constFun);


%     numGCsScaled = [numATs*probYoyoBasedAT']';
% numGCsScaled = arrayfun(@(x) numATs(x,:)*netrConst*exp(-(0:4)'/sigma(k)),1:px-1);


if isC
    barLong = imgaussfilt(numGCsScaled,pxPsf,'Padding','circular');
else
    EXTRA_PX = 1000;
    barLong = imgaussfilt([zeros(1,EXTRA_PX) numGCsScaled zeros(1,EXTRA_PX)],pxPsf,'Padding','circular');
    barLong = barLong(EXTRA_PX+1:end-EXTRA_PX);
end



% theoryStr{1}.rawBarcode = barLong;
% theoryStr{1}.rawBitmask =[];
% theoryStr{1}.length = length(barLong);
% theoryStr{1}.name = 'theory';
% theoryStr{1}.isLinearTF = ~isC;
% theoryStr{1}.nmpx

end

