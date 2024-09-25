function [f,brkPNT,leftExt,p, posBreakPointA, posBreakPointB, bar2,pcc,bar1Name] = plot_breakpoint(ix,stridx,mpI,LOCS,pksUniquePos,k,baridx,sF,A,h,B)
%   plot_breakpoint - finds and plots a breakpoint (both left and right) of one barcode vs the
%   other
%
%   Args:
%       ix,stridx,mpI,LOCS,pksUniquePos,k,baridx,sF,A,h,B
%
%   Returns:
%
    if nargin < 11 % if less arguments, we're performing the self-join. Otherwise A-B join
        B = A; % A=B
    end
    % ix = 4
    % In short: A -  time series of all possible barcodes
    % B - time series of a single barcode with all possible re-scale factors

    pos = LOCS(pksUniquePos(ix)); % position on MP.
    ts2 = k;

    import Core.get_two_bars_alignment_params;
    [bar1Name, pA, pB, rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{ts2},mpI, pos, baridx,sF);

% print out  
% [rescaleFactor, orSign]

% maybe consider more reasonable names
b2 = A(bar1Name).rawBarcode(logical(A(bar1Name).rawBitmask));
b = B(ts2).rawBarcode(logical(B(ts2).rawBitmask));

    f=figure;
tiledlayout(3,1);
nexttile
     hold on

% now we rescale b
bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);
if orSign == -1
%     bRescaled = fliplr(bRescaled);
    
%     bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);
    pAnew = pA; % remove pA and h, this should give correct position ( check for +1 )
%     pAnew = length(bRescaled)-pA-h; % remove pA and h, this should give correct position ( check for +1 )

    plot(-pAnew+1:-pAnew+length(bRescaled),zscore(fliplr(bRescaled))+5,'black')
    pBnew = pB; % remove pA and h, this should give correct position ( check for +1 )
%     pBnew = length(b2)-pB-h; % remove pA and h, this should give correct position ( check for +1 )

    plot(-pBnew+1:-pBnew+length(b2),zscore(b2)+7,'red')
%     pA
%     pB
    bar1 = fliplr(bRescaled);
else
%     f=figure;hold on
    plot(-pA+1:-pA+length(bRescaled),zscore(bRescaled)+5,'black')
    plot(-pB+1:-pB+length(b2),zscore(b2)+7,'red')
    bar1 = bRescaled;
end
    % positions of both barcodes
    pos1 = -pA+1:-pA+length(bRescaled);
    pos2 = -pB+1:-pB+length(b2);
    
    
    bar2 = b2;

    fP1 = find(pos1==1,1,'first');
    fP2 = find(pos2==1,1,'first');
    plot(zscore(bar1(fP1:fP1+h-1)),'black')
    plot(zscore(bar2(fP2:fP2+h-1)),'red')
%     
    % if -1 then we switch bar1 with bar2... for consistency
    pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';
    pcc
    pBest = 1-Zeromodel.beta_ev_cdf(pcc, 0.08*h, 1,  2*(length(bar1)-h).*(length(bar2)-h)/100, 1)

%     pvalPar2 = 2*(len1-pvalPar1).*(len2-pvalPar1)/100;

    %% Now breakpoint detection. Calculate extensions. Possibility : use VALMOD to extend?
    
    % now calculate PCC for total overlap
    
    [C,IA,IB] = intersect(pos1,pos2); % all intersecting positions between pos1 and pos2
    pcc2 = 1/length(IA) * zscore(bar1(IA),1)*zscore(bar2(IB),1)';
    pcc2 
    % PCC for full overlap. We first test if this is significant. If yes,
    % then no breakpoint detected.
    
    if pBest > 0.0001
        posBreakPointLeft = [];
        posBreakPointRight = [];
        brkPNT =[];
        leftExt =[];
        p = pBest;
        return;
    end

    
    % calculate all possible PCC extensions (so far a bit slow.).
    leftExt = [];
    curLen = [];
    curPtsA = [];
    curPtsB = [];
%     Cpts = [];
    for j=IA(1):fP1
        leftExt(j) = 1/length(j:fP1+h-1) * zscore(bar1(j:fP1+h-1),1)*zscore(bar2(IB(j)-IA(1)+1:fP2+h-1),1)';
        curLen(j) = length(j:fP1+h-1);
        curPtsA{j} = j:fP1+h-1; % points on A
        curPtsB{j} = IB(j)-IA(1)+1:fP2+h-1; %points on B
%         Cpts{j} = C(1)
    end
    % now to the right
    for j=fP1+1:IA(end)-h+1
        leftExt(j) = 1/length(fP1+1:j+h-1) * zscore(bar1(fP1+1:j+h-1),1)*zscore(bar2(fP2+1:IB(j)-IA(1)+h),1)';
        curLen(j) = length(fP1+1:j+h-1);
        curPtsA{j} =fP1+1:j+h-1;
        curPtsB{j} =fP2+1:IB(j)-IA(1)+h;
    end
%     leftExt

    % possible method 2 (no calculation needed): check the MP score and position for that location!

    % calculate MP-pval. TODO: make sure this corresponds to same for
    % synthetic data.
    pvalPar1 = 0.06*curLen; % pval distribution
    len1 = length(bar1);
    len2 = length(bar2);
    pvalPar2 = 2*(len1-curLen).*(len2-curLen)/100;
    import Zeromodel.beta_ev_pdf;
%     pval =

%     [p] = beta_ev_pdf(allPCCvals', pvalPar1, 1, pvalPar2);%y_est(k)
    p = zeros(1,length(leftExt));
    for j=1:length(leftExt)
        [p(j)] =  1-Zeromodel.beta_ev_cdf(leftExt(j), pvalPar1(j), 1, pvalPar2(j), 1);
    end

    %find first and last points where p-val is more than thresh
    pthresh = 0.01;

    brkPNT = find(p< pthresh,1,'first');
    brkPNTlast = find(p< pthresh,1,'last');

%     brkPNTlast = find(p> 0.5,1,'first');

    posBreakPointA = [find(pos2==IA(brkPNT)) find(pos2==IA(brkPNTlast))]; % IA or IB?
    posBreakPointB = [find(pos1==IB(brkPNT)) find(pos1==IB(brkPNTlast))]; % IA or IB?

%     plot(C(posBreakPointLeft),5,'blackx','MarkerSize',10)
%     plot(C(posBreakPointRight),5,'blackx','MarkerSize',10)

    posA = C(brkPNT);
    posB = C(brkPNTlast)+ h-1;%curLen(brkPNTlast);
    posMin = 3;
    posMax = 9;
    rectangle('position',[posA posMin posB-posA posMax- posMin])

    nexttile
    hold on

    plot(leftExt)
legend({'PCC score'});
nexttile
    plot(p)
legend({'p-value'});

%    plot(C(1:pos2),zscore(bar1(IA(1):IA(pos2)))-6,'black')
%     plot(C(1:pos2), zscore(bar2(IB(1):IB(pos2)))-6,'red')
    

%     % to p-val.
%     a = 0.13*curLen;
%     n = 2*(lengthsThry(i)-min_overlap);
%     pval(i) = 1-Zeromodel.beta_ev_cdf(maxcoef(i), a, 1, n, 1);
    %
    
    
%    bar2(IB)
%     figure;hold on
 
%     text(h+50,0,'Overlap')
%     text(max(C)+50,-6,'Overlap')


%     if orSign==-1
%         legend({num2str(ts2),strcat(['$$\bar ' num2str(bar1Name) '$$'])},'location','eastoutside','Interpreter','latex')
%     else 
%         legend({num2str(ts2),num2str(bar1Name)},'location','eastoutside')
%     end
%     stopPos = min(pos1(end),pos2(end));
%     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';


    saveas(f,'figs/fig3.png')


end

