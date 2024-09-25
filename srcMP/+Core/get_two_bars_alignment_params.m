function [bar1, pA, pB, rescaleFactor, orSign] = get_two_bars_alignment_params(stridx,mpI, pos, baridx,sF)
    % function to get two barcode orientation parameters
    
    %       Args:
    %    stridx,mpI, pos, baridx,sF
    
    %       Returns:
    %    bar1, pA, pB, rescaleFactor, orSign
    
    
    posA = mpI(pos)+1;

    % rescale factor on time series B
    rescaleFactor = stridx(posA); 
    % orientation
    orSign = sign(rescaleFactor);
        
    % barcode on time series A
    bar1 = baridx(pos);
    
    % first position of barcode
    v1 = find(baridx == bar1,1,'first');
    % first position of rescale factor
    v2 = find(stridx==rescaleFactor,1,'first');

    % position of overlap on A
    pA = posA-v2;
    % position of overlap on B
    pB = pos-v1;
    
    % final scale factor
    rescaleFactor = sF(abs(rescaleFactor));
    


end

% %% speed test
% tic
%     v1 = find(baridx == bar1,1,'first');
%     v2 = find(stridx==rescaleFactor,1,'first');
% toc
% 
% [a,b] = unique(stridx);
% [c,d] = unique(baridx);
% 
% for i=1:100
% tic
% v1=d(bar1);
% v2=b(rescaleFactor-a(1));
% t(i)= toc
% end

