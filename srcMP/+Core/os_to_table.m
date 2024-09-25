function [tableS] = os_to_table(synthS)
    % Converts os structure (two elements match) to a table


    % Output: tableS - posStart posStop orFactor bestBarStretch

    tableS = [-synthS.pB+1 -synthS.pB+1+synthS.lenB-1  1 1;...
        -synthS.pA+1 -synthS.pA+1+synthS.lenA-1 (synthS.or==-1)+1 synthS.bestBarStretch];
    tableS(:,1:2) = tableS(:,1:2)-tableS(1,1); % first row alw


 
end

