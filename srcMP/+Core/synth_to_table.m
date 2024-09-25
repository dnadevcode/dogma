function [tableS] = synth_to_table(synthS)
    % Converts synth structure to a table


    % Output: tableS - posStart posStop orFactor bestBarStretch

    tableS = zeros(length(synthS),4);

    % Initially, bestBarStretch in synthS table is w.r.t. the theory
    % barcode

    for i=1:length(synthS)
        tableS(i,:) = [synthS{i}.pos(1) synthS{i}.pos(1)+synthS{i}.lengthMatch-1 synthS{i}.or(1) synthS{i}.bestBarStretch];
%         sF = bbS;
% 
%         newLens = tableS(i,2)-tableS(i,1)+1;    
%     
%         tableS(i,1) = round(tableS(i,1)/sF);
%         tableS(i,4) = tableS(i,4)*sF;
%         tableS(i,2) = round(newLens*sF+tableS(i,1)-1);
    end

    % 

end

