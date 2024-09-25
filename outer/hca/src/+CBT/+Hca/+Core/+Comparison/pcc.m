function [pcc_table] = pcc(a,b)
   % pcc
    % Simplest pcc for two sequences a b of the same lengths
    %     Args:
    %         a
    %         b
    % 
    %     Returns:
    %         pcc_table
    % 
    seq1 = zscore(a,1);
    seq2 = zscore(b,1);

    pcc_table =1./length(seq1)*seq1*seq2;

end

