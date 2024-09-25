function [mpIpos,mpI] = convert_subfragment_pos_to_query(mpI,q,r,mpOr)
    % converts positions outputed by MP - which is best position of
    % subfragment, to relative orientations of query vs data (which are
    % more useful for plotting and overall output, unless we're discarding 
    % the rest of the barcode)
    
    
    %   Args:
    %       mpI - mpI positions
    %       q - length of query
    %       r - length of alignment
    %       mpOr - orientation of query
    
    %   Returns:
    %       mpIpos - overall positions of query vs data
  
    mpIpos = zeros(length(mpI),1);
    for i=1:length(mpI)
        if mpOr(i) == 0
            mpIpos(i) = mpI(i)-i+1;
        else
            mpIpos(i) = mpI(i)-(q-r+1-i);
            mpI(i) = mpI(i)-(q-r); % add to correct idxpos
        end
    end
    
%     mpIpos = zeros(length(mpI),1);
%     for i=1:length(mpI)
%         if mpOr(i) == 0
%             mpIpos(i) = mpI(i)-i+1;
%         else
%             mpIpos(i) = mpI(i)-i+(q-r-(i-1));
%         end
%     end
%     
    
    

end

