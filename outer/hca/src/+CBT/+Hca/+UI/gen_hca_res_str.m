function [ resultStrucC ] = gen_hca_res_str( len1,xcorr,StretchedBar,markers)

  
        resultStrucC.maxcoef = [];
        resultStrucC.pos =[];
        resultStrucC.or = [];

        for i=len1+1
            [f,s] =max(xcorr{i});
            [mV,id]=max(f);

             [ b, ix ] = sort( f(:), 'descend' );
             indx = b(1:3)' ;
             resultStrucC.maxcoef = [resultStrucC.maxcoef;indx];
           %  resultStruc.pos = [resultStruc.pos; ix(1:3)'];
             resultStrucC.or = [resultStrucC.or; s(ix(1:3)')];

             if s(ix(1:3)') == 1
                resultStrucC.pos = [resultStrucC.pos; ix(1:3)'];
             else
                resultStrucC.pos = [resultStrucC.pos; ix(1:3)'-length(StretchedBar{i})];
             end

        end
            p3 = plot(resultStrucC.maxcoef(1,:),len1+1,'ob');
            p3(1).Marker = markers(1);
            p3(2).Marker = markers(2);
            p3(3).Marker = markers(3);
            
            p4 = plot(resultStrucC.maxcoef(2,:),len1+1,'or');
            p4(1).Marker = markers(4);
            p4(2).Marker = markers(5);
            p4(3).Marker = markers(6);


end

