function [bar1, bar2,matchTable,lengths] = linear_sv(len1, lenVar, svType, sigma,circ,prePsf)
    % Generates a simple linear SV, when both barcodes
    % are assumed to be linear
    
    % Args:
    %
    % len1, lenVar, svType, sets
    %
    % Returns:
    %
    % bar1, bar2,matchTable
    
    if nargin < 6
        prePsf = 0;
    end

    bar1 = normrnd(0,1,1,len1);
    if prePsf
        bar1 =  imgaussfilt(bar1,sigma);
    end

    switch  svType
        case 0 % random
            bar2 = normrnd(0,1,1,lenVar);
%             if prePsf % if circ might want to circularly convolvewith
%   %          gaussian
            bar2 =  imgaussfilt(bar2,sigma);
%             end
            matchTable = [];
            lengths  = [len1 lenVar];
            return;

        case 1  % insertion
             % rand position for insertion
             randPosStart = lenVar+randi(len1-2-2*lenVar)+1;
             % variation 
             variation = normrnd(0,1,1,lenVar);
             if prePsf
                 variation = imgaussfilt(variation,sigma);
             end
             % bar2
             bar2 = [bar1(1:randPosStart) variation bar1(randPosStart+1:end)];
             
             % matchtable - 2 parts
             matchTable = [1 randPosStart 1 randPosStart 1;...
                           randPosStart+1 len1 randPosStart+lenVar+1 len1+lenVar 1];
           lengths  = [len1 len1+lenVar];
                       
        case 2 % invertion
             randPosStart = randi(len1-lenVar-2)+1;
             
             bar2 = [bar1(1:randPosStart) bar1(randPosStart+lenVar:-1:randPosStart+1) bar1(randPosStart+lenVar+1:len1)];
             
             % matchtable - 3 parts
             matchTable = [1 randPosStart 1 randPosStart 1;...
                           randPosStart+1 randPosStart+lenVar randPosStart+lenVar randPosStart+1 2;...
                           randPosStart+lenVar+1 len1 randPosStart+lenVar+1 len1 1];
           lengths  = [len1 len1];

        case 3 % repetition
             randPosStart = randi(len1-lenVar-2)+1; % position of random 
             
             % bar2 of three parts: 1:r
             bar2 = [bar1(1:randPosStart+lenVar-1) bar1(randPosStart:randPosStart+lenVar-1) bar1(randPosStart+lenVar:len1)];
             
             matchTable = [1 randPosStart+lenVar-1 1 randPosStart+lenVar-1 1;...
                           randPosStart randPosStart+lenVar-1 randPosStart+lenVar randPosStart+2*lenVar-1 1;...
                           randPosStart+lenVar len1 randPosStart+2*lenVar len1+lenVar 1];
              lengths  = [len1 len1+lenVar];

        case 4 % translocation
            % two random numbers, one bar1 and the other bar2
            % How far away we move should not be smaller than lenVar
            % itself!
            % take rand1 < rand2 
%             st=[];
%             st2 = [];
%             for i=1:10000
% %                 randPos1 = lenVar+randi(len1-3*lenVar+1);
%              randPos1 = lenVar+randi(len1-3*lenVar);
%             randPos2 = randPos1+lenVar+randi(len1-randPos1-2*lenVar+1);
%             
% %                 st= [st (len1-randPos1+1)];
%                   st= [st (len1-randPos2+1)];
%                 st2= [st2 (randPos2-randPos1+1)];
% 
%             end
%             min(st)
%             max(st2)

            % We don't allow the SV be in the first lenVar pixels, and we
            % also want the smaller part not to be smaller than lenVar
            % rand1 has to leave 2*lenVar at the end
            randPos1 = lenVar+randi(len1-4*lenVar);
            randPos2 = randPos1+lenVar+randi(len1-randPos1-3*lenVar+1);
          
%             randPos1+lenVar+randi(len1-randPos1-lenVar-2)-1;
%             st= [st (randPos2-randPos1)];
%             end
            % it could as well be that translocation is 0
            barW = [bar1(1:randPos1) bar1(randPos1+lenVar+1:len1)];
                
            bar2 = [barW(1:randPos2) bar1(randPos1+1:randPos1+lenVar) barW(randPos2+1:len1-lenVar)];
            
%             if randPos1 < randPos2
            matchTable = [1 randPos1 1 randPos1 1;...
            randPos1+1 randPos1+lenVar randPos2+1  randPos2+lenVar 1;...
            randPos1+lenVar+1 randPos2+lenVar randPos1+1 randPos2 1;
            randPos2+lenVar+1 len1 randPos2+lenVar+1 len1 1];
%             else
%                 % this case randPos1 > randPos2. But is randPos1>
%                 % randPos2+lenVar?
%                 matchTable = [1 randPos2 1 randPos2 1;...
%                 randPos2+1 randPos1 randPos2+lenVar+1 randPos1+lenVar 1;...
%                 randPos1+1 randPos1+lenVar randPos2+1 randPos2+lenVar 1;...
%                 randPos1+lenVar+1 len1 randPos1+lenVar+1 len1 1];        
%             end
            
           lengths  = [len1 len1];

        case 5 % random
            bar2 = bar1;
            matchTable = [1 len1 1 len1 1];
            lengths  = [len1 len1];

        otherwise
            error('No available structural variation method was selected');
            
    end
    
    % if need circular, merge first and last row
    if circ(1)
        % merge last and first row. 
        
        matchTable(1,:) = [matchTable(end,1) matchTable(1,2) matchTable(end,3) matchTable(1,4) matchTable(1,5)];
        matchTable(end,:) = [];
        % for brevity (to make it more difficult) also circularly shift both barcodes
    end
    
    
    if ~prePsf
        bar1 =  imgaussfilt(bar1,sigma);
        bar2 =  imgaussfilt(bar2,sigma);
    end
    
    % If we want to add noise to this barcode
%     import Rand.add_noise_to_barcode;
%     [bar2] = add_noise_to_barcode(bar2, sets);


    
end

