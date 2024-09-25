function [bar1, bar2,matchTable,lengths] = gensv(len1, lenVar, svType, kernelSigma,islinear)
    % Generates a simple linear SV, when both barcodes
    % are assumed to be linear
    
    % Args:
    %
    % len1, lenVar, svType, sets
    %
    % Returns:
    %
    % bar1, bar2,matchTable

    bar1 = normrnd(0,1,1,len1);

    switch  svType
        case 0 % random
            bar2 = normrnd(0,1,1,lenVar);
            matchTable = [];
           lengths  = [len1 lenVar];

        case 1  % insertion
             % rand position for insertion
             randPosStart = randi(len1-2)+1;
             % variation 
             variation = normrnd(0,1,1,lenVar);
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
             randPosStart = randi(len1-lenVar-2)+1;
             
             bar2 = [bar1(1:randPosStart+lenVar) bar1(randPosStart+1:randPosStart+lenVar) bar1(randPosStart+lenVar+1:len1)];
             
             matchTable = [1 randPosStart+lenVar 1 randPosStart+lenVar 1;...
                           randPosStart+1 randPosStart+lenVar randPosStart+lenVar+1 randPosStart+2*lenVar 1;...
                           randPosStart+lenVar+2 len1 randPosStart+2*lenVar+1 len1+lenVar 1];
              lengths  = [len1 len1+lenVar];

        case 4 % translocation
            % two random numbers, one bar1 and the other bar2
            randPos1 = randi(len1-lenVar-2)+1;
            randPos2 = randi(len1-lenVar-2)+1;
            
            % it could as well be that translocation is 0
            barW = [bar1(1:randPos1) bar1(randPos1+lenVar+1:len1)];
                
            bar2 = [barW(1:randPos2) bar1(randPos1+1:randPos1+lenVar) barW(randPos2+1:len1-lenVar)];
            
            if randPos1 < randPos2
                matchTable = [1 randPos1 1 randPos1 1;...
                randPos1+1 randPos1+lenVar randPos2+1  randPos2+lenVar 1;...
                randPos1+lenVar+1 randPos2+lenVar randPos1+1 randPos2 1;
                randPos2+lenVar+1 len1 randPos2+lenVar+1 len1 1];
            else
                % this case randPos1 > randPos2. But is randPos1>
                % randPos2+lenVar?
                matchTable = [1 randPos2 1 randPos2 1;...
                randPos2+1 randPos1 randPos2+lenVar+1 randPos1+lenVar 1;...
                randPos1+1 randPos1+lenVar randPos2+1 randPos2+lenVar 1;...
                randPos1+lenVar+1 len1 randPos1+lenVar+1 len1 1];        
            end
            
           lengths  = [len1 len1];

        case 5 % random
            bar2 = bar1;
            matchTable = [1 len1 1 len1 1];
           lengths  = [len1 len1];

        otherwise
            error('No available structural variation method was selected');
            
    end
    
    if islinear(1)
        bar1 =  imgaussfilt(bar1,kernelSigma);
    else
        bar1T = [bar1 bar1(1:10)];
        bar1T  =  imgaussfilt(bar1T,kernelSigma);
        bar1 = bar1T(1:length(bar1));
    end
    
    if islinear(2)
        bar2 =  imgaussfilt(bar2,kernelSigma);
    else
        bar2T = [bar2 bar2(1:10)];
        bar2T =  imgaussfilt(bar2T,kernelSigma);
        bar2 = bar2T(1:length(bar2));
    end
        
    % If we want to add noise to this barcode
%     import Rand.add_noise_to_barcode;
%     [bar2] = add_noise_to_barcode(bar2, sets);


    
end

