function [bar1,bar2,matchTable] = gen_synthetic_sv(len1,len2, svTypes, circ, sigma,noise)
    % function to generate synthetic structural variations
    %
    %   Args: 
    %           len1 - length of barcode without structural variation
    %           len2 - length of structural variation (is there is one
    %           svTypes - types of structural variations 
    %           circ - if the synthetic barcodes should be circular
    %           noise - what noise level should synthetic barcodes have
    %           sigma - point spread function sigma for synthetic barcodes
    
    if nargin < 6
        noise = 0; % no noise
    end
    
    if nargin < 5
        sigma = 230/159.2; % same as 159.2 nm/px experiments
    end
    
    if nargin < 4
        circ = [0 0]; % not circular
    end
    
    import CBT.Hca.Core.Comparison.pcc;

                
    numSamples = length(svTypes);
    bar1 = cell(1,numSamples);
    bar2 = cell(1,numSamples);
    matchTable = cell(1,numSamples);
%     structBar =  cell(1,numSamples);
    for i=1:numSamples
        switch svTypes(i)
            case 1
                sets.svType = 'Insertion';
            case 2
                sets.svType = 'Invertion';
            case 3
                sets.svType = 'Repeat';
            case 4
                sets.svType = 'Translocation';
            case 5
                sets.svType = 'Deletion';
            case 0
                sets.svType = 'Random';
            otherwise
                error('not a valid sv type');
        end

        % in case we wan tto compute idealized barcodes, psf is computed
        % before adding a structural variation. For noisified examples, we
        % compute this after adding structural variation
        import CBT.Hca.Import.linear_sv;
        [bar1{i}, bar2{i},matchTable{i},lengths] = linear_sv(len1, len2, svTypes(i), sigma,circ,1);

        if noise~= 0
            barRnd = normrnd(0,1,1,length(bar2{i}));        
            barRnd =  imgaussfilt(barRnd,sigma);
            tf = @(x) pcc( bar2{i},(1-x)* bar2{i}+barRnd*x)-(1-noise);
            x0 = 0.5;  x = fzero(tf,x0);
            bar2{i} = (1-x)* bar2{i}+barRnd*x;
        end
%         
        % circularly shift
        if circ(1)
            % shift circularly
            % now shift one of these. This should be tested via tests before
            % generating final figures for the paper
            shift = randi(length(bar1{i}));
            bar1{i}  = circshift(bar1{i}, [0,-shift]);
           for k=1:size(matchTable{i},1)
                matchTable{i}(k,1:2) = mod( matchTable{i}(k,1:2)-1-shift,length(bar1{i}))+1; 
           end
           shift2 = randi(length(bar2{i}));

            bar2{i}  = circshift(bar2{i}, [0,-shift2]);
           for k=1:size(matchTable{i},1)
                matchTable{i}(k,3:4) = mod( matchTable{i}(k,3:4)-1-shift2,length(bar2{i}))+1; 
           end
        end
        
        %% Could also add some stretching to barcode
        
%                 % save DATA
%         rS = strcat([folder num2str(i) '_data_seq.txt']);
%         fid = fopen(rS,'w');
%         fprintf(fid, '%5.4f ', structBar{i}.bar2  );
%         fprintf(fid, '\n');
%         fclose(fid);
%         structBar{i}.name1 = rS;
% 
%        % save QUERY (this includes the variation)
%         rS = strcat([folder num2str(i) '_query_seq.txt']);
%         fid = fopen(rS,'w');
%         fprintf(fid, '%5.4f ',structBar{i}.bar1  );
%         fprintf(fid, '\n');
%         fclose(fid);
%         structBar{i}.name2 = rS;

    end


end

