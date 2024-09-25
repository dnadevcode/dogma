function [cache] = display_results_ui(ts, cache)
     if nargin < 2
        cache = containers.Map();
     end
%     
% Display_results_UI
    hcaSessionStruct = cache('hcaSessionStruct');
    
    len1=length(hcaSessionStruct.barcodeGen);

    markers = ['o';'s';'x';'+';'d';'v'];
    maxCorCoefs = [];
    resultStruc.maxcoef = [];
    resultStruc.pos =[];
    resultStruc.or = [];
    
    resultStruc2.maxcoef = [];
    resultStruc2.pos =[];
    resultStruc2.or = [];
    
    for i=1:length(hcaSessionStruct.rawBarcodes)
        
    	[f,s] =max(hcaSessionStruct.comparisonStructure.xcorr{i});
        [mV,id]=max(f);
        
         [ b, ix ] = sort( f(:), 'descend' );
         indx = b(1:3)' ;
         resultStruc.maxcoef = [resultStruc.maxcoef;indx];
       %  resultStruc.pos = [resultStruc.pos; ix(1:3)'];
         resultStruc.or = [resultStruc.or; s(ix(1:3)')];
         
         if s(ix(1:3)') == 1
         	resultStruc.pos = [resultStruc.pos; ix(1:3)'];
         else
         	resultStruc.pos = [resultStruc.pos; ix(1:3)'-length(hcaSessionStruct.comparisonStructure.StretchedBar{i})];
         end
         
         if hcaSessionStruct.filterSettings.filter==1 
             [f,s] =max(hcaSessionStruct.comparisonStructure2.xcorr{i});
             [mV,id]=max(f);
        
             [ b, ix ] = sort( f(:), 'descend' );
             indx = b(1:3)' ;
             resultStruc2.maxcoef = [resultStruc2.maxcoef;indx];
            % resultStruc2.pos = [resultStruc2.pos; ix(1:3)'];
             resultStruc2.or = [resultStruc2.or; s(ix(1:3))];
             if s(ix(1:3)') == 1
                resultStruc2.pos = [resultStruc2.pos; ix(1:3)'];
             else
                resultStruc2.pos = [resultStruc2.pos; ix(1:3)'-length(hcaSessionStruct.comparisonStructure.StretchedBar{i})];
             end
             
         end
         maxCorCoefs =[maxCorCoefs; mV id s(id)];
    end
    
    fig1 = figure;
    subplot(2,2,1)

%     import CA.CombAuc.Core.Comparison.compute_p_value;
%     pvals = compute_p_value(resultStruc.maxcoef, hcaSessionStruct.comparisonStructure.evdPar, 'exact' );
 
    p = plot(resultStruc.maxcoef(1:length(hcaSessionStruct.barcodeGen),:),1:length(hcaSessionStruct.barcodeGen),'ob');
    
    p(1).Marker = markers(1);
    p(2).Marker = markers(2);
    p(3).Marker = markers(3);

  %  plot(maxCorCoefs(:,1),1:length(hcaSessionStruct.comparisonStructure),'*')
    hold on
    if hcaSessionStruct.filterSettings.filter==1
    	p2 = plot(resultStruc2.maxcoef(1:length(hcaSessionStruct.barcodeGen),:),1:length(hcaSessionStruct.barcodeGen),'or');
        p2(1).Marker = markers(4);
        p2(2).Marker = markers(5);
        p2(3).Marker = markers(6);
    end
    
    if length(hcaSessionStruct.barcodeGen) < length(hcaSessionStruct.comparisonStructure.xcorr)
        plot([0.1:0.1:1], 0.5+repmat(len1,10,1))
        p3 = plot(resultStruc.maxcoef(end,:),length(hcaSessionStruct.barcodeGen)+1,'ob');
        p3(1).Marker = markers(1);
        p3(2).Marker = markers(2);
        p3(3).Marker = markers(3);
        p4 = plot(resultStruc2.maxcoef(end,:),length(hcaSessionStruct.barcodeGen)+1,'or');
        p4(1).Marker = markers(4);
        p4(2).Marker = markers(5);
        p4(3).Marker = markers(6);
% 
%         import CBT.Hca.UI.gen_hca_res_str;
%         resultStrucC = gen_hca_res_str(len1,hcaSessionStruct.comparisonStructure.xcorr,hcaSessionStruct.comparisonStructure.StretchedBar,markers);
    end
    
   % plot(repmat(mean(maxCorCoefs(:,1)),1,length(hcaSessionStruct.comparisonStructure)),1:length(hcaSessionStruct.comparisonStructure),'r')
   % legend({'barcodes max match score','mean match score'},'Location','sw')
    ylabel('Barcode nr.','Interpreter','latex')
    xlabel('max match score','Interpreter','latex')
    xlim([0.5 1])
    ylim([0,len1+2])
    legend({'$\hat C$','$C_2$','$C_3$','filtered $\hat C$','filtered $C_2$','filtered $C_3$','consensus line'},'Location','sw','Interpreter','latex')

    subplot(2,2,2),
    p3 = plot(resultStruc.pos(1:len1,:),1:len1,'ob');
    p3(1).Marker = markers(1);
    p3(2).Marker = markers(2);
    p3(3).Marker = markers(3);
    hold on
    
    if hcaSessionStruct.filterSettings.filter==1
    	p4 = plot(resultStruc2.pos(1:len1,:),1:len1,'or');
        p4(1).Marker = markers(4);
        p4(2).Marker = markers(5);
        p4(3).Marker = markers(6);
    end
    
    if length(hcaSessionStruct.comparisonStructure.xcorr) > len1
        plot([0:100:100000], 0.5+repmat(len1,length(0:100:100000),1))
        
        p3 = plot(resultStruc.pos(end,:),len1+1,'ob');
        p3(1).Marker = markers(1);
        p3(2).Marker = markers(2);
        p3(3).Marker = markers(3);

        p4 = plot(resultStruc2.pos(end,:),len1+1,'or');
        p4(1).Marker = markers(4);
        p4(2).Marker = markers(5);
        p4(3).Marker = markers(6);
    end
    
    ylabel('Barcode nr.','Interpreter','latex')
    xlabel('best position (pixel)','Interpreter','latex')
    legend({'$\hat C$','$C_2$','$C_3$','filtered $\hat C$','filtered $C_2$','filtered $C_3$'},'Location','sw','Interpreter','latex')
    ylim([0,len1+2])

    

    [dd,ii] =max(resultStruc.maxcoef(:,1));
    
%     ii = 4;
%     [dd,~] =max(resultStruc.maxcoef(ii,1));

    bar = hcaSessionStruct.theoryGen.theoryBarcodes{1};
    %b1=hcaSessionStruct.rawBarcodes{ii};
    b1=hcaSessionStruct.comparisonStructure.StretchedBar{ii};

    if ii <= length(hcaSessionStruct.barcodeGen)
        b1=b1(hcaSessionStruct.comparisonStructure.StretchedBitmask{ii});


        if resultStruc.or(ii,1) == 2
            b1=fliplr(b1);
            % if there is no shift?
            shift=find(fliplr(hcaSessionStruct.comparisonStructure.StretchedBitmask{ii})==1,1);
        else
            shift=find(hcaSessionStruct.comparisonStructure.StretchedBitmask{ii}==1,1);
        end


        cutB = bar(shift+resultStruc.pos(ii,1)-1:shift+resultStruc.pos(ii,1)+length(b1)-2);
        m1=mean(cutB);
        s1= std(cutB);

        subplot(2,2,3);

        plot(bar)
        hold on
        plot([resultStruc.pos(ii,1)+shift:resultStruc.pos(ii,1)+shift+length(b1)-1],zscore(b1)*s1+m1)
        import CBT.Hca.Export.consistency_check;
        display('Running consistency check for the comparison between theory and exp plot 1...')
        consistency_check(b1,cutB,dd);
    else
        posInd = find(hcaSessionStruct.comparisonStructure.StretchedBitmask{ii});
        bitm=hcaSessionStruct.comparisonStructure.StretchedBitmask{ii};
        if resultStruc.or(ii,1) == 2
            b1=fliplr(b1);
            posInd = fliplr(posInd);
            bitm=fliplr(bitm);
        end
    	cutB = bar(resultStruc.pos(ii,1):resultStruc.pos(ii,1)+length(b1)-1);
        m1=mean(cutB);
        s1= std(cutB);

        subplot(2,2,3);

        plot(bar)
        hold on
        plot(resultStruc.pos(ii,1)+posInd,zscore(b1(posInd))*s1+m1)
        
        import CBT.Hca.Export.consistency_check;
        display('Running consistency check for the comparison between theory and exp plot 1...')
        consistency_check(b1,cutB,dd,bitm);
    end
    
    xlabel('pixel nr.')
    ylabel('Rescaled intensity')
    legend({'Theory barcode', strcat([hcaSessionStruct.displayNames{ii}(1:min(20,length(hcaSessionStruct.displayNames{ii}))),'...'])})
    xlim([resultStruc.pos(ii,1)-400 resultStruc.pos(ii,1)+400 ])
    
     % Consistency check that the cc for the cut out barcodes are the same

  
     
    if hcaSessionStruct.filterSettings.filter==1
    	subplot(2,2,4);
        title('Filtered one frame barcode');


    [dd,ii] =max(resultStruc2.maxcoef(:,1));
    
    bar = hcaSessionStruct.theoryGen.theoryBarcodes{1};
    %b1=hcaSessionStruct.rawBarcodes{ii};
    b1=hcaSessionStruct.comparisonStructure2.StretchedBar{ii};
    
    % support for older session files loading
    if ~isfield(hcaSessionStruct.comparisonStructure2, 'StretchedBitmask')
    	hcaSessionStruct.comparisonStructure2.StretchedBitmask = hcaSessionStruct.comparisonStructure.StretchedBitmask;
    end
    
    b1=b1(hcaSessionStruct.comparisonStructure2.StretchedBitmask{ii});
    if resultStruc2.or(ii,1) == 2
        b1=fliplr(b1);
        shift=find(fliplr(hcaSessionStruct.comparisonStructure2.StretchedBitmask{ii})==1,1);

       % resultStruc2.pos(ii,1) = resultStruc2.pos(ii,1)-length(hcaSessionStruct.comparisonStructure2.StretchedBar{ii});
    else
        shift=find(hcaSessionStruct.comparisonStructure2.StretchedBitmask{ii}==1,1);

    end

    cutB = bar(shift+resultStruc2.pos(ii,1)-1:shift+resultStruc2.pos(ii,1)+length(b1)-2);
    m1=mean(cutB);
    s1= std(cutB);
    
    plot(bar)
    hold on
    plot([resultStruc2.pos(ii,1)+shift:resultStruc2.pos(ii,1)+shift+length(b1)-1],zscore(b1)*s1+m1)
    xlabel('pixel nr.')
    ylabel('rescaled intensity')
    legend({'Theory barcode', strcat([hcaSessionStruct.displayNamesFiltered{ii}(1:min(20,length(hcaSessionStruct.displayNamesFiltered{ii}))),'...'])})
    xlim([resultStruc2.pos(ii,1)-400 resultStruc2.pos(ii,1)+400 ])
     
    display('Running consistency check for the comparison between theory and exp plot 2...')
     import CBT.Hca.Export.consistency_check;
     consistency_check(b1,cutB,dd);
%    legend({'Theory barcode','filtered one frame experiment'})

    end
    
    hcaSessionStruct.resultStruc = resultStruc;
    hcaSessionStruct.resultStruc2 = resultStruc2;
    assignin('base','hcaSessionStruct',hcaSessionStruct)
    
    cache('hcaSessionStruct') = hcaSessionStruct ;

    
    
%     subplot(2,2,4);
% 
%     plot(zscore(cutB))
%     hold on
%     plot(zscore(b1))
%    
%     legend({'Theory sequence','Molecule from experiment'})
%     
%     
    
%     
%     m = load('meanF.mat');
%     meanFFTest = interp1(m.meanFFT,linspace(1,length(m.meanFFT),m.len));
%     
%     % THESE SHOULD BE TAKEN FROM CBT.ini/ move them to comparison_settings.
%     import CBT.get_default_barcode_gen_settings;
%     barcodeGenSettings = get_default_barcode_gen_settings();
%     barcodeGenSettings = defaultBarcodeGenSettings;
%     barcodeGenSettings.pixelWidth_nm = 130;
%     barcodeGenSettings.psfSigmaWidth_nm=300;
%     barcodeGenSettings.isLinearTF = 0;
%     barcodeGenSettings.concNetropsin_molar = 6e-6; % Netropsin concentration, units molar
%     barcodeGenSettings.concYOYO1_molar = 2e-8; % YOYO-1 concentration, units molar
%    
%     barcodeGenSettings.meanBpExt_nm = 0.225;
% 

  
    
% 
% 
%     import CA.CombAuc.Core.Zeromodel.generate_random_sequences;
%     import CA.CombAuc.Core.Comparison.compute_correlation;
%     import CA.CombAuc.Core.Comparison.generate_evd_par;
%     import CA.CombAuc.Core.Comparison.cc_fft;
%     import CA.CombAuc.Core.Comparison.compute_p_value;
% 
%     %contigItems/aka fragments
%     evdPar = cell(1,length(hcaSessionStruct.comparisonStructure));
%     pValueMatrix = cell(1,length(hcaSessionStruct.comparisonStructure));
%    
%     % only compute for the first contig/molecule since it takes some time..
%     % later let user choose which ones..
%     for contigNum=1   
% 
%         %should we save randomsequences to resutls? since this has to be
%         %reproducible. but don't want too much things in the result
%         %structure..
%         [ randomSequences ] = generate_random_sequences(2*length(hcaSessionStruct.rawBarcodes{contigNum}),1000,meanFFTest, barcodeGenSettings.psfSigmaWidth_nm/(barcodeGenSettings.pixelWidth_nm/barcodeGenSettings.meanBpExt_nm ),'phase');
%       %  expBit = ones(1,length(hcaSessionStruct.rawBarcodes{contigNum}));
%         ccMax = zeros(1,length(randomSequences));
%         for i=1:length(randomSequences)
%             randomSequences{i} = zscore(randomSequences{i}(1:end/2));
%             
%             [cc1,cc2] =CA.CombAuc.Core.Comparison.cc_fft(randomSequences{i}, hcaSessionStruct.theoryGen.theoryBarcodes{1}');
%             ccMax(i) = max([cc1(:);cc2(:)]);
%             % the comp with bitweights is 20 times slower for this...
% %             import SignalRegistration.XcorrAlign.get_no_crop_lin_circ_xcorrs;
% %             [xcorrs, coverageLens, firstOffset] = get_no_crop_lin_circ_xcorrs(randomSequences{i}, hcaSessionStruct.theoryGen.theoryBarcodes{1}, expBit,hcaSessionStruct.theoryGen.bitmask{1});
%         end
%        	[~,rSq,evdPar{contigNum}] = generate_evd_par( ccMax,[],length(randomSequences{1})/5, 'exact2' );
%       	%[cc1,cc2] = cc_fft(zscore(contigItems{contigNum}.barcode),refBarcode);
% 
%         %ccValueMatrix(contigNum,:) = [cc1 cc2];
%         import CBT.ExpComparison.Core.exact_full_PDF;
%         
%         pValueMatrix{contigNum} = compute_p_value(hcaSessionStruct.comparisonStructure{2}, evdPar{contigNum}, 'exact' );
% 
%     end
    
end