function [cache] = display_additional_results_theory_ui(ts,cache)
    if nargin < 2
        cache = containers.Map();
    end
    
    hcaSessionStruct = cache('hcaSessionStruct') ;
     
	tabTitle = 'Additional results II';
    [hTabTheoryImport] = ts.create_tab(tabTitle);
    hPanelTheoryImport = uipanel(hTabTheoryImport);
    ts.select_tab(hTabTheoryImport);

	import Fancy.UI.FancyList.FancyListMgr;
    lm = FancyListMgr();
    lm.set_ui_parent(hPanelTheoryImport);
    lm.make_ui_items_listbox();
%      
%     % Change to filtered and unfiltered names..
%     lm.add_list_items(hcaSessionStruct.names, hcaSessionStruct.names);
%     
    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    
    flmbs1 = FancyListMgrBtnSet();
    flmbs1.NUM_BUTTON_COLS = 2;
    flmbs1.add_button(add_theory(ts,'add theory'));
    flmbs1.add_button(add_theory(ts,'add fragment'));
    
    
        import Fancy.UI.FancyList.FancyListMgrBtnSet;
    
    flmbs3 = FancyListMgrBtnSet();
    flmbs3.NUM_BUTTON_COLS = 1;
    flmbs3.add_button(make_remove_consensus_btn());   
    
    flmbs2 = FancyListMgrBtnSet();
    flmbs2.NUM_BUTTON_COLS = 2;
	flmbs2.add_button(FancyListMgr.make_select_all_button_template());
    flmbs2.add_button(FancyListMgr.make_deselect_all_button_template());
    
    
    flmbs4 = FancyListMgrBtnSet();
    flmbs4.NUM_BUTTON_COLS = 1;
    flmbs4.add_button(theory_vs_theory(ts));

       function [btnAddKymos] =add_theory(ts,addwhat)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnAddKymos = FancyListMgrBtn(...
            addwhat, ...
            @(~, ~, lm) on_add_theory(lm, ts));
        
        function [] = on_add_theory(lm, ts)
            
             [barcodeFilenames{1}, barcodeFilenamesDirpath{1}] = uigetfile(...
                    {'*.mat;'}, ...
                    'Select fragment theory barcodes to import', ...
                    pwd, ...
                    'MultiSelect','on');
             lm.add_list_items(barcodeFilenames,barcodeFilenamesDirpath);
        end
	end
       
%     
       % add barcodes
    function [btnAddKymos] =theory_vs_theory(ts)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnAddKymos = FancyListMgrBtn(...
            'Plot theory vs theory comparison(s)', ...
            @(~, ~, lm) on_theory_vs_theory(lm, ts));
        
        function [] = on_theory_vs_theory(lm, ts)
            [selected, ~] = get_selected_list_items(lm);
            barcodeFilenames =  selected(:,1);
            barcodeFilenamesDirpath = selected(:,2);
         
             barcodeData = load(strcat(barcodeFilenamesDirpath{1},barcodeFilenames{1}));
             barcode = barcodeData.hcaSessionStruct.theoryBarcodes{1};
             bitmask = ones(1,length(barcode));


             barcodeData = load(strcat(barcodeFilenamesDirpath{2},barcodeFilenames{2}));
             theory = barcodeData.hcaSessionStruct.theoryBarcodes{1};
             theoryBitmask = ones(1,length(theory));
             
             barcode = barcode(logical(bitmask));
            import SignalRegistration.XcorrAlign.get_no_crop_lin_circ_xcorrs;
            [xcorrs, ~, ~] = get_no_crop_lin_circ_xcorrs(barcode,theory, bitmask, theoryBitmask);
%            
            [f,s] =max(xcorrs);
        
            [ b, ix ] = sort( f(:), 'descend' );
        	indx = b(1:3)' ;
            maxcoef = [indx];

            or = [s(ix(1:3)')];
         
             if s(ix(1:3)') == 1
                pos = [ix(1:3)'];
             else
                pos = [ix(1:3)'-length(barcode)];
             end
           
            if or(1) == 2
                barcode=fliplr(barcode);
                % if there is no shift?
                shift=find(fliplr(bitmask)==1,1);
            else
                shift=find(bitmask==1,1);
            end

        disp(strcat(['Best position is starting at ' num2str(pos) ' along the chromosome']));
        cutB = theory(shift+pos(1)-1:shift+pos(1)+length(barcode)-2);
        m1=mean(cutB);
        s1= std(cutB);

        figure
        plot(theory)
        hold on
        plot([pos(1)+shift-1:pos(1)+shift+length(barcode)-2],zscore(barcode)*s1+m1)
        import CBT.Hca.Export.consistency_check;
        display('Running consistency check for the comparison between theory and exp plot 1...')
        consistency_check(barcode,cutB,maxcoef(1));  

        xlabel('pixel nr.')
        ylabel('Rescaled intensity')
        legend({ 'Theory','Fragment'})
        xlim([pos(1,1)-400 pos(1,1)+400 ])
        end
    end
       
    function [btnRemoveConsensus] = make_remove_consensus_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnRemoveConsensus = FancyListMgrBtn(...
            'Remove selected sequence(s)', ...
            @(~, ~, lm) on_remove_selected_consensus(lm));
        function [] = on_remove_selected_consensus(lm)
            lm.remove_selected_items();
        end
    end

    lm.add_button_sets(flmbs1,flmbs2,flmbs3,flmbs4);
  
    cache('hcaSessionStruct') = hcaSessionStruct;
end