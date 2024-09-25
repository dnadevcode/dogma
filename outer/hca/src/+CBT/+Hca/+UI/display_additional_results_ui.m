function [cache] = display_additional_results_ui(ts,cache)
    % display_additional_results_ui
    % Displays options for additional results of hca
    %
    %     Args:
    %         ts
    %         cache (container): cached results
    % 
    %     Returns:
    %         cache (container): cached results
    %

% 

% 
%    function [btnAddKymos] = compute_p_plot(ts,title)
%         import Fancy.UI.FancyList.FancyListMgrBtn;
%         btnAddKymos = FancyListMgrBtn(title,@(~, ~, lm) on_compute_p_plot(lm, ts,title));
%         
%         function [] = on_compute_p_plot(lm, ts,title)
%             try 
%                 [file,path] = uigetfile({'*.txt'},'load pre-computed p-value file');
%                 fullPath = strcat([path,file]);  
%                 addpath(genpath(path));
%                 import CBT.Hca.Import.load_pval_struct;
%                 [ pvalData.len1, pvalData.data ] = load_pval_struct(fullPath);
%             catch
%                  disp('No pre-computed p-value database chosen. Running precompute method... ');
%                  import CBT.Hca.Core.Pvalue.pregenerate_pvalue_db;
%                  pregenerate_pvalue_db('pval.txt',pwd);
%             end
%             
% 
%             hcaSessionStruct.sets = sets;
%             import CBT.Hca.Core.Pvalue.compute_p_val;
%             [ hcaSessionStruct ] = compute_p_val(pvalData, hcaSessionStruct );
%          
%             cache('hcaSessionStruct') =hcaSessionStruct ;
%             assignin('base','hcaSessionStruct',hcaSessionStruct)
%         end
%    end
% 
%    function [btnAddKymos] = plot_p_plot(ts,title)
%         import Fancy.UI.FancyList.FancyListMgrBtn;
%         btnAddKymos = FancyListMgrBtn(title, @(~, ~, lm) on_plot_p_plot(lm, ts,title));
%       
%         function [] = on_plot_p_plot(lm, ts,title)
%         pvals = hcaSessionStruct.pValueResults.pValueMatrix;
%         savePvals=1;
%         if savePvals==1
%             defaultMatFilename ='savepvalshere';
%             [~, matDirpath] = uiputfile('*.txt', 'Save pvals data as', defaultMatFilename);
%             CBT.Hca.Export.export_pvals_txt(hcaSessionStruct.pValueResults,pvals,matDirpath)
%         end
%         
%         fig1 = figure;
%         plot(pvals,1:length(pvals),'rx')
%         hold on
% 
%         if sets.filterSettings.filter==1
%             pvalsFiltered = hcaSessionStruct.pValueResults.pValueMatrixFiltered;
%             plot(pvalsFiltered,1:length(pvalsFiltered),'bo')
%         end
%         
%         ylabel('Barcode nr.','Interpreter','latex')
%         xlabel('p-value','Interpreter','latex')
%         legend({'Unfiltered p-value','Filtered p-value'},'Interpreter','latex')
%       
%         cache('hcaSessionStruct') =hcaSessionStruct ;
% 
%         assignin('base','hcaSessionStruct',hcaSessionStruct)
%         end
%    end
% 
%    function [btnAddKymos] = true_positive_statistics(ts,title)
%         import Fancy.UI.FancyList.FancyListMgrBtn;
%         btnAddKymos = FancyListMgrBtn(...
%             title, ...
%             @(~, ~, lm) on_true_positive_statistics(lm, ts,title));
%         
%         function [] = on_true_positive_statistics(lm, ts,title)
%             
%             defaultPlace = 0;
%             defaultError = 20;
%             defaultChromosome = 1;
%            % defaultThresh = 0.01;
%             titleText = 'Correct place';
%             import CBT.Hca.UI.choose_correct_place;
%             [ sets.correctChromosome, sets.correctPlace,sets.allowedError,sets.pvaluethresh] = choose_correct_place(defaultChromosome,defaultPlace,defaultError,sets.pvaluethresh,titleText); 
% 
%             import CBT.Hca.UI.compute_true_positives;
%             hcaSessionStruct = compute_true_positives(hcaSessionStruct, sets);
%             
%             cache('hcaSessionStruct') =hcaSessionStruct ;
% 
%             assignin('base','hcaSessionStruct',hcaSessionStruct)
% 
%         end
%    end
% % 
%     lm.add_button_sets(flmbs1,flmbs2,flmbs3,flmbs31,flmbs4,flmbs5,flmbs7);
%   
%     cache('hcaSessionStruct') = hcaSessionStruct;
end