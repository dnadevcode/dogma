function  [fullmat] = pval_par_full(length1,length2,numBars,islinear,kernelSigma,name,w)
    % COMPUTES PVALUES dist FOR the most significant fragment    

    %
    % Args:
    %  length1:  length of barcode 1
    %  length2: length of barcode 2
    %  w: window width for running mp
    import CBT.Hca.Core.Pvalue.gensv;
    import mp.mpfull;
    % pval_par_full(300:50:600,300:50:600,100,[0,0],2.88)
%      pval_par_full(300,300,100,[0,0],2.88)
    % import mp.mp_profile_stomp_dna;
    kk = 2^16;
    svList = zeros(1,numBars);
    % import mp.mp_full_masked_profile_stomp_dna;
    comparisonFun = @(x,y,z1,z2,w,u) mpfull(x,y,z1,z2,w,kk,u);
    
    if nargin < 7
        w = 40:10:min(length1(end),length2(end));
    end
    
    fullmat = cell(length(length1),length(length1),length(w));
    
    for i = 1:length(length1)
        bit1=  ones(length1(i),1);
        l1 = length1(i);
        i
        for j = 1:length(length2)
            l2 = length2(j);
            bit2 = ones(length2(j),1);
            [bar1,bar2,~,~]  = arrayfun(@(x) gensv(l1,l2, x,kernelSigma,islinear),svList,'UniformOutput',false);
%             pvec= cell(1,length(w));
            par1= nan(1,length(w));
            par2= nan(1,length(w));

            parfor k= 1:length(w)
                if w(k) < min(l1,l2)
                    vals = zeros(1,numBars);

                    for idx=1:numBars
                            [mp,~]=  comparisonFun(bar1{idx}', bar2{idx}',bit1, bit2, w(k),islinear);
                            vals(idx) =   max(mp)  ;                  
    %                         fullmat(i,j,k) = max(maxL);
                    end
                    pars = pvalpar(vals(:),w(k)/3); % why did I set x0 to 100 here?
                    par1(k) = pars(1);
                    par2(k) = pars(2);
%                     pvec2(k,) = pvalpar(vals(:),w(k)/3); 
                end
            end
            for k = 1:length(w)
                fullmat{i,j,k} =[ par1(k) par2(k)];
            end;
        end
    end
    % [cthreshFun,params,sol,leng] = pval_par(100:10:200,50,100,1,[1 0],2.3)
    
    save(strcat(['pval_' name '.mat']),'fullmat','length1','length2','numBars','islinear','kernelSigma','w','-v7.3');

    
% sol2 = sol;
% 
% sol2 = sol2 + transpose(triu(sol2))-tril(sol2);
% 
% cthreshFun = @(x,y) interp2(lengw, lengL, sol2', x,y);
% 
% par1 = cellfun(@(x) x(1),params);
% par2 = cellfun(@(x) x(2),params);
% par1Fun = @(x,y) interp2(lengw, lengL, par1', x,y);
% par2Fun = @(x,y) interp2(lengw, lengL, par2', x,y);
% 
%     
% 
% save(strcat(['pval_' num2str(w) '.mat']),'cthreshFun','params','sol','par1Fun','par2Fun','lengw','lengL','-v7.3');
%    


%     % how to generate: still need to test correctness..
% %     [cthreshFun,params,sol,leng] = pval_par(100:10:200,10,50,1000,1,[1 0],2.3)
% % [t,cthreshFun,params,sol,lengw,lengL] = pval_par_simple(100:10:200,200:100:500,50,1000,[1 0],2.3)
%     %% structural variations list / probably does not matter too much if it's circular or not for (0)
%     sets.svList = zeros(1,numBars);
% 
% 
%     tic
% 
% %  stretch = sets.stretch;
% % import Rand.generate_linear_sv;
% 
% sol = zeros(length(length1),length(length1));
% params = cell(length(length1),length(length1));
% for j = 1:min(numws,length(length1))
%     j
%     for k=j:length(length1)
%         [bar1,bar2,matchTable,lengths]  = arrayfun(@(x) gensv(length1(j), length1(k), x,sets,islinear),sets.svList,'UniformOutput',false);
%         ccMax = zeros(1,length(bar1));
%         for idx =1:length(bar1)
% 
% %             ccMax2 = zeros(1,length(bar1));
%             for idx =1:length(bar1)
%                 maxL = zeros(1,length(stretch));
%                 for l=1:length(stretch)
%                     bar1stretch = imresize(bar1{idx},[1 round(length1(j)*stretch(l))]);
%                     % stretch 
%                     % run this for lengthJ
% 
%                      [mp,mpI]=  comparisonFun(bar1stretch', bar2{idx}', ones(length(bar1stretch),1), ones(length(bar2{idx}),1), length1(j),islinear);
% %                   
%                     maxL(l) = max(mp);
%                 end
%                 ccMax(idx) = max(maxL);
%     
%             end
% 
% 
% %         import Pvalue.compute_evd_params; % x0 should be better bound, make an analysis
%         params{j,k} = pvalpar(ccMax(:),100); % why did I set x0 to 100 here?
%         p = @(x) 1-(0.5+0.5*(1-betainc((x).^2,0.5,params{j,k}(1)/2-1,'upper'))).^params{j,k}(2) ;
%         % we want to have approx region
%         sol(j,k) = fzero(@(x) p(x)-0.01, [0 1]);
%     end
% 
% end
% 
% sol2 = sol;
% 
% sol2 = sol2 + transpose(triu(sol2))-tril(sol2);
% 
% cthreshFun = @(x,y) interp2(length1, length1, sol2, x,y);
% 
% 
% save(strcat(['pval_' num2str(w) '.mat']),'cthreshFun','params','sol','leng','-v7.3');
% %     
% 
% t= toc;
% end
% 
