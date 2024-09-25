function  [t,cthreshFun,params,sol,lengw,lengL] = pval_par_simple(lengw,lengL, w,numBars,islinear,kernelSigma)

    % Args:
    %   low, high, w, numBars, sets
    
    % [cthreshFun,params,sol,leng] = pval_par_simple(100:10:200,50,100,1,[1 0],2.3)
    
    % compare to import CBT.Hca.Core.Pvalue.precompute_pvalue_files_multi_theory;
    import CBT.Hca.Core.Pvalue.gensv;

    % how to generate: still need to test correctness..
%     [t,cthreshFun,params,sol,lengw,lengL] = pval_par_simple(100:100:200,200:100:500,50,1000,[1 0],2.3)
    sets.kernelsigma = kernelSigma;
    %% structural variations list / probably does not matter too much if it's circular or not for (0)
    sets.svList = zeros(1,numBars);

    % import mp.mp_profile_stomp_dna;
    kk = 2^16;

    % import mp.mp_full_masked_profile_stomp_dna;
    comparisonFun = @(x,y,z1,z2,w,u) massfull(x,y,z1,z2,w,kk,u);
    tic
    
%     comparisonFun = @(x,y,z,w) unmasked_MASS_PCC(y,x,z,2^(4+nextpow2(length(x))));

%  stretch = sets.stretch;
% import Rand.generate_linear_sv;

sol = zeros(length(lengw),length(lengL));
params = cell(length(lengw),length(lengL));
par1 = zeros(1,length(lengL));
par2 = zeros(1,length(lengL));

for j = 1:length(lengw)
    j


    lenS =lengw(j);
    for k=1:length(lengL)
        [bar1,bar2,~,~]  = arrayfun(@(x) gensv(lenS, lengL(k), x,kernelSigma,islinear),sets.svList,'UniformOutput',false);
        ccMax = zeros(1,length(bar1));

        for idx=1:length(bar1)
            [profile]=  comparisonFun(bar1{idx}', bar2{idx}', ones(length(bar1{idx}),1), ones(length(bar2{idx}),1), lenS,islinear);
            ccMax(idx) = max(profile(:));
        end

        %         import Pvalue.compute_evd_params; % x0 should be better bound, make an analysis
        pars = pvalpar(ccMax(:),100); % why did I set x0 to 100 here?
        p = @(x) 1-(0.5+0.5*(1-betainc((x).^2,0.5, pars(1)/2-1,'upper'))).^pars(2);
        % we want to have approx region
        sol(j,k) = fzero(@(x) p(x)-0.01, [0 1]);
        par1(k) = pars(1);
        par2(k) = pars(2);

        % generate barcodes of long length
%         parfor idx =1:length(bar1)

%             maxL = zeros(1,length(stretch));
%             for l=1:length(stretch)
%                 bar1stretch = imresize(bar1{idx},[1 round(lengw(j)*stretch(l))]);
%             maxL(l) = max(mp(:));
%             end
    



    end
    for k = 1:length(lengL)
        params{j,k} =[ par1(k) par2(k)];
    end;

 
end

sol2 = sol;

% sol2 = sol2 + transpose(triu(sol2))-tril(sol2);

cthreshFun = @(x,y) interp2(lengw, lengL, sol2', x,y);

par1 = cellfun(@(x) x(1),params);
par2 = cellfun(@(x) x(2),params);
par1Fun = @(x,y) interp2(lengw, lengL, par1', x,y);
par2Fun = @(x,y) interp2(lengw, lengL, par2', x,y);

    

save(strcat(['pval_' num2str(w) '.mat']),'cthreshFun','params','sol','par1Fun','par2Fun','lengw','lengL','-v7.3');
%    
% 
% % % 
% f = figure;imagesc(lengL,lengw,sol);colorbar
% xlabel("Long length")
% ylabel("Short length")
% title("dist thresh for p-value'thresh 0.01")
% saveas(f,'pvalthresh','epsc')
% 
% f = figure;imagesc(lengL,lengw,par1);colorbar
% xlabel("Long length")
% ylabel("Short length")
% title("par1")
% saveas(f,'par1','epsc')
% 
% par2(par2>10^10) = nan;
% f = figure;imagesc(lengL,lengw,par2);colorbar
% xlabel("Long length")
% ylabel("Short length")
% title("par2")
% saveas(f,'par2','epsc')


t= toc;
end

