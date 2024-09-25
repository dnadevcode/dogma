function [compI,rezI,m,st,sS] = run_comparison(barGen,fastaFile,parlistcell,parchangeId,par,sF,sets)

import Helper.get_theory_twostate;

% parlistcell
m = zeros(1,length(par));
st = zeros(1,length(par));
sS = zeros(1,length(par));
compI = cell(1,length(par));
rezI = cell(1,length(par));

parfor k=1:length(par)
%     k
%     parlistcell{parchangeId} = par(k);
    pars = [parlistcell(1:parchangeId-1) par(k) parlistcell(parchangeId+1:end)]; 
    [~,theoryStr] = get_theory_twostate(fastaFile,pars{:});

    % todo: change to hca_compare_distance
    [compI{k},rezI{k},~] = compare_to_t(barGen,theoryStr,sF,sets);

    m(k) =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI{k});
    st(k) = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI{k});

%     m(k)
    stoufferScoresSimple =cellfun(@(x) double(norminv(1-x.pval)),compI{k},'un',true);
    sS(k ) = mean(stoufferScoresSimple);%(stoufferScoresSimple>1.6))
end

end

