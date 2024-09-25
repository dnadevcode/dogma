% compare to len_based_null


% Generate randomized data
import Nullmodel.gen_random;

NUM_RAND_FRAGMENTS = 1;
PSF_WIDTH_PIXELS = 300/110;
RAND_LENGTH_2 = 4000;
minOverlapLen = RAND_LENGTH_MIN;
% overlap length in general will be fixed for the analysis.
RAND_LENGTH_MIN = 600;

% sF = 1;%
sF = 1;
gap = 100;

%  calculates all MP and MPI
NN = 30;
out ='output';

TT = 10;

nufit  =zeros(TT,NN);
nufitGaussian = zeros(TT,NN);
for kk=1:TT
    tic
    import Nullmodel.gen_scores_pcc;
    [maxPCC,allPCC] = gen_scores_pcc(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN+(kk-1)*gap,RAND_LENGTH_2,NN,sF,out,gap);
    toc
    
    
    import Nullmodel.full_dist_mle;
    % compare with actual data
    for ii=1:NN
    
        scores = cell2mat(allPCC{ii}{1}); % mat = cell2mat(cellfun(@(x) reshape(x, 1, []), allPCC{1}, 'UniformOutput', false));
    
        f = @(x0) abs(full_dist_mle(scores(:),0,x0));
    
        [nufit(kk,ii),fval,exitflag,output] = fminbnd(f, 2, RAND_LENGTH_MIN/2,optimset('TolX',1e-10));

        % alternatively Gaussian
        params = mle(scores(:), 'distribution', 'norm');
        nufitGaussian(kk,ii) = 1/(params(2)).^2+2;
    end

end

shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
longL = RAND_LENGTH_2+(0:NN-1)*gap;


f=figure
tiledlayout(2,2)
nexttile
imagesc(longL,shortL,nufit )
ylabel('Length short')
xlabel('Length long')
colorbar
title('Fit PCC distribution')
nexttile
imagesc(longL,shortL,nufitGaussian )
title('Fit Gaussian distribution')
ylabel('Length short')
xlabel('Length long')
% figure,plot(nufit)
nexttile
plot(shortL, mean(nufit,2))
hold on
errorbar(shortL, mean(nufit,2),std(nufit,1,2))
plot(shortL, mean(nufitGaussian,2))
hold on
errorbar(shortL, mean(nufitGaussian,2),std(nufitGaussian,1,2))
lgd = legend({'Fit PCC','Fit Gaussian'})
lgd.Location = 'eastoutside'
xlabel('Length short')
ylabel('Parameter fitted value ')

saveas(f,fullfile(foldsave,'figTS.png'))

