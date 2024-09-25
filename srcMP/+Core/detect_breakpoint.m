function [subThry1,subThry2,exBzscored,pcc,C,b1,bCB,bDLE,chr4Start,chr11Start] = detect_breakpoint(barid,numWorkers,expBar,cs,theoryStruct,theoryStruct2,minLen, timestamp,...
    MIN_OVERLAP_PIXELS,gapW,bpPx)

    %   Returns:
    %       subThry1 - subtheory of Chr 4
    %       subThry2 - subtheory of Chr 11
    %       pcc - pcc scores
    %       C - combined score
    
    % step1: run comparison with averaged re-scaling factor, so that
    % theoriess are without scaling/shifting.

    for i=1:length(theoryStruct)
        thry{i} = importdata(theoryStruct{i}.filename);
    end

        for i=1:length(theoryStruct)
        thry2{i} = importdata(theoryStruct2{i}.filename);
    end
    import Core.compare_mp_all;

%     % here: gen smaller version fo thry?
%     fastaFile = {'/proj/snic2022-5-384/users/x_albdv/data/chr4.fna','/proj/snic2022-5-384/users/x_albdv/data/chr11.fna'}
% 
% 
%     gapW = 10000;
%     limT1save = round([limitsT1(1)+gapW-lenB+1 limitsT1(1)+gapW]*theoryStruct{1}.pixelWidth_nm/theoryStruct{1}.meanBpExt_nm);
%     limT2save = [limitsT2(2)-gapW limitsT2(2)-gapW+lenB-1];
% 
%     delete('test1.fasta','test2.fasta');
%     fastaF = fastaread(fastaFile{1});
%     fastaF2 = fastaread(fastaFile{1});
% 
%         fastawrite('test1.fasta',fastaF.Sequence(limT1save(1):limT1save(2)));
% 
%     fastawrite('test2.fasta',fastaF2.Sequence(limT2save(1):limT2save(2)));
%     fastaFile2 = {'test1.fasta','test2.fasta'};
% 
%     [tS{1},tG{1}] = gen_theor(fastaFile2,[]);
%     
%     [tS{2},tG{2}] = gen_theor(fastaFile2,2);



    bar = [];
    bar{1}.rawBarcode =  expBar{3};
    bar{1}.rawBitmask = ones(1,length(expBar{3}));

    meanRF = (cs{1}.compStr{barid}.bestBarStretch+cs{2}.compStr{barid}.bestBarStretch)/2;

    compStr = [];
    [mpI,mp,~,~,compStr{1}] = compare_mp_all(theoryStruct,bar,minLen,1, timestamp,meanRF,MIN_OVERLAP_PIXELS,numWorkers);
    [mpI,mp,~,~,compStr{2}] = compare_mp_all(theoryStruct,bar,minLen,2, timestamp,meanRF,MIN_OVERLAP_PIXELS,numWorkers);

%     import Core.plot_pair_simple;     
%     fig=figure
%     tiledlayout(2,2,'TileSpacing','tight','Padding','tight')
%     nexttile
%     ix=1;CH=1;
%     plot_pair_simple(compStr{ix}{1}.pA,compStr{ix}{1}.pB,compStr{ix}{1}.bestBarStretch,compStr{ix}{1}.or(1),...
%     bar{1}.rawBarcode,bar{1}.rawBitmask,thry{ix},MIN_OVERLAP_PIXELS,fig,[]);
%     nexttile
%     ix=1;CH=3;
%     plot_pair_simple(compStr{ix}{1}.pA,compStr{ix}{1}.pB,compStr{ix}{1}.bestBarStretch,compStr{ix}{1}.or(1),...
%     bar{1}.rawBarcode,bar{1}.rawBitmask,thry2{ix},MIN_OVERLAP_PIXELS,fig,[]);
%       nexttile
%     ix=2;CH=1;
%     plot_pair_simple(compStr{ix}{1}.pA,compStr{ix}{1}.pB,compStr{ix}{1}.bestBarStretch,compStr{ix}{1}.or(1),...
%     bar{1}.rawBarcode,bar{1}.rawBitmask,thry2{ix},MIN_OVERLAP_PIXELS,fig,[]);
%     nexttile
%     ix=2;CH=3;
%     plot_pair_simple(compStr{ix}{1}.pA,compStr{ix}{1}.pB,compStr{ix}{1}.bestBarStretch,compStr{ix}{1}.or(1),...
%     bar{1}.rawBarcode,bar{1}.rawBitmask,thry{ix},MIN_OVERLAP_PIXELS,fig,[]);
%     saveas(fig,'figs/fig8.png')

%%
    ix = 1;
    limitsE1 = [compStr{ix}{1}.pA compStr{ix}{1}.pA+MIN_OVERLAP_PIXELS];
    limitsT1 = [compStr{ix}{1}.pB compStr{ix}{1}.pB+MIN_OVERLAP_PIXELS];
    % limitsT1 = [cs{ix}.compStr{barid}.pA cs{ix}.compStr{barid}.pA+MIN_OVERLAP_PIXELS];
    ix=2;
    limitsE2 = [compStr{ix}{1}.pA compStr{ix}{1}.pA+MIN_OVERLAP_PIXELS];
    limitsT2 = [compStr{ix}{1}.pB compStr{ix}{1}.pB+MIN_OVERLAP_PIXELS];

% gapW = 100;
if limitsE1(1)<limitsE2(1)
    theoryLeftIdx = 1;
    limitsE = [limitsE1(2)-gapW limitsE2(1)+gapW];
    lenB = limitsE(2)-limitsE(1)+1;
    limT1 = [limitsT1(2)-gapW limitsT1(2)-gapW+lenB-1];
    limT2 = [limitsT2(1)+gapW-lenB+1 limitsT2(1)+gapW];
else
    theoryLeftIdx = 2;
    limitsE = [limitsE2(2)-gapW limitsE1(1)+gapW];
    lenB = limitsE(2)-limitsE(1)+1;
    limT1 = [limitsT1(1)+gapW-lenB+1 limitsT1(1)+gapW];
    limT2 = [limitsT2(2)-gapW limitsT2(2)-gapW+lenB-1];
end

subBar=[];subThry2=[];subThry1=[];
channels = [1 3];
for idx=channels
    expB =     imresize( expBar{idx},'Scale' ,[1 meanRF]);
    
    subBar{idx} =  expB(limitsE(1):limitsE(2));
    if idx==3
        subThry1{idx} = thry{1}(limT1(1):limT1(2));
        subThry2{idx} = thry{2}(limT2(1):limT2(2));
    else
        subThry1{idx} = thry2{1}(limT1(1):limT1(2));
        subThry2{idx} = thry2{2}(limT2(1):limT2(2));
    end
end

% idx = 1;
% figure,plot(zscore(subThry1{idx})+3,'blue')
% hold on
% plot(zscore(subThry2{idx})+6,'red')
% plot(zscore(subBar{idx})+9,'black')
% legend({'THR1','THR2','EXP'})

% step2: PCC for all modified versions

% expB =     imresize( expBar{idx},'Scale' ,[1 meanRF]);

pcc = [];
for idx=channels
    exBzscored{idx} = zscore(subBar{idx},1);
    pcc{idx} = [];
    for j = 1:length(subThry1{idx})
        if theoryLeftIdx == 1
            barTest = [subThry1{idx}(1:j) subThry2{idx}(j+1:end)];
        else
            barTest = [subThry2{idx}(1:j) subThry1{idx}(j+1:end)];
        end
        pcc{idx}(j) = zscore(barTest,1)*exBzscored{idx}'/length(barTest);
    end
end
%%
idx=1;
C =(pcc{1}+pcc{3})/2;
[a1,b1] = max(C);
[aCB,bCB] = max(pcc{1});
[aDLE,bDLE] = max(pcc{3});

ebar = zscore(subBar{idx})+9;

fig = figure

xax = [1:length(ebar)]*bpPx/10^6;
tiledlayout(3,1,'TileSpacing','none','Padding','none')
for idx=[1 3]
    nexttile
    % idx = 1;
    ebar = zscore(subBar{idx})+9;

    plot(xax,zscore(subThry1{idx})+3,'blue')
    hold on
    plot(xax,zscore(subThry2{idx})+6,'red')
    if theoryLeftIdx==1
        plot(xax(1:b1),ebar(1:b1),'blue-.')
        plot(xax(b1+1:length(ebar)),ebar(b1+1:end),'red-.')
    else
        plot(xax(1:b1),ebar(1:b1),'red-.')
        plot(xax(b1+1:length(ebar)),ebar(b1+1:end),'blue-.')
    end
end
xlabel('Position (Mbp)')
% legend({strcat(['Theoretical Chr.4 from ' num2str(limT1(1)*bpPx/10^6) 'Mbp']),strcat(['Theoretical Chr.11, from ' num2str(limT2(1)*bpPx/10^6) ' Mbp']),'Experimental barcode'},'Location','southeastoutside')
nexttile
plot(xax,pcc{1})
hold on
plot(xax,pcc{3})
plot(xax,C,'black')
plot(xax(b1),a1,'redx')
legend({'Score CB','Score DLE-1','Combined score','Breakpoint'},'Location','southoutside')
    saveas(fig,strcat(['figs/fig11bar' num2str(barid) '.png']))

    chr4Start = limT1(1)*bpPx/10^6;
        chr11Start = limT2(1)*bpPx/10^6;




