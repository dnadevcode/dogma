function [] = plot_any_bar(idx,bars,rezMax,theoryStruct,thryIdx, f)
%   plots any barcode against thry

% if nargin < 6
%     f = figure;
% end


    lenBarTested = length(bars{idx}.rawBarcode);
    curBar = interp1(bars{idx}.rawBarcode, linspace(1,lenBarTested,lenBarTested*rezMax{thryIdx}{idx}.bestBarStretch));
    curBit = bars{idx}.rawBitmask(round(linspace(1,lenBarTested,lenBarTested*rezMax{thryIdx}{idx}.bestBarStretch)));
    curBar = curBar(curBit);

% curBar = imresize(bars{idx}.rawBarcode,'Scale' ,[1 comparisonStruct{idx}.bestBarStretch]);
% curBit = imresize(bars{idx}.rawBitmask,'Scale' ,[1 comparisonStruct{idx}.bestBarStretch]);
% curBar = curBar(curBit);
% curBar = imresize(bars{idx}.rawBarcode(bars{idx}.rawBitmask),'Scale' ,[1 comparisonStruct{idx}.bestBarStretch]) ;
% curBit = ones(1,length(curBar));

try % also define stdbar
    stdbar = nanstd(bars{idx}.alignedKymo);
    stdbar = stdbar(bars{idx}.lE:bars{idx}.rE);
    stdbar = interp1(stdbar, linspace(1,lenBarTested,lenBarTested*rezMax{thryIdx}{idx}.bestBarStretch));
    stdbar = stdbar(curBit);
catch
    stdbar = nanstd(bars{idx}.data.alignedKymo);
    stdbar = stdbar(bars{idx}.data.lE:bars{idx}.data.rE);
    stdbar = interp1(stdbar, linspace(1,lenBarTested,lenBarTested*rezMax{thryIdx}{idx}.bestBarStretch));
    stdbar = stdbar(curBit);
end

if rezMax{thryIdx}{idx}.or(1)~=1
    curBar = fliplr(curBar);
    try
       stdbar = fliplr(stdbar);
    end
end

if ~iscell(theoryStruct) % enable possibility
    thr = theoryStruct(thryIdx).rawBarcode;
else
    thr = importdata(theoryStruct{thryIdx}.filename);
end

if ~iscell(theoryStruct) 
    isLinearTF = theoryStruct(thryIdx).isLinearTF;
else
    isLinearTF = theoryStruct{thryIdx}.isLinearTF;
end

if ~isLinearTF
    thr = [thr thr];
else
    thr = [thr nan(1,length(thr))];
end
% thrLambda = importdata(theoryStruct2{1}.filename);

% rawBg = barcodeGen{idx}.rawBg; %mode(cellfun(@(x) x.rawBg,barcodeGenLambda));
% meanlambda = mean(cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGenLambda));
% meanbar = mean(cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGen));
pos = find(curBit==1,1,'first');
% figure,plot(pos+[comparisonStruct{idx}.pos:comparisonStruct{idx}.pos+length(curBar)-1],zscore(curBar))
% hold on
% plot(comparisonStruct{idx}.pos:comparisonStruct{idx}.pos+length(curBar)-1,zscore(thr(comparisonStruct{idx}.pos:comparisonStruct{idx}.pos+length(curBar)-1)))
% legend({'experiment','thry'})
% 
% pccScore = zscore(curBar,1)*zscore(thr(comparisonStruct{idx}.pos:comparisonStruct{idx}.pos+length(curBar)-1),1)'/length(curBar)
% title(['pcc = ' ,num2str(pccScore)])
if nargin <6 || isempty(f)
    f = figure;
end
x = pos+[rezMax{thryIdx}{idx}.pos(1):rezMax{thryIdx}{idx}.pos(1)+length(curBar)-1];

xconf = [x x(end:-1:1)] ;  
y =  zscore(curBar);
% stdbarscore = 1.96*stdbar/sqrt(size(bars{idx}.alignedKymo,1));
stdbarscore = 3*stdbar;

yconf = [y+stdbarscore/std(curBar,1) y(end:-1:1)-stdbarscore/std(curBar,1)];
% 
% figure
p = fill(xconf,yconf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';           
% 
hold on
plot(x,y,'r-')
% hold off
plot(rezMax{thryIdx}{idx}.pos(1):rezMax{thryIdx}{idx}.pos(1)+length(curBar)-1,zscore(thr(rezMax{thryIdx}{idx}.pos(1):rezMax{thryIdx}{idx}.pos(1)+length(curBar)-1)))

pccScore = zscore(curBar,1)*zscore(thr(pos+rezMax{thryIdx}{idx}.pos(1)-1:pos+rezMax{thryIdx}{idx}.pos(1)+length(curBar)-1-1),1)'/length(curBar);
title(['pcc = ' ,num2str(pccScore)])
xlabel('location x (pixels)')
hold off
% legend({'B(x) 95% CI','B(x)','T(x)'},'location','southoutside')
[~,st,~] =fileparts(bars{idx}.name)
legend({'B(x) +-3sigma',['B(x) ',strrep(st,'_','\_')] ,['T(x)',strrep(theoryStruct(thryIdx).name,'_','\_')]},'FontSize',7,'location','southoutside')


end

