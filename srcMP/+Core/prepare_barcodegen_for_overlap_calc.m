function [SCAMP_LINE, barStruct, namesBar, stridx, names2, baridx2] = prepare_barcodegen_for_overlap_calc(barcodeGen, numWorkers, sF, foldSynth)

%   Args:
%       barcodeGen, numWorkers, sF, foldSynth

%   Returns:
%       SCAMP_LINE, barStruct, namesBar, stridx, names2, baridx2
%   
%
%       stridx - saves stretch parameter for barcodes starting at
%       particular locations
%% ASSEMBLY MP
delete(gcp('nocreate'))
parpool('local',numWorkers)

% numWorkers = 28;

% addpath(genpath('C:\Users\Lenovo\git\bargroupingprototype'))
% cd C:\Users\Lenovo\git\bargroupingprototype

if ispc
    SCAMP_LINE = '.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else
    SCAMP_LINE = '~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
end

% MIN_OVERLAP_PIXELS = 150;
% sF = 0.85:0.01:1.15;

%%

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
% have to save separately..
% foldSynth= 'barcoli';
% have to save separately..
[namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);



[names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);


end

