function [] = plot_kymo_masks(kymoStructs,ii)
%   plot_kymo_masks - plots kymographs with masks
%
%   Args:
%       kymoStructs - kymostructures
%       ii - which dataset to plot

    figure,tiledlayout(ceil(sqrt(length(kymoStructs{ii}))),ceil(length(kymoStructs{ii})/sqrt(length(kymoStructs{ii}))),'TileSpacing','none','Padding','none')
    for i=1:length(kymoStructs{ii})
    %         nexttile;        imagesc(imresize(kymoStructs{ii}{i}.unalignedKymo,[200 500]));    title(num2str(i));
            
            nexttile;       
            if ~isempty(kymoStructs{ii}{i}.unalignedBitmask)
                imshowpair(imresize(kymoStructs{ii}{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{ii}{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  );    title(num2str(i));
            else
                axis off
            end
    end
end

