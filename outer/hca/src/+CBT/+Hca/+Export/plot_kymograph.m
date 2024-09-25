function [] = plot_kymograph(idx, kymoStructs,barcodeGen )
    % plot_kymograph
    % plots kymograph for selected idx
    % Plots kymograph for selected barcode
    %     Args:
    %         idx, kymoStructs,barcodeGen 
    % 
    %     Returns:
    % 
    %     Example:
    %

	figure,hold on

    subplot(4,2,[1 2])
    imshow(kymoStructs{idx}.unalignedKymo,[])
    subplot(4,2,[3 4])
    imshow(kymoStructs{idx}.alignedKymo,[]) 
    subplot(4,2,[5 6])
    plot(barcodeGen{idx}.rawBarcode)
    hold on
    plot(find(barcodeGen{idx}.rawBitmask),barcodeGen{idx}.rawBarcode(barcodeGen{idx}.rawBitmask))
    ylabel('intensity')
    xlim([0,length(barcodeGen{idx}.rawBitmask)])
    subplot(4,2,[7 8])
    imshow(repmat(barcodeGen{idx}.rawBarcode,30,1),[])

end

