function [] = plot_features_overlay(hAxisFeaturePaths, imgMat, pathLabelsMat)
    % PLOT_FEATURES_OVERLAY - Plots features on top of a grayscale image
    %
    % Inputs:
    %  hAxisFeaturePaths
    %   handle for axis in which to plot image with paths
    %  imgMat
    %   normalized kymograph, preferably the one in which the features have
    %   been found
    %  pathLabelsMat
    %   matrix of the same size as imgMat where the features are
    %    labeled with a different positive integer and everything else is
    %    zero
    % 
    % Authors:
    %  Saair Quaderi

    
    %Each pixel in inpImg that is crossed by a path, is set to 2


%     %normalize so it's all visible
%     dispArrRGB = imgMat;
%     dispArrRGB = dispArrRGB - min(dispArrRGB(:));
%     dispArrRGB = dispArrRGB./max(dispArrRGB(:));
%     dispArrRGB(pathLabelsMat > 0) = 0;
%     dispArrRGB = repmat(dispArrRGB, [1 1 3]);
% 
%     dispArrRGB = dispArrRGB + label2rgb(pathLabelsMat);
% 
%     %An image is shown, where the paths are colored
%     imshow(dispArrRGB, 'Parent', hAxisFeaturePaths);
%     colormap(hAxisFeaturePaths, lines());
  dispArrRGB = double(imgMat);
    dispArrRGB = dispArrRGB - min(dispArrRGB(:));
    dispArrRGB = dispArrRGB./max(dispArrRGB(:));
    dispArrRGB = round(dispArrRGB*65535);
    dispArrRGB = uint16(dispArrRGB);
    dispArrRGB(pathLabelsMat > 0) = 0;
    dispArrRGB = repmat(dispArrRGB, [1 1 3]);
    %Luis: Added zerocolor as black, so you can actually see something
    dispArrRGB = dispArrRGB + uint16(double(label2rgb(pathLabelsMat,'jet','k'))*(65535/255));

    %An image is shown, where the paths are colored
    imshow(dispArrRGB, 'Parent', hAxisFeaturePaths);
    colormap(hAxisFeaturePaths, lines());
end