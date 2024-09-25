function [] = path_track(Kymo, pathsColIdxs)
% Takes any paths and any kymo and sticks them together    
     szImg = size(Kymo);
     pathLabelsMat = zeros(szImg);
     numPaths = size(pathsColIdxs, 2);
     pathRowIdxs = 1:szImg(1);
     for pathNum = 1:numPaths
         pathColIdxs = pathsColIdxs(:, pathNum);
         pathLabelsMat(sub2ind(szImg, pathRowIdxs, pathColIdxs.')) = pathNum;
     end
     
     hFig = figure();
     hPanel = uipanel('Parent', hFig);
     hAxis = axes('Parent', hPanel);
     import Featuretrack.plot_features_overlay;
     plot_features_overlay(hAxis, Kymo, pathLabelsMat);
end
