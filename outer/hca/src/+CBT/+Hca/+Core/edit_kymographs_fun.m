function [ kymoStructs,passingKymos ] = edit_kymographs_fun(kymoStructs, timeFramesNr,timeFrameIdx )
    % edit_kymographs_fun
    % Function for putting the kymographs into a session file
    %
    %     Args:
    %         kymoStructs: Kymograph structure
    %         timeFramesNr: number of timeframes
    % 
    %     Returns:
    %         kymoStructs: kymoStructs
    %

    if nargin < 3
        timeFrameIdx = 1; % start of time-frame idx
    end
    % If the setting for number of timeframes is non-zero, we have to
    % remove some rows from the kymographs, and remove the kymographs with
    % insuficient amount of rows altogether
    try
        timeframeTotal = cellfun(@(x) size(x.unalignedKymo, 1), kymoStructs)-timeFrameIdx+1;
    catch
        timeframeTotal = 0; % in case no kymos
    end

    if timeFramesNr ~= 0
        % number of timegrames
        timeframeTotal = cellfun(@(x) size(x.unalignedKymo, 1), kymoStructs);
        % remove kymo's that have not enough timeframes
        kymoStructs(timeframeTotal < timeFramesNr) = [];
        % dispplay which kymographs were removed
        disp(strcat(['Kymographs nr ' num2str(find(timeframeTotal < timeFramesNr)) ' were removed because they do not have enough time-frames'])); 
        % save only the preselected number of rows
        for i=1:length(kymoStructs)
%             try
                kymoStructs{i}.unalignedKymo = kymoStructs{i}.unalignedKymo(timeFrameIdx:timeFramesNr+timeFrameIdx-1,:);
%             catch
%                 kymoStructs{i} = [];
%             end
        end        
    end    
    passingKymos = timeframeTotal >= timeFramesNr;
    
end

