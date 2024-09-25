function [kymoStructsUpdated,kymoKeep] = shrink_finder_fun( hcaSets, kymoStructs, saveOutput)
    %% shrink_finder_fun
    % by Luis Mario Leal Garza
    % refactored by AD
    %
    %   Args:
    %       hcaSets
    %       kymoStructs
 
    shrink_threshold = hcaSets.default.shrinkThreshold; % low=strict, high=loose, good=25-40 Minimum improvement in total residual error for each changepoint 
    restrict = hcaSets.default.restrict; % 0-1 value, fraction of time elapsed to find stable region
    gd = hcaSets.default.gd; %slope below this is considered as "shrinking"
    p_feat = hcaSets.default.pFeat; %percentile for the first/last feature, can be 0 if "phantom features" are ignored
    around_feat = hcaSets.default.aroundFeat; %how many more features included in mean
    minStableRegion = hcaSets.default.minStableRegion; % minimum length in timeframes of stable region

    import CBT.Hca.Import.add_kymographs_fun;
    import CBT.Hca.Core.align_kymos;
    import Featuretrack.path_track;

    if nargin < 2
        hcaSets.kymosets.filenames = hcaSets.kymofolder;
        hcaSets.kymosets.kymofilefold = cell(1,length(   hcaSets.kymosets.filenames));
 
       % add kymographs
        [kymoStructs] = add_kymographs_fun(hcaSets);
    end

    if nargin < 3
        saveOutput = 1;
    end

    kymoStructsUpdated = kymoStructs;



%         %  put the kymographs into the structure
%         import CBT.Hca.Core.edit_kymographs_fun;
%         kymoStructs = edit_kymographs_fun(kymoStructs,hcaSets.timeFramesNr);

    % align kymos
    [kymoStructs] = align_kymos(hcaSets,kymoStructs);

    if saveOutput
        %% Set Paths and Make Directories
        hcaSets.output.matDirpath = fileparts(hcaSets.kymofolder{1});
        [~,~] = mkdir(fullfile(hcaSets.output.matDirpath,'goodKymos')); 
        [~,~] = mkdir(fullfile(hcaSets.output.matDirpath,'unfixableKymos')); 
        [~,~] = mkdir(fullfile(hcaSets.output.matDirpath,'fixedKymos'));
        % maybe: can also save folders for "check" things, removed for
        % now
    end

        %% Edge Track / for shift alignet kymo
%         f=1;
%         fn=1;
%         for i=f:fn 
%             feature_paths=kymoStructs{1,i}.featuresIdxs;
%             path_track(uint16(kymoStructs{1,i}.shiftalignedKymo),feature_paths);
%         end
        
%          i=1; 
%           l_edge=kymoStructs{1,i}.leftEdgeIdxs;
%           r_edge=kymoStructs{1,i}.rightEdgeIdxs;
%           feature_paths=[l_edge',r_edge'];
%          %  feature_paths=[l_edge,m_edge,r_edge];
%          % feature_paths=kymoStructs{1,i}.featuresIdxs;
%           plot(barcodeGen{1,i}.rawBarcode)
%           path_track(uint16(kymoStructs{1,i}.alignedKymo),feature_paths);
% %           figure(3)
%           findchangepts(kymosize_overtime{1,i}.used,'Statistic','mean','MinThreshold',shrink_threshold)
%           axis([-inf inf 0.9*mean(kymosize_overtime{1,i}.used) 1.10*mean(kymosize_overtime{1,i}.used)])
%          

     %% REdge and LEdge approach
    import Featuretrack.path_track;
    nKymos=length(kymoStructs);
    kymosize_overtime=cell(1,nKymos);
    fixedkymo_counter=0;
    unfix_counter=0;
    good_counter=0;
    fixed=string();
    unfixables=string();
    goodkymos=string();
    background_thres=2:-0.1:-1; % ?? 

    kymoKeep = nan(nKymos,2); % which rows to keep
% loop through all kymos
    for i = 1:nKymos
        [~,name1,name2] = fileparts(kymoStructs{1,i}.name);

        %  To check for size variation, several segments of the kymo can be
        %  compared (i.e. feature in q1 vs feature in q3 vs median feature)
         nfeat=size(kymoStructs{1,i}.featuresIdxs,2);
        
        % This uses several features at once % First you group a certain amount from the left
        % starting at p_feat and extending to around_feat
        leftFeatures=kymoStructs{1,i}.featuresIdxs(:,max(1,round(nfeat*p_feat)):round(nfeat*around_feat+max(1,round(nfeat*p_feat))));
        [l_edge ] = find_slopes(leftFeatures); % then you calculate their slopes
  
        % Now for right edge
        rightFeatures=kymoStructs{1,i}.featuresIdxs(:,round(nfeat*(1-p_feat)-nfeat*around_feat):round(nfeat*(1-p_feat)));
        [r_edge ] = find_slopes(rightFeatures);
         
        % Now for mid edge
        centerFeatures=kymoStructs{1,i}.featuresIdxs(:,max(1,round(nfeat*(0.5-around_feat/2))):round(nfeat*(0.5+around_feat/2)));
        [m_edge] = find_slopes(centerFeatures);
        
        feature_paths=[l_edge,m_edge,r_edge];
        % This helps you visualize detected edges
        % path_track(uint8(kymoStructs{1,i}.shiftalignedKymo),feature_paths);


        left = m_edge - l_edge; % left side
        right = r_edge - m_edge; % right side
        
        left_slope=polyfit(1:length(left),left,1);
        left_slope=left_slope(1);
        right_slope=polyfit(1:length(right),right,1);
        right_slope=right_slope(1);

        % choose the part with greater negative slope
        if left_slope<right_slope
            used =  left;
            going_down = left_slope;
        else
            used = right;
            going_down = right_slope;
        end

        % This thing here creates a stability plot. Divides the sizeplot in
        % segments according to stability. Gets more sensitive with a lower
        % "shrink_threshold". I also refer to this as "change analysis" later
        % down in the code
        [ipt,residual]=findchangepts(used,'Statistic','mean','MinThreshold',shrink_threshold);

        %% Kymo Rescue % maybe kill them if its going smaller from both sides
        if size(ipt,1)>0 && going_down<gd  % After change analysis, if there more than one "segment" then it is already not that good
            ipt=[1,transpose(ipt),size(used,1)];
            diffs=diff(ipt);
            result=diffs(diffs>=minStableRegion); % This looks for the latest stable (>10 frames) section
            if (size(result,2) > 0) %if there was any stable (>10 section) it checks if its within the restricted limits (towards the end)
                stablemin=find(diffs==result(size(result,2))); % This is the index were the stable region starts
                stablemin=stablemin(size(stablemin,2));
                stablemax=stablemin+1; %so if there is a most stable region towards the end, lets take it and do a slope check
                slope_check=polyfit((ipt(stablemin)+1):ipt(stablemax),used((ipt(stablemin)+1):ipt(stablemax)),1);
                slope_check=slope_check(1)>gd; %if the last "stable" segment, is still shrinking, dump it
                end_check=ipt(stablemax)>=(size(used,1)*restrict);% check if stable min is on the last 10% frames? if not, unfixable
            else % If there was no stable region, just dump it
                end_check=0;
                slope_check=0;
            end
            if (end_check && slope_check) % if stable region within the restricted end and slope check passed, extract them
                fixedkymo_counter=fixedkymo_counter+1;
                fixed(fixedkymo_counter)=kymoStructs{1,i}.name;
                kymoKeep(i,:) = [(ipt(stablemin)+1) ipt(stablemax)];
                unalignedKymo =double(kymoStructs{1,i}.unalignedKymo((ipt(stablemin)+1):ipt(stablemax),:));
                try
                unalignedBitmask =kymoStructs{1,i}.unalignedBitmask((ipt(stablemin)+1):ipt(stablemax),:);
                catch
                end
                enhanced  =  imadjust(unalignedKymo/max(unalignedKymo(:)),[0.1 1]);
                if saveOutput % todo: move this outside
                    folder='fixedKymos';
                    save_path=fullfile(hcaSets.output.matDirpath, folder,strcat(name1,name2));
                    imwrite(uint16(round(double(enhanced)./max(enhanced(:))*2^16)),save_path,'tiff') % Enhanced
                    imwrite(uint16(unalignedKymo),save_path,'tiff', 'writemode', 'append') % raw
                    imwrite(uint16(unalignedBitmask),save_path,'tiff', 'writemode', 'append') % raw  
                end
             else
                unfix_counter=unfix_counter+1;
                unfixables(unfix_counter)=kymoStructs{1,i}.name;
                if saveOutput
                    folder='unfixableKymos/';
                    save_path=fullfile(hcaSets.output.matDirpath,folder,strcat(name1,name2));
                    copyfile(kymoStructs{1,i}.name,save_path);   
                end
            end
    else  %Allocate good kymos somewhere else, also with QC
        good_counter=good_counter+1;
        goodkymos(good_counter)=kymoStructs{1,i}.name;
        kymoKeep(i,:)  = [1 size(kymoStructs{1,i}.unalignedKymo,1)];
        if saveOutput
            folder='goodKymos';
            save_path=fullfile(hcaSets.output.matDirpath,folder,strcat(name1,name2));
            copyfile(kymoStructs{1,i}.name,save_path);
        end
    end

    if ~isnan(kymoKeep(i,1))
        kymoStructsUpdated{i}.unalignedKymo = kymoStructsUpdated{i}.unalignedKymo(kymoKeep(i,1):kymoKeep(i,2),:);
        if isfield(kymoStructsUpdated{i},'unalignedBitmask')
            kymoStructsUpdated{i}.unalignedBitmask = kymoStructsUpdated{i}.unalignedBitmask(kymoKeep(i,1):kymoKeep(i,2),:);
        end
        if isfield(kymoStructsUpdated{i},'leftEdgeIdxs')
            kymoStructsUpdated{i}.leftEdgeIdxs = kymoStructsUpdated{i}.leftEdgeIdxs(kymoKeep(i,1):kymoKeep(i,2));
            kymoStructsUpdated{i}.rightEdgeIdxs = kymoStructsUpdated{i}.rightEdgeIdxs(kymoKeep(i,1):kymoKeep(i,2));
        end
    else
    kymoStructsUpdated{i}.unalignedKymo = [];
    kymoStructsUpdated{i}.unalignedBitmask = [];
    kymoStructsUpdated{i}.leftEdgeIdxs =[];
    kymoStructsUpdated{i}.rightEdgeIdxs = [];
   
    end
    
    kymoStructsUpdated{i}.feature_paths = feature_paths;
    kymoStructsUpdated{i}.featuresIdxs = kymoStructs{1,i}.featuresIdxs;
    kymoStructsUpdated{i}.shiftalignedKymo = kymoStructs{1,i}.shiftalignedKymo;
    kymoStructsUpdated{i}.left_slope = left_slope;
    kymoStructsUpdated{i}.right_slope = right_slope;


%% Save for QC
% Put after analysis so folder gets sorted by criteria
%  save_path=[sets.output.matDirpath,folder,'Edge_check/',"Featurewise_edges_",kymoStructs{1,i}.name];
%  path_track(uint16(kymoStructs{1,i}.unalignedKymo),feature_paths);
%  imwrite(getframe(gcf).cdata,strjoin(save_path,""))
%  close
    end

    
display(['Fixed ', num2str(fixedkymo_counter), ' kymos'])
display(['Total ', num2str(good_counter), ' good kymos'])
display(['Total ', num2str(unfix_counter), ' unfixable kymos'])



    


    function [l_edge ] = find_slopes(group_of_kymos)
        b=[];
        for m=1:size(group_of_kymos,2)
            a=polyfit(1:length(group_of_kymos(:,m)),group_of_kymos(:,m),1);
            b(m)=a(1);
        end
        % Then you cluster them by slope (to leave out weird features)
        epsilon = mean(diff(sort(b))); %how similar slopes should be to belong together
        idx =dbscan(b',epsilon,1);
        group_of_kymos = group_of_kymos(:,(idx==mode(idx)));
        % Now we got the similar features (leaving out weird ones) we get average
        l_edge = round(mean(group_of_kymos,2));
    end
    
end
% end

