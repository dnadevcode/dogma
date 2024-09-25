function [] = shrink_finder_featurwise()
%% Shrink_finder_edge_remover
% by Luis Mario Leal Garza

% Keep HCA in path
% In short, this script removes kymos that touch any of the edges
% also removes manually cropped kymos.
% Then it detects edges by contrast (not using HCA output)
% Finally checks if Kymos shrink and try to fix

% The input for this script is:
% 1. Raw Kymos
% 2. Output Folder
% These three parameters are the only "movables"
shrink_threshold=40; % low=strict, high=loose, good=25-40 Minimum improvement in total residual error for each changepoint 
consecutive_bgthreshold=100; %consecutivepx over background to consider contrast as ok
consecutive_edge_threshold=[4 6]; %1st number is consecutive px in a window of size (2nd number) to consider the start/end of a molecule
restrict=0.9; % 0-1 value, fraction of time elapsed to find stable region
% 0=stable region can be anywhere
% 0.6=stable region must be until the last 40% of time
% 0.9=stable region must be until the last 10% of time
gd=-0.2; %slope below this is considered as "shrinking"
p_feat=0; %percentile for the first/last feature
% Can be 0 if "phantom features" are ignored
around_feat=0.10; %how many more features included in mean

% The output of this script is three directories and 2 subdirectories:
% --------Directories------------------------------------------------
% GOOD KYMOS:Kymos that are considered stable
% FIXED KYMOS: Shrinking kymos which had the latest >10 frames stable part
% extracted.
% UNFIXABLE KYMOS: Kymos that were continuously shrinking and have no 
% stable part
% --------SubDirectories------------------------------------------------
% MANUAL_CHECK: Inside all directories, shows the stable section
% selection plot.
% EDGE_CHECK: Shows the features used as left, right and mid edge.
% This script uses HCA importing functions and kymograph processing
% to measure each kymo and get the "size vs time" curve of each.
% Then findchangepts() is used to find stable sections to classify
% each kymograph.


%% Import kymos
% It is just a copy-paste from HCA_Gui
        import CBT.Hca.Import.import_hca_settings;
        [sets] = import_hca_settings('hca_settings.txt');
         
        import CBT.Hca.Settings.get_user_settings;
        sets = get_user_settings(sets);
         
        % timestamp for the results
        sets.timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
 
        % add kymographs
        import CBT.Hca.Import.add_kymographs_fun;
        [kymoStructs] = add_kymographs_fun(sets);

        %  put the kymographs into the structure
        import CBT.Hca.Core.edit_kymographs_fun;
        kymoStructs = edit_kymographs_fun(kymoStructs,sets.timeFramesNr);

        % align kymos
        import CBT.Hca.Core.align_kymos;
        [kymoStructs] = align_kymos(sets,kymoStructs);
  
 %Include this after everything to be compatible with pipeline
        % generate barcodes
         import CBT.Hca.Core.gen_barcodes;
         barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);
 %% Set Paths and Make Directories
mkdir([sets.output.matDirpath,'goodKymos']); 
mkdir([sets.output.matDirpath,'unfixableKymos']); 
mkdir([sets.output.matDirpath,'fixedKymos']);
mkdir([sets.output.matDirpath,'fixedKymos/Manual_check']);
mkdir([sets.output.matDirpath,'unfixableKymos/Manual_check']);
mkdir([sets.output.matDirpath,'unfixableKymos/Edge_check']);
mkdir([sets.output.matDirpath,'fixedKymos/Edge_check']);
mkdir([sets.output.matDirpath,'goodKymos/Edge_check']);
mkdir([sets.output.matDirpath,'goodKymos/Manual_check']);
%% Edge Track
% f=1;
% fn=50;
% for i=f:n 
% % This helps you visualize edges in kymo (by traditional edge detect)
%  l_edge=kymoStructs{1,i}.leftEdgeIdxs;
%  r_edge=kymoStructs{1,i}.rightEdgeIdxs;
% % feature_paths=[l_edge;r_edge].';
% %  feature_paths=[l_edge,m_edge,r_edge];
%  feature_paths=kymoStructs{1,i}.featuresIdxs;
%  path_track(uint16(kymoStructs{1,i}.unalignedKymo),feature_paths);
% % imagesc(kymoStructs{1,i}.unalignedKymo)
% % imagesc(kymoStructs{1,i}.alignedKymo)
% % path_track(uint8(kymoStructs{1,i}.shiftalignedKymo),feature_paths);
% end
%% Test
%  i=2; 
%   l_edge=kymoStructs{1,i}.leftEdgeIdxs;
%   r_edge=kymoStructs{1,i}.rightEdgeIdxs;
%   feature_paths=[l_edge',r_edge'];
%  %  feature_paths=[l_edge,m_edge,r_edge];
%  % feature_paths=kymoStructs{1,i}.featuresIdxs;
%   plot(barcodeGen{1,i}.rawBarcode)
%   path_track(uint16(kymoStructs{1,i}.unalignedKymo),feature_paths);
%   figure(3)
%   findchangepts(kymosize_overtime{1,i}.used,'Statistic','mean','MinThreshold',shrink_threshold)
%   axis([-inf inf 0.9*mean(kymosize_overtime{1,i}.used) 1.10*mean(kymosize_overtime{1,i}.used)])
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
background_thres=2:-0.1:-1;

for i = 1:nKymos
%% Ignore phantom features %Outdated
% % This is introduced because many of the first and last features are just
% % noise. Lets subset for the features "inside" the kymo
% % Edges are not easy to find, but how about "bright spots"
% bright_idx=kymoStructs{1,i}.shiftalignedKymo(1,:)>nanmean(kymoStructs{1,i}.shiftalignedKymo(1,:));
% bright_idx=find(bright_idx==1);
% bright_idx=min(bright_idx):max(bright_idx);
% kymoStructs{1,i}.featuresIdxs
% bright_paths=ismember(kymoStructs{1,i}.featuresIdxs(1,:),bright_idx);
% kymoStructs{1,i}.featuresIdxs=kymoStructs{1,i}.featuresIdxs(:,bright_paths);

%% Feature wise edges
%  Taking first and last identified feature should be better than "edges"
%  Features identified are sometimes just noise, so first and last several
%  features should be ignored.
%  It it then a parameter to choose which percentile of the features should
%  be ignored on each side.
%  To check for size variation, several segments of the kymo can be
%  compared (i.e. feature in q1 vs feature in q3 vs median feature)
 nfeat=size(kymoStructs{1,i}.featuresIdxs,2);

% This uses several features at once
% l_edge=kymoStructs{1,i}.featuresIdxs(:,round(nfeat*p_feat));
% r_edge=kymoStructs{1,i}.featuresIdxs(:,round(nfeat*(1-p_feat)));

% This uses several features at once
% First you group a certain amount from the left
% starting at p_feat and extending to around_feat
if p_feat==0
    group_of_kymos=kymoStructs{1,i}.featuresIdxs(:,1:round(nfeat*around_feat+1));
else
    group_of_kymos=kymoStructs{1,i}.featuresIdxs(:,round(nfeat*p_feat):round(nfeat*around_feat+round(nfeat*p_feat)));
end
% then you calculate their slopes
b=[];
for m=1:size(group_of_kymos,2)
    a=polyfit(1:length(group_of_kymos(:,m)),group_of_kymos(:,m),1);
    b(m)=a(1);
end
% Then you cluster them by slope (to leave out weird features)
epsilon=mean(diff(sort(b))); %how similar slopes should be to belong together
idx=dbscan(b',epsilon,1);
group_of_kymos=group_of_kymos(:,(idx==mode(idx)));
% Now we got the similar features (leaving out weird ones) we get average
l_edge=round(mean(group_of_kymos,2));

% Now for right edge
group_of_kymos=kymoStructs{1,i}.featuresIdxs(:,round(round(nfeat*(1-p_feat))-nfeat*around_feat):round(nfeat*(1-p_feat)));
% then you calculate their slopes
b=[];
for m=1:size(group_of_kymos,2)
    a=polyfit(1:length(group_of_kymos(:,m)),group_of_kymos(:,m),1);
    b(m)=a(1);
end
% Then you cluster them by slope (to leave out weird features)
epsilon=mean(diff(sort(b))); %how similar slopes should be to belong together
idx=dbscan(b',epsilon,1);
group_of_kymos=group_of_kymos(:,(idx==mode(idx)));
% Now we got the similar features (leaving out weird ones) we get average
r_edge=round(mean(group_of_kymos,2));
 
% Now for mid edge
group_of_kymos=kymoStructs{1,i}.featuresIdxs(:,round(nfeat*(0.5-around_feat/2)):round(nfeat*(0.5+around_feat/2)));
% then you calculate their slopes
b=[];
for m=1:size(group_of_kymos,2)
    a=polyfit(1:length(group_of_kymos(:,m)),group_of_kymos(:,m),1);
    b(m)=a(1);
end
% Then you cluster them by slope (to leave out weird features)
epsilon=mean(diff(sort(b))); %how similar slopes should be to belong together
idx=dbscan(b',epsilon,1);
group_of_kymos=group_of_kymos(:,(idx==mode(idx)));
% Now we got the similar features (leaving out weird ones) we get average
m_edge=round(mean(group_of_kymos,2));

% This helps you visualize edges in kymo (by traditional edge detect)
% l_edge=kymoStructs{1,i}.leftEdgeIdxs;
% r_edge=kymoStructs{1,i}.rightEdgeIdxs;
 feature_paths=[l_edge,m_edge,r_edge];

% This helps you visualize all features paths
% feature_paths=kymoStructs{1,i}.featuresIdxs;
% path_track(uint8(kymoStructs{1,i}.alignedKymo),feature_paths);
% path_track(uint8(kymoStructs{1,i}.shiftalignedKymo),feature_paths);

% I have to put all of these above into kymoStructure so it can be
% combined with the other code.
kymoStructs{1,i}.rightEdgeIdxs=r_edge;
kymoStructs{1,i}.leftEdgeIdxs=l_edge;
kymoStructs{1,i}.midEdgeIdxs=m_edge;
%% Get Size over time
% Now we have two sections L and R
kymosize_overtime{1,i}.left=kymoStructs{1,i}.midEdgeIdxs-kymoStructs{1,i}.leftEdgeIdxs;
kymosize_overtime{1,i}.right=kymoStructs{1,i}.rightEdgeIdxs-kymoStructs{1,i}.midEdgeIdxs;
    %% Is it shrinking?
    % We check here both segments and overall, and if  choose the one with greater negative
    % slope, and keep the most negative overall
    left_slope=polyfit(1:length(kymosize_overtime{1,i}.left),kymosize_overtime{1,i}.left,1);
    left_slope=left_slope(1);
    right_slope=polyfit(1:length(kymosize_overtime{1,i}.right),kymosize_overtime{1,i}.right,1);
    right_slope=right_slope(1);
    if left_slope<right_slope
        kymosize_overtime{1,i}.used=kymosize_overtime{1,i}.left;
        going_down=left_slope;
    else
        kymosize_overtime{1,i}.used=kymosize_overtime{1,i}.right;
        going_down=right_slope;
    end
    % This thing here creates a stability plot. Divides the sizeplot in
    % segments according to stability. Gets more sensitive with a lower
    % "shrink_threshold". I also refer to this as "change analysis" later
    % down in the code
    [ipt,residual]=findchangepts(kymosize_overtime{1,i}.used,'Statistic','mean','MinThreshold',shrink_threshold);
    %% Kymo Rescue % maybe kill them if its going smaller from both sides
    if size(ipt,1)>0 && going_down<gd  % After change analysis, if there more than one "segment" then it is already not that good
        ipt=[1,transpose(ipt),size(kymosize_overtime{1,i}.used,1)];
        diffs=diff(ipt);
       %result=max(diffs(diffs>=0)); %this looks for the biggest stable (>10 frames) section
        result=diffs(diffs>=10); % This looks for the latest stable (>10 frames) section
         if (size(result,2) > 0) %if there was any stable (>10 section) it checks if its within the restricted limits (towards the end)
          stablemin=find(diffs==result(size(result,2))); % This is the index were the stable region starts
          stablemin=stablemin(size(stablemin,2));
          stablemax=stablemin+1; %so if there is a most stable region towards the end, lets take it and do a slope check
          slope_check=polyfit((ipt(stablemin)+1):ipt(stablemax),kymosize_overtime{1,i}.used((ipt(stablemin)+1):ipt(stablemax)),1);
          slope_check=slope_check(1)>gd; %if the last "stable" segment, is still shrinking, dump it
          end_check=ipt(stablemax)>=(size(kymosize_overtime{1,i}.used,1)*restrict);% check if stable min is on the last 10% frames? if not, unfixable
         else % If there was no stable region, just dump it
          end_check=0;
          slope_check=0;
         end
            if (end_check && slope_check) % if stable region within the restricted end and slope check passed, extract them
                
                fixedkymo_counter=fixedkymo_counter+1;
                fixed(fixedkymo_counter)=kymoStructs{1,i}.name;
                fixedkymoStructs{1,fixedkymo_counter}.name =kymoStructs{1,i}.name;
                fixedkymoStructs{1,fixedkymo_counter}.unalignedKymo =kymoStructs{1,i}.unalignedKymo((ipt(stablemin)+1):ipt(stablemax),:);
                fixedkymoStructs{1,fixedkymo_counter}.alignedKymo =kymoStructs{1,i}.alignedKymo((ipt(stablemin)+1):ipt(stablemax),:);
                fixedkymoStructs{1,fixedkymo_counter}.leftEdgeIdxs =kymoStructs{1,i}.leftEdgeIdxs((ipt(stablemin)+1):ipt(stablemax),:);
                fixedkymoStructs{1,fixedkymo_counter}.rightEdgeIdxs =kymoStructs{1,i}.rightEdgeIdxs((ipt(stablemin)+1):ipt(stablemax),:);
                try
                    fixedkymoStructs{1,fixedkymo_counter}.unalignedBitmask =kymoStructs{1,i}.unalignedBitmask((ipt(stablemin)+1):ipt(stablemax),:);
                end    
                            % enhanced
                sampIm = mat2gray(fixedkymoStructs{1,fixedkymo_counter}.unalignedKymo);
                minInt = min(sampIm(:));
                medInt = median(sampIm(:));
%               maxInt = max(sampIm(:));
                try
                    J = imadjust(sampIm,[minInt 4*medInt]);
                catch
                    J =  imadjust(sampIm,[0.1 0.9]);
                end
                fixedkymoStructs{1,fixedkymo_counter}.enhanced = J;
                
                %help see were stable section starts and ends
                %index_ref(fixedkymo_counter,1)=ipt(stablemin);
                %index_ref(fixedkymo_counter,2)=ipt(stablemax);
                %index_ref(fixedkymo_counter,3)=size(kymosize_overtime{1,i},1);
                folder='fixedKymos/';
                save_path=[sets.output.matDirpath,folder,fixedkymoStructs{1,fixedkymo_counter}.name];
                imwrite(uint16(fixedkymoStructs{1,fixedkymo_counter}.enhanced*(power(2,16)-1)),save_path,'tiff') % Enhanced
                imwrite(uint16(fixedkymoStructs{1,fixedkymo_counter}.unalignedKymo),save_path,'tiff', 'writemode', 'append') % raw
                try
                 imwrite(uint16(fixedkymoStructs{1,fixedkymo_counter}.unalignedBitmask),save_path,'tiff', 'writemode', 'append') % bitmask
                end
                
                findchangepts(kymosize_overtime{1,i}.used,'Statistic','mean','MinThreshold',shrink_threshold)
                hold on
                rectangle('Position',[(ipt(stablemin)-0.5) min(kymosize_overtime{1,i}.used((ipt(stablemin)+1):ipt(stablemax))) (ipt(stablemax)-ipt(stablemin)) range(kymosize_overtime{1,i}.used(ipt(stablemin):ipt(stablemax)))],'Curvature',0.2,'EdgeColor','blue')
                title(fixedkymoStructs{1,fixedkymo_counter}.name, 'Interpreter', 'none')
                text(ipt(stablemin),mean([max(kymosize_overtime{1,i}.used), max(kymosize_overtime{1,i}.used(ipt(stablemin):ipt(stablemax)))]),"Selected Region")
                save_path=[sets.output.matDirpath,folder,'Manual_check/',"Analysis_of_",fixedkymoStructs{1,fixedkymo_counter}.name];
                imwrite(getframe(gcf).cdata,strjoin(save_path,""))
                close
      % Unfixables dump        
             else
                unfix_counter=unfix_counter+1;
                unfixables(unfix_counter)=kymoStructs{1,i}.name;
                unfixableskymoStructs{1,unfix_counter}.name =kymoStructs{1,i}.name;
                unfixableskymoStructs{1,unfix_counter}.unalignedKymo =kymoStructs{1,i}.unalignedKymo;
                unfixableskymoStructs{1,unfix_counter}.alignedKymo =kymoStructs{1,i}.alignedKymo;
                unfixableskymoStructs{1,unfix_counter}.leftEdgeIdxs =kymoStructs{1,i}.leftEdgeIdxs;
                unfixableskymoStructs{1,unfix_counter}.rightEdgeIdxs =kymoStructs{1,i}.rightEdgeIdxs;
                try
                    unfixableskymoStructs{1,unfix_counter}.unalignedBitmask =kymoStructs{1,i}.unalignedBitmask;
                end
                % enhanced
                sampIm = mat2gray(unfixableskymoStructs{1,unfix_counter}.unalignedKymo);
                minInt = min(sampIm(:));
                medInt = median(sampIm(:));
%               maxInt = max(sampIm(:));
                try
                    J = imadjust(sampIm,[minInt 4*medInt]);
                catch
                    J =  imadjust(sampIm,[0.1 0.9]);
                end
                unfixableskymoStructs{1,unfix_counter}.enhanced = J;
                
                folder='unfixableKymos/';
                save_path=[sets.output.matDirpath,folder,unfixableskymoStructs{1,unfix_counter}.name];
                imwrite(uint16(unfixableskymoStructs{1,unfix_counter}.enhanced*(power(2,16)-1)),save_path,'tiff') % Enhanced
                imwrite(uint16(unfixableskymoStructs{1,unfix_counter}.unalignedKymo),save_path,'tiff', 'writemode', 'append') % raw
                try
                    imwrite(uint16(unfixableskymoStructs{1,unfix_counter}.unalignedBitmask),save_path,'tiff', 'writemode', 'append') % bitmask
                end
                findchangepts(kymosize_overtime{1,i}.used,'Statistic','mean','MinThreshold',shrink_threshold)
                title(unfixableskymoStructs{1,unfix_counter}.name, 'Interpreter', 'none')
                save_path=[sets.output.matDirpath,folder,'Manual_check/',"Analysis_of_",unfixableskymoStructs{1,unfix_counter}.name];
                imwrite(getframe(gcf).cdata,strjoin(save_path,""))
                close        
            end
        % To remove bad kymo from struct
        % kymoStructs{1,i}=[];
    else  %Allocate good kymos somewhere else, also with QC
    good_counter=good_counter+1;
    goodkymos(good_counter)=kymoStructs{1,i}.name;
    goodkymoStructs{1,good_counter}.name =kymoStructs{1,i}.name;
    goodkymoStructs{1,good_counter}.unalignedKymo =kymoStructs{1,i}.unalignedKymo;
    goodkymoStructs{1,good_counter}.alignedKymo =kymoStructs{1,i}.alignedKymo;
    goodkymoStructs{1,good_counter}.leftEdgeIdxs =kymoStructs{1,i}.leftEdgeIdxs;
    goodkymoStructs{1,good_counter}.rightEdgeIdxs =kymoStructs{1,i}.rightEdgeIdxs;
    try    
        goodkymoStructs{1,good_counter}.unalignedBitmask =kymoStructs{1,i}.unalignedBitmask;
    end
    % enhanced
    sampIm = mat2gray(goodkymoStructs{1,good_counter}.unalignedKymo);
    minInt = min(sampIm(:));
    medInt = median(sampIm(:));
%   maxInt = max(sampIm(:));
       try
          J = imadjust(sampIm,[minInt 4*medInt]);
       catch
          J =  imadjust(sampIm,[0.1 0.9]);
       end
    goodkymoStructs{1,good_counter}.enhanced = J;
    folder='goodKymos/';
    save_path=[sets.output.matDirpath,folder,goodkymoStructs{1,good_counter}.name];
    imwrite(uint16(goodkymoStructs{1,good_counter}.enhanced*(power(2,16)-1)),save_path,'tiff') % Enhanced
    imwrite(uint16(goodkymoStructs{1,good_counter}.unalignedKymo),save_path,'tiff', 'writemode', 'append') % raw
    try
        imwrite(uint16(goodkymoStructs{1,good_counter}.unalignedBitmask),save_path,'tiff', 'writemode', 'append') % bitmask
    end            
    findchangepts(kymosize_overtime{1,i}.used,'Statistic','mean','MinThreshold',shrink_threshold)
    title(goodkymoStructs{1,good_counter}.name, 'Interpreter', 'none')
    save_path=[sets.output.matDirpath,folder,'Manual_check/',"Analysis_of_",goodkymoStructs{1,good_counter}.name];
    imwrite(getframe(gcf).cdata,strjoin(save_path,""))
    close
    end
%% Save for QC
% Put after analysis so folder gets sorted by criteria
 save_path=[sets.output.matDirpath,folder,'Edge_check/',"Featurewise_edges_",kymoStructs{1,i}.name];
 path_track(uint16(kymoStructs{1,i}.unalignedKymo),feature_paths);
 imwrite(getframe(gcf).cdata,strjoin(save_path,""))
 close
end
    

%% Show final results
disp("Out of the "+i+" Kymos, "+good_counter+" are good, "+fixedkymo_counter+" were bad but fixed, and "+unfix_counter+" are unfixable.")
disp('')
if fixed=="" 
   disp('No Kymos were fixed') 
else
    disp('The following Kymos were shrinking but were fixed')
    disp(transpose(fixed))
end
%disp(index_ref)
disp('')
if unfixables==""
    disp('No Kymos were unfixable')
else
    disp('The following Kymos were shrinking and were impossible to fix')
    disp(transpose(unfixables))
end
end