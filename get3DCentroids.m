function centroidsOUT=get3DCentroids(centroidsIN,areasIN,pixflagsIN,zStepSize,twoDstereoflag,exclusionLimits)
% Takes 2D centroids (output from segmentObjectsByClustering) and groups by
% proximity in XYZ to identify clusters of centroids that are probably a
% single object. Assigns 3D object numbers and finds XYZ centroid position.
% Author: J. Kelly
%
% Input: 
% centroidsIN: value assigned to "centroidsOut" by segmentObjectsByClustering 
%   (XY centroids of 2D objects in each z plane)
% areasIN: value assigned to "areasOut" by segmentObjectsByClustering
%   (number of pixels of 2D objects; indexing matches centroidsIN)
% pixflagsIN: value assigned to "pixflagsOut" by segmentObjectsByClustering
%   (1 if any pixels in a 2D object are on rightmost or lowermost edge of
%   image, else 0; indexing matches centroidsIN)
% zStepSize: z step size in image acquisition
% 2Dstereoflag: 1 if user wants to exclude objects that intersect rightmost
%   and/or lowermost edges of image (applied after 3D grouping)
% exclusionLimits: 2-vector with z planes to include (first) and exclude
%   (second); if not assigned, no exclusion in z is applied 
%
% Output:
% centroidsOUT is a structure with four fields: 'twoD','threeD','nonemptyZ','stereo'
%   centroidsOUT.twoD: Nx11 matrix, with each row representing a 2D
%       centroid from centroidsIN, and columns as follows:
%       1: X position
%       2: Y position
%       3: Z plane
%       4: 1 if from centroidsIN.replaced, 2 if from centroidsIN.segmented
%       5: object label number; initialized such that each 2D centroid has
%           its own object label, but replaced with group label once 
%           centroid is added to a group
%       6: 0 if from centroidsIN.replaced, 1 if from centroidsIN.segmented
%       7: # times grouped; initialized to 0, should end as 0 or 1
%       8: XY distance to centroid grouped with this one in z plane above it
%       9: centroid index within its z plane, separate for replaced/seg.
%       10: area of 2D object (from areasIN)
%       11: intersects right/lower edge, y/n? (from pixflagsIN)
%   centroidsOUT.threeD: Nx7 matrix
%       one row per object, columns as follows: 
%       1: object label number (matches column 5 from centroidsOUT.twoD)
%       2: X position of 3D centroid 
%       3: Y position of 3D centroid
%       4: Z position of 3D centroid
%       5: object type 0-4; 0=truncated with widest plane at boundary of
%           stack; 1=bottom contained within stack; 2=top contained within
%           stack; 3=widest point is within stack and within object;
%           4=truncated by stack and has more than one area "zero crossing"
%       6: pixflags (max)
%       7: maximum area of 2D objects in group
%   centroidsOUT.nonemptyZ: z planes with 2D centroids assigned
%   centroidsOUT.stereo: same as centroidsOUT.threeD, but exclusion
%       boundaries have been applied
fullNargin=6;

% % centroidsIN.segmented{z}=new centroids after 2D segmenting 
% % centroidsIN.replaced{z}=centroids output of regionprops that were replaced with segmented centroids
nonemptyZ=find(~cellfun(@isempty,centroidsIN.segmented));
areaMeans=cellfun(@median,areasIN.original(nonemptyZ));
meanParticleDiam=median(sqrt(areaMeans/pi))*2;  % mean across all z, scalar
pix_um_Conversion=1024/246.03; 

zStepSize=zStepSize*pix_um_Conversion;%/(5.25/NA);
centroidDistCutoff=0.5*pix_um_Conversion;% centroidDistCutoff=0.2*pix_um_Conversion; %centroidDistCutoff=meanParticleDiam*5/12;
centroidDistCutoff3d=meanParticleDiam*10/12;%topBottomMatchCutoff=meanParticleDiam*testIN(2)/12;
centroidDistCutoff3d_2=meanParticleDiam*8/12;
centroidDistCutoff3dXY=meanParticleDiam*4/12;centroidDistCutoff3dZ= meanParticleDiam*12/12;
centroidDistCutoff3dXY_2=meanParticleDiam*8/12;centroidDistCutoff3dZ_2= 0.5;


if nargin<fullNargin   % no inclusion/exclusion planes entered, do not apply any
    threeDstereoflag=0;
else
    threeDstereoflag=1;
end

% convert z information to matrix column from cell struct
segmentedXYZ=nan(sum(cellfun(@length,centroidsIN.segmented)),11);
segmentedXYZ(1:size(centroidsIN.segmented{nonemptyZ(1)},1),5)=(1:size(centroidsIN.segmented{nonemptyZ(1)},1))';  % object labels
segmentedXYZ(:,6)=1;    % 0/1 for non-/segmented
segmentedXYZ(:,7)=0;    % number of times grouped
segmentedXYZ(1:size(centroidsIN.segmented{nonemptyZ(1)},1),8)=0;  % distances
% replacedXYZ(1:size(centroidsIN.replaced{nonemptyZ(1)},1),9)=(1:size(centroidsIN.replaced{nonemptyZ(1)},1))';
segmentedXYZ(1:size(centroidsIN.segmented{nonemptyZ(1)},1),9)=(1:size(centroidsIN.segmented{nonemptyZ(1)},1))'; % redundant but don't leave as nan

for z=nonemptyZ
    unfilledSegmentedXYZ=find(isnan(segmentedXYZ(:,1)));
    segmentedXYZ(unfilledSegmentedXYZ(1):unfilledSegmentedXYZ(1)+size(centroidsIN.segmented{z},1)-1,1:4)=horzcat(centroidsIN.segmented{z},z*ones(size((centroidsIN.segmented{z}),1),1),1*ones(size((centroidsIN.segmented{z}),1),1));
%     replacedXYZ(unfilledReplacedXYZ(1):unfilledReplacedXYZ(1)+size(centroidsIN.replaced{z},1)-1,9)=(1:size(centroidsIN.replaced{z},1))';
    segmentedXYZ(unfilledSegmentedXYZ(1):unfilledSegmentedXYZ(1)+size(centroidsIN.segmented{z},1)-1,9)=(1:size(centroidsIN.segmented{z},1))'; % redundant but don't leave as nan
    segmentedXYZ(unfilledSegmentedXYZ(1):unfilledSegmentedXYZ(1)+size(centroidsIN.segmented{z},1)-1,10)=areasIN.segmented{z};
    segmentedXYZ(unfilledSegmentedXYZ(1):unfilledSegmentedXYZ(1)+size(centroidsIN.segmented{z},1)-1,11)=pixflagsIN.segmented{z};
end

fullXYZ=segmentedXYZ;
[~,sortInds]=sort(fullXYZ(:,3));
fullXYZ=fullXYZ(sortInds,:);
nextGroupNum=max(fullXYZ(:,5))+1;
% ID vertical clusters
groupingInds=find(fullXYZ(:,3)>nonemptyZ(1));

for evalCentroid=groupingInds'
    thisCentroidZ=fullXYZ(evalCentroid,3);
    centroidsAboveInds=find(fullXYZ(:,3)==(thisCentroidZ-1));
    centroidsAbove=fullXYZ(centroidsAboveInds,:);
    XYdists=sqrt((centroidsAbove(:,1)-fullXYZ(evalCentroid,1)).^2+(centroidsAbove(:,2)-fullXYZ(evalCentroid,2)).^2);
    [minDist,minDistInd]=min(XYdists);
%     minDistZDiff=(thisCentroidZ-centroidsAbove(minDistInd,3))*zstepsize;
    minDistInd=centroidsAboveInds(minDistInd);  % index into fullXYZ instead
    if minDist<centroidDistCutoff           %%%%%%minDistZDiff>minDist
        if fullXYZ(minDistInd,7)>0
            competitorInd=find(fullXYZ(:,5)==fullXYZ(minDistInd,5) & fullXYZ(:,3)==thisCentroidZ);
            if length(competitorInd)>1
                disp(fullXYZ(competitorInd,:))
                competitorInd=competitorInd(1);
            end
            
            if minDist<fullXYZ(competitorInd,8) % if this point is closer, reassign competitor and add this one
                fullXYZ(competitorInd,5)=nextGroupNum;
                nextGroupNum=nextGroupNum+1;
                fullXYZ(competitorInd,8)=0;
                
                fullXYZ(evalCentroid,5)=fullXYZ(minDistInd,5);
                fullXYZ(evalCentroid,8)=minDist;
            
            else        % if competitor is closer, put this one in a new group
                fullXYZ(evalCentroid,5)=nextGroupNum;
                nextGroupNum=nextGroupNum+1;
                fullXYZ(evalCentroid,8)=minDist;                
            end
        else
            fullXYZ(evalCentroid,5)=fullXYZ(minDistInd,5);
            fullXYZ(minDistInd,7)=fullXYZ(minDistInd,7)+1;
            fullXYZ(evalCentroid,8)=minDist;
        end
    elseif any(fullXYZ(:,3)==(thisCentroidZ-2))  % try plane 2 above
        centroidsAboveInds=find(fullXYZ(:,3)==(thisCentroidZ-2));
        centroidsAbove=fullXYZ(centroidsAboveInds,:);
        XYdists=sqrt((centroidsAbove(:,1)-fullXYZ(evalCentroid,1)).^2+(centroidsAbove(:,2)-fullXYZ(evalCentroid,2)).^2);
        [minDist,minDistInd]=min(XYdists);
    %     minDistZDiff=(thisCentroidZ-centroidsAbove(minDistInd,3))*zstepsize;
        minDistInd=centroidsAboveInds(minDistInd);  % index into fullXYZ instead
        if minDist<centroidDistCutoff           %%%%%%minDistZDiff>minDist
            if fullXYZ(minDistInd,7)>0
                competitorInd=find(fullXYZ(:,5)==fullXYZ(minDistInd,5) & fullXYZ(:,3)==thisCentroidZ);
                if length(competitorInd)>1
                    disp(fullXYZ(competitorInd,:))
                    competitorInd=competitorInd(1);
                end

                if minDist<fullXYZ(competitorInd,8) % if this point is closer, reassign competitor and add this one
                    fullXYZ(competitorInd,5)=nextGroupNum;
                    nextGroupNum=nextGroupNum+1;
                    fullXYZ(competitorInd,8)=0;

                    fullXYZ(evalCentroid,5)=fullXYZ(minDistInd,5);
                    fullXYZ(evalCentroid,8)=minDist;

                else        % if competitor is closer, put this one in a new group
                    fullXYZ(evalCentroid,5)=nextGroupNum;
                    nextGroupNum=nextGroupNum+1;
                    fullXYZ(evalCentroid,8)=minDist;                
                end
            else
                fullXYZ(evalCentroid,5)=fullXYZ(minDistInd,5);
                fullXYZ(minDistInd,7)=fullXYZ(minDistInd,7)+1;
                fullXYZ(evalCentroid,8)=minDist;
            end
        else  % new centroid group
                fullXYZ(evalCentroid,5)=nextGroupNum;
                nextGroupNum=nextGroupNum+1;
                fullXYZ(evalCentroid,8)=minDist;
        end
    else  % new centroid group
        fullXYZ(evalCentroid,5)=nextGroupNum;
        nextGroupNum=nextGroupNum+1;
        fullXYZ(evalCentroid,8)=minDist;
    end
end

centroidsOUT=struct('twoD',[],'threeD',[],'nonemptyZ',[],'stereo',[]);
centroidsOUT.twoD=fullXYZ;


% find midpoint of vertical cluster
uniqueObjectNums=unique(centroidsOUT.twoD(:,5));
uniqueObjectNums=horzcat(uniqueObjectNums,nan(length(uniqueObjectNums),6));

for objNum=1:length(uniqueObjectNums(:,1))
    centroids2d=centroidsOUT.twoD(centroidsOUT.twoD(:,5)==uniqueObjectNums(objNum,1),:);
    centroids2d(:,10)=smooth(centroids2d(:,10));
    
    if size(centroids2d,1)==1
        uniqueObjectNums(objNum,5)=5;
        uniqueObjectNums(objNum,4)=centroids2d(1,3);
    elseif any(centroids2d(:,3)==1) || any(centroids2d(:,3)==max(nonemptyZ))% object is truncated by image stack boundary
        
        sizeChangeZeroCrossings=diff(sign(diff(centroids2d(:,10))));
        if ~any(sizeChangeZeroCrossings)  
            [~,maxsizeind]=max(centroids2d(:,10));
            if maxsizeind==min(nonemptyZ) || maxsizeind==max(nonemptyZ)   % object's widest point maybe outside stack boundary
                uniqueObjectNums(objNum,4)=centroids2d(maxsizeind,3);
                uniqueObjectNums(objNum,5)=0;
            else               
                if centroids2d(maxsizeind,3)==min(centroids2d(:,3))
                    uniqueObjectNums(objNum,5)=1;
                    uniqueObjectNums(objNum,4)=centroids2d(maxsizeind,3);
                elseif centroids2d(maxsizeind,3)==max(centroids2d(:,3))
                    uniqueObjectNums(objNum,5)=2;
                    uniqueObjectNums(objNum,4)=centroids2d(maxsizeind,3);
                else
                    uniqueObjectNums(objNum,4)=sum(centroids2d(:,10).*centroids2d(:,3))/sum(centroids2d(:,10));  %% not sure what to do in this case
                end
            end
        elseif length(find(sizeChangeZeroCrossings))<2
            [~,maxsizeind]=max(centroids2d(:,10));
            flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
            centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
            uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
            uniqueObjectNums(objNum,5)=3;
        elseif length(find(sizeChangeZeroCrossings))==2
            SCZCinds=find(sizeChangeZeroCrossings);
            if SCZCinds(2)==SCZCinds(1)+1
                if centroids2d(SCZCinds(2),10)==centroids2d(SCZCinds(2)+1,10)
                    [~,maxsizeind]=max(centroids2d(:,10));
                    flankLength=min([size(centroids2d,1)-(maxsizeind+1),(maxsizeind-1)]);
                    centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength+1),:);
                    uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                    uniqueObjectNums(objNum,5)=3;
                else 
                    [~,maxsizeind]=max(centroids2d(:,10));
                    flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
                    centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
                    uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                    uniqueObjectNums(objNum,5)=4;
                end
            else
                [~,maxsizeind]=max(centroids2d(:,10));
                flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
                centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
                uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                uniqueObjectNums(objNum,5)=4;
            end
        else
            [~,maxsizeind]=max(centroids2d(:,10));
            flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
            centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
            uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
            uniqueObjectNums(objNum,5)=4;
        end
    
    else   % entire object contained within image stack
        sizeChangeZeroCrossings=diff(sign(diff(centroids2d(:,10))));
        if ~any(sizeChangeZeroCrossings)
            [~,maxsizeind]=max(centroids2d(:,10));
            if centroids2d(maxsizeind,3)==min(centroids2d(:,3))
                uniqueObjectNums(objNum,4)=min(centroids2d(:,3)); 
                uniqueObjectNums(objNum,5)=1;
            elseif centroids2d(maxsizeind,3)==max(centroids2d(:,3))
                uniqueObjectNums(objNum,4)=max(centroids2d(:,3));
                uniqueObjectNums(objNum,5)=2;
            else
                flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
                centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
                uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
            end
        elseif length(find(sizeChangeZeroCrossings))<2
            [~,maxsizeind]=max(centroids2d(:,10));
            flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
            centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
            uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
            uniqueObjectNums(objNum,5)=3;
        elseif length(find(sizeChangeZeroCrossings))==2
            SCZCinds=find(sizeChangeZeroCrossings);
            if SCZCinds(2)==SCZCinds(1)+1
                if centroids2d(SCZCinds(2),10)==centroids2d(SCZCinds(2)+1,10)
                    [~,maxsizeind]=max(centroids2d(:,10));
                    flankLength=min([size(centroids2d,1)-(maxsizeind+1),(maxsizeind-1)]);
                    centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength+1),:);
                    uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                    uniqueObjectNums(objNum,5)=3;
                else 
                    [~,maxsizeind]=max(centroids2d(:,10));
                    flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
                    centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
                    uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                    uniqueObjectNums(objNum,5)=4;
                end
            else
                [~,maxsizeind]=max(centroids2d(:,10));
                flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
                centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
                uniqueObjectNums(objNum,4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                uniqueObjectNums(objNum,5)=4;
            end
        else
            uniqueObjectNums(objNum,4)=sum(centroids2d(:,10).*centroids2d(:,3))/sum(centroids2d(:,10)); 
            uniqueObjectNums(objNum,5)=3;
        end
    end
    if isnan(uniqueObjectNums(objNum,5))
        error
    end
    if isnan(uniqueObjectNums(objNum,4))
        uniqueObjectNums(objNum,4)=sum(centroids2d(:,10).*centroids2d(:,3))/sum(centroids2d(:,10));
    end
    uniqueObjectNums(objNum,2)=median(centroids2d(:,1));
    uniqueObjectNums(objNum,3)=median(centroids2d(:,2));
    uniqueObjectNums(objNum,6)=max(centroids2d(:,11));  % pixflag
    uniqueObjectNums(objNum,7)=max(centroids2d(:,10));  % area

end

if zStepSize/pix_um_Conversion<1      %if round(meanParticleDiam/zStepSize/8)>1
uniqueObjectNums=uniqueObjectNums(uniqueObjectNums(:,5)<5,:);   % delete singletons
end

% size(uniqueObjectNums)
for rpts=1:3
    uniqueObjectNumsCopy=uniqueObjectNums;
    objSkipList=[];
    for objNum=1:size(uniqueObjectNumsCopy(:,1),1)     % identify objects that are too close in 3d and merge them
        if ~ismember(uniqueObjectNumsCopy(objNum),objSkipList)
            inds=1:size(uniqueObjectNums(:,1),1);
            delInd=find(uniqueObjectNums(:,1)==uniqueObjectNumsCopy(objNum));
            inds=setdiff(inds,delInd)';
            if uniqueObjectNumsCopy(objNum,5)==1 

                distsXYZ=sqrt((uniqueObjectNums(inds,2)-uniqueObjectNumsCopy(objNum,2)).^2+(uniqueObjectNums(inds,3)-uniqueObjectNumsCopy(objNum,3)).^2+((uniqueObjectNums(inds,4)-uniqueObjectNumsCopy(objNum,4))*zStepSize).^2);  
                closeXYZObjInds=inds(distsXYZ<centroidDistCutoff3d);
                closeXYZObjInds=closeXYZObjInds(find(uniqueObjectNums(closeXYZObjInds,4)-uniqueObjectNumsCopy(objNum,4)>=0));
  
                distsXY=sqrt((uniqueObjectNums(inds,2)-uniqueObjectNumsCopy(objNum,2)).^2+(uniqueObjectNums(inds,3)-uniqueObjectNumsCopy(objNum,3)).^2);
                closeXYObjInds=inds(distsXY<centroidDistCutoff3dXY);
                closeZDists=(uniqueObjectNums(closeXYObjInds,4)-uniqueObjectNumsCopy(objNum,4))*zStepSize;
                closeXYObjInds=closeXYObjInds(find((closeZDists<centroidDistCutoff3dZ).*(closeZDists>=0)));
            elseif uniqueObjectNumsCopy(objNum,5)==2
                distsXYZ=sqrt((uniqueObjectNums(inds,2)-uniqueObjectNumsCopy(objNum,2)).^2+(uniqueObjectNums(inds,3)-uniqueObjectNumsCopy(objNum,3)).^2+((uniqueObjectNums(inds,4)-uniqueObjectNumsCopy(objNum,4))*zStepSize).^2);  
                closeXYZObjInds=inds(distsXYZ<centroidDistCutoff3d);
                closeXYZObjInds=closeXYZObjInds(find(uniqueObjectNumsCopy(objNum,4)-uniqueObjectNums(closeXYZObjInds,4)>=0));
  
                distsXY=sqrt((uniqueObjectNums(inds,2)-uniqueObjectNumsCopy(objNum,2)).^2+(uniqueObjectNums(inds,3)-uniqueObjectNumsCopy(objNum,3)).^2);
                closeXYObjInds=inds(distsXY<centroidDistCutoff3dXY);
                closeZDists=(uniqueObjectNumsCopy(objNum,4)-uniqueObjectNums(closeXYObjInds,4))*zStepSize;
                closeXYObjInds=closeXYObjInds(find((closeZDists<centroidDistCutoff3dZ).*(closeZDists>=0)));
            else
                closeXYObjInds=[];
                closeXYZObjInds=[];
            end
            closeObjInds=unique(vertcat(closeXYZObjInds,closeXYObjInds));   
            if ~isempty(closeObjInds) 
                tagsToGroup = getTagsToGroup(closeObjInds,objNum,uniqueObjectNums,uniqueObjectNumsCopy,centroidsOUT,nonemptyZ,zStepSize);
                if ~isempty(tagsToGroup)
                    objSkipList=vertcat(objSkipList,tagsToGroup);
                    [uniqueObjectNums,centroidsOUT] = combineTaggedGroups(tagsToGroup,uniqueObjectNums,centroidsOUT,nonemptyZ);
                end
            end
        end
    end
end

for rpts=1:3
    uniqueObjectNumsCopy=uniqueObjectNums;
    objSkipList=[];
    for objNum=1:size(uniqueObjectNumsCopy(:,1),1)     % identify objects that are too close in 3d and merge them
        if ~ismember(uniqueObjectNumsCopy(objNum),objSkipList)
            inds=1:size(uniqueObjectNums(:,1),1);
            delInd=find(uniqueObjectNums(:,1)==uniqueObjectNumsCopy(objNum));
            inds=setdiff(inds,delInd)';
            distsXYZ=sqrt((uniqueObjectNums(inds,2)-uniqueObjectNumsCopy(objNum,2)).^2+(uniqueObjectNums(inds,3)-uniqueObjectNumsCopy(objNum,3)).^2+((uniqueObjectNums(inds,4)-uniqueObjectNumsCopy(objNum,4))*zStepSize).^2);  
            closeObjInds=inds(distsXYZ<centroidDistCutoff3d_2);     
            if ~isempty(closeObjInds) 
                dists=sqrt((uniqueObjectNums(closeObjInds,2)-uniqueObjectNumsCopy(objNum,2)).^2+(uniqueObjectNums(closeObjInds,3)-uniqueObjectNumsCopy(objNum,3)).^2+((uniqueObjectNums(closeObjInds,4)-uniqueObjectNumsCopy(objNum,4))*zStepSize).^2);  
                [~,mindistind]=min(dists);
                
                tagsToGroup = getTagsToGroup(closeObjInds(mindistind),objNum,uniqueObjectNums,uniqueObjectNumsCopy,centroidsOUT,nonemptyZ,zStepSize);
                if ~isempty(tagsToGroup)
                    objSkipList=vertcat(objSkipList,tagsToGroup);
                    [uniqueObjectNums,centroidsOUT] = combineTaggedGroups(tagsToGroup,uniqueObjectNums,centroidsOUT,nonemptyZ);
                end
                
            end
        end
    end

end

for rpts=1:3
    uniqueObjectNumsCopy=uniqueObjectNums;
    objSkipList=[];
    for objNum=1:size(uniqueObjectNumsCopy(:,1),1)     % identify objects that are too close in 3d and merge them
        if ~ismember(uniqueObjectNumsCopy(objNum),objSkipList)
            inds=1:size(uniqueObjectNums(:,1),1);
            delInd=find(uniqueObjectNums(:,1)==uniqueObjectNumsCopy(objNum));
            inds=setdiff(inds,delInd)';
            distsXY=sqrt((uniqueObjectNums(inds,2)-uniqueObjectNumsCopy(objNum,2)).^2+(uniqueObjectNums(inds,3)-uniqueObjectNumsCopy(objNum,3)).^2);
            closeXYObjInds=inds(distsXY<centroidDistCutoff3dXY_2);
            closeZDists=abs(uniqueObjectNums(closeXYObjInds,4)-uniqueObjectNumsCopy(objNum,4));
            closeObjInds=closeXYObjInds(closeZDists<=centroidDistCutoff3dZ_2);     
            if ~isempty(closeObjInds)        
                tagsToGroup=vertcat(uniqueObjectNums(closeObjInds,1),uniqueObjectNumsCopy(objNum,1));
                objSkipList=vertcat(objSkipList,tagsToGroup);
                [uniqueObjectNums,centroidsOUT] = combineTaggedGroups(tagsToGroup,uniqueObjectNums,centroidsOUT,nonemptyZ);
            end
        end
    end
end


% eliminate erroneous centroids: singles (no neighbor matches), clusters entirely of "replaced" centroids
keepObjectTag=ones(size(centroidsOUT.twoD,1),1);

for entryNum = 1:size(centroidsOUT.twoD,1)
    points=centroidsOUT.twoD(centroidsOUT.twoD(:,5)==centroidsOUT.twoD(entryNum,5),:);
    if min(points(:,3))~=min(nonemptyZ) && max(points(:,3))~=max(nonemptyZ)
        if (max(points(:,3))-min(points(:,3)))*zStepSize<meanParticleDiam/4 && zStepSize/pix_um_Conversion<1      %round(meanParticleDiam/zStepSize/8)>1
            keepObjectTag(entryNum)=0;
        elseif size(points,1)<1.5/(zStepSize/pix_um_Conversion) && zStepSize/pix_um_Conversion<1      %round(meanParticleDiam/zStepSize/8)>1  % ~3 for 0.5 step size
            keepObjectTag(entryNum)=0;
        end
    elseif size(points,1)<2 && zStepSize/pix_um_Conversion<1      %round(meanParticleDiam/zStepSize/8)>1
        keepObjectTag(entryNum)=0;
    end    
end
keepObjectTags=centroidsOUT.twoD(find(keepObjectTag),5);
[~,inds]=setdiff(uniqueObjectNums(:,1),keepObjectTags);
excludeRows=uniqueObjectNums(inds,:);
uniqueObjectNums=setdiff(uniqueObjectNums,excludeRows,'rows');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excludeRows=[];

if ~isempty(excludeRows)
    uniqueObjectNums=setdiff(uniqueObjectNums,excludeRows,'rows');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%

centroidsOUT.threeD=uniqueObjectNums;
centroidsOUT.nonemptyZ=nonemptyZ;

if threeDstereoflag
    centroidsOUT.stereo=applyExclusion(centroidsOUT.threeD,twoDstereoflag,threeDstereoflag,exclusionLimits);
else
    centroidsOUT.stereo=applyExclusion(centroidsOUT.threeD,twoDstereoflag,threeDstereoflag,[]);
end

function tagsToGroup = getTagsToGroup(candidates,objNum,uniqueObjectNums,uniqueObjectNumsCopy,centroidsOUT,nonemptyZ,zStepSize)

if length(candidates)>1
    disp('')
end

candidates=horzcat(candidates,nan(length(candidates),3));
origObjType=uniqueObjectNumsCopy(objNum,5);
for candNum=1:size(candidates,1)
    candidates(candNum,2)=uniqueObjectNums(candidates(candNum),5);
    minTagNum=min(vertcat(uniqueObjectNums(candidates(candNum),1),uniqueObjectNumsCopy(objNum,1)));
    [temp,~] = combineTaggedGroups(vertcat(uniqueObjectNums(candidates(candNum),1),uniqueObjectNumsCopy(objNum,1)),uniqueObjectNums,centroidsOUT,nonemptyZ);   
    candidates(candNum,3)=temp(temp(:,1)==minTagNum,5);
    
    switch candidates(candNum,3)
        case {3,0}   % makes a complete object, or widest at stack bound
            candidates(candNum,4)=1;            
        case {1,2}  % 1 bottom, 2 top
            if ismember(candidates(candNum,2),[1 2 5]) && ismember(origObjType,[1 2 5])
                candidates(candNum,4)=1;
            else
                candidates(candNum,4)=0;
            end
        case 4  % multiple zero crossings
            if candidates(candNum,2)==4 && origObjType==4
                candidates(candNum,4)=1;  
            else
                candidates(candNum,4)=0;      
            end
    end
end
candidates=candidates(candidates(:,4)==1,1);

if isempty(candidates)
    tagsToGroup=[];
else
    dists=sqrt((uniqueObjectNums(candidates,2)-uniqueObjectNumsCopy(objNum,2)).^2+(uniqueObjectNums(candidates,3)-uniqueObjectNumsCopy(objNum,3)).^2+((uniqueObjectNums(candidates,4)-uniqueObjectNumsCopy(objNum,4))*zStepSize).^2);  
    [~,mindistind]=min(dists);
    tagsToGroup=vertcat(uniqueObjectNums(candidates(mindistind),1),uniqueObjectNumsCopy(objNum,1));
end


function [uniqueObjectNums,centroidsOUT] = combineTaggedGroups(tagsToGroup,uniqueObjectNums,centroidsOUT,nonemptyZ)
centroidsToGroup=[];
for ind=tagsToGroup'
    centroidsToGroup=vertcat(centroidsToGroup,centroidsOUT.twoD(centroidsOUT.twoD(:,5)==ind,:));
    uniqueObjectNums=setdiff(uniqueObjectNums,uniqueObjectNums(uniqueObjectNums(:,1)==ind,:),'rows');
end
centroidsOUT.twoD=setdiff(centroidsOUT.twoD,centroidsToGroup,'rows');
centroidsToGroup(:,5)=min(centroidsToGroup(:,5));
centroidsOUT.twoD=vertcat(centroidsOUT.twoD,centroidsToGroup);  

uniqueZ=unique(centroidsToGroup(:,3));
centroids2d=horzcat(nan(length(uniqueZ),2),uniqueZ,nan(length(uniqueZ),8));
for zind=1:length(uniqueZ)
    thisZ=centroidsToGroup(centroidsToGroup(:,3)==uniqueZ(zind),:);
    centroids2d(zind,1)=mean(thisZ(:,1));
    centroids2d(zind,2)=mean(thisZ(:,2));
    centroids2d(zind,5)=centroidsToGroup(1,5);
    centroids2d(zind,10)=sum(thisZ(:,10));
    centroids2d(zind,11)=max(thisZ(:,11));
end

centroid3dinfo=nan(1,7);
[~,sortzind]=sort(centroids2d(:,3),'ascend');
centroids2d=centroids2d(sortzind,:);
centroids2d(:,10)=smooth(centroids2d(:,10));
if size(centroids2d,1)==1
    centroid3dinfo(5)=5;
    centroid3dinfo(4)=centroids2d(1,3);
elseif any(centroids2d(:,3)==1) || any(centroids2d(:,3)==max(nonemptyZ))% object is truncated by image stack boundary
    sizeChangeZeroCrossings=diff(sign(diff(centroids2d(:,10))));
    if ~any(sizeChangeZeroCrossings)  
        [~,maxsizeind]=max(centroids2d(:,10));
        if maxsizeind==min(nonemptyZ) || maxsizeind==max(nonemptyZ)   % object's widest point maybe outside stack boundary
            centroid3dinfo(4)=centroids2d(maxsizeind,3);
            centroid3dinfo(5)=0;
        else
            centroid3dinfo(4)=sum(centroids2d(:,10).*centroids2d(:,3))/sum(centroids2d(:,10));  %% not sure what to do in this case
            if centroids2d(maxsizeind,3)==min(centroids2d(:,3))
                centroid3dinfo(5)=1;
            elseif centroids2d(maxsizeind,3)==max(centroids2d(:,3))
                centroid3dinfo(5)=2;
            end
        end
    elseif length(find(sizeChangeZeroCrossings))<2
        [~,maxsizeind]=max(centroids2d(:,10));
        flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
        centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
        centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
        centroid3dinfo(5)=3;
    elseif length(find(sizeChangeZeroCrossings))==2
        SCZCinds=find(sizeChangeZeroCrossings);
        if SCZCinds(2)==SCZCinds(1)+1
            if centroids2d(SCZCinds(2),10)==centroids2d(SCZCinds(2)+1,10)
                [~,maxsizeind]=max(centroids2d(:,10));
                flankLength=min([size(centroids2d,1)-(maxsizeind+1),(maxsizeind-1)]);
                centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength+1),:);
                centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                centroid3dinfo(5)=3;
            else 
                [~,maxsizeind]=max(centroids2d(:,10));
                flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
                centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
                centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                centroid3dinfo(5)=4;
            end
        else
            [~,maxsizeind]=max(centroids2d(:,10));
            flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
            centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
            centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
            centroid3dinfo(5)=4;
        end
    else
        [~,maxsizeind]=max(centroids2d(:,10));
        flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
        centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
        centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
        centroid3dinfo(5)=4;
    end

else   % entire object contained within image stack
    sizeChangeZeroCrossings=diff(sign(diff(centroids2d(:,10))));
    if ~any(sizeChangeZeroCrossings)
        [~,maxsizeind]=max(centroids2d(:,10));
        if centroids2d(maxsizeind,3)==min(centroids2d(:,3))
            centroid3dinfo(4)=min(centroids2d(:,3)); 
            centroid3dinfo(5)=1;
        elseif centroids2d(maxsizeind,3)==max(centroids2d(:,3))
            centroid3dinfo(4)=max(centroids2d(:,3));
            centroid3dinfo(5)=2;
        else
            flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
            centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
            centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
        end
    elseif length(find(sizeChangeZeroCrossings))<2
        [~,maxsizeind]=max(centroids2d(:,10));
        flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
        centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
        centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
        centroid3dinfo(5)=3;
    elseif length(find(sizeChangeZeroCrossings))==2
        SCZCinds=find(sizeChangeZeroCrossings);
        if SCZCinds(2)==SCZCinds(1)+1
            if centroids2d(SCZCinds(2),10)==centroids2d(SCZCinds(2)+1,10)
                [~,maxsizeind]=max(centroids2d(:,10));
                flankLength=min([size(centroids2d,1)-(maxsizeind+1),(maxsizeind-1)]);
                centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength+1),:);
                centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                centroid3dinfo(5)=3;
            else 
                [~,maxsizeind]=max(centroids2d(:,10));
                flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
                centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
                centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
                centroid3dinfo(5)=4;
            end
        else
            [~,maxsizeind]=max(centroids2d(:,10));
            flankLength=min([size(centroids2d,1)-maxsizeind,maxsizeind-1]);
            centeredSeries=centroids2d((maxsizeind-flankLength):(maxsizeind+flankLength),:);
            centroid3dinfo(4)=sum(centeredSeries(:,10).*centeredSeries(:,3))/sum(centeredSeries(:,10)); 
            centroid3dinfo(5)=4;
        end
    else
        centroid3dinfo(4)=sum(centroids2d(:,10).*centroids2d(:,3))/sum(centroids2d(:,10)); 
        centroid3dinfo(5)=3;
    end
    
end

centroid3dinfo(1)=centroids2d(1,5);
centroid3dinfo(2)=median(centroids2d(:,1));
centroid3dinfo(3)=median(centroids2d(:,2));
centroid3dinfo(6)=max(centroids2d(:,11));
centroid3dinfo(7)=max(centroids2d(:,10));
if isnan(centroid3dinfo(4))
    centroid3dinfo(4)=sum(centroids2d(:,10).*centroids2d(:,3))/sum(centroids2d(:,10));
end
uniqueObjectNums=vertcat(uniqueObjectNums,centroid3dinfo);

function stereoCentroids=applyExclusion(threeDCentroids,twoDstereoflag,threeDstereoflag,exclusionLimits)

if threeDstereoflag
    excludeInds=find(threeDCentroids(:,4)>=exclusionLimits(2));
    excludeInds=vertcat(excludeInds,find(threeDCentroids(:,4)<exclusionLimits(1)));
    if twoDstereoflag
        excludeInds=vertcat(excludeInds,find(threeDCentroids(:,6)));
    end
elseif twoDstereoflag
    excludeInds=find(threeDCentroids(:,6));
else
    excludeInds=[];
end
excludeRows=threeDCentroids(excludeInds,:);
if ~isempty(excludeRows)
    stereoCentroids=setdiff(threeDCentroids,excludeRows,'rows');
else
    stereoCentroids=[];
end
