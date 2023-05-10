function [imagesOut,centroidsOut,areasOut,pixflagsOut] = segmentObjectsByClustering(filename,zsteps,threshznum,zStepSize,printflag)
% updated 12/09/14 
% requires bradley.m and averagefilter.m by Jan Motl shared on Matlab
% Central file exchange: https://www.mathworks.com/matlabcentral/fileexchange/40854-bradley-local-image-thresholding

global zstepsize
zstepsize=zStepSize;
tic;
imagesOut=struct([]);
centroidsOut=struct('original',[],'segmented',[],'replaced',[]);pixflagsOut=struct('original',[],'segmented',[],'replaced',[]);areasOut=struct('original',[],'segmented',[],'replaced',[]);
imgSeries=struct([]);
centroidsForComparison=struct([]);allAreas=struct([]);
    
img=double(imread(filename,threshznum));
img=medfilt2(img);
lowPlaneMinLevel=0.5*graythresh(img/max(max(img))); 
numZ=numel(imfinfo(filename));
se=strel('disk',2);
for z = 1:numZ 
    img=double(imread(filename,z));
    img=medfilt2(img);
    if exist('lowPlaneMinLevel','var')
        img=logical(bradley(img.*im2bw(img/max(max(img)),max(0.5*graythresh(img/max(max(img))),lowPlaneMinLevel)),[600 600],-70,'symmetric'));
    else
        img=logical(bradley(img.*im2bw(img/max(max(img)),0.5*graythresh(img/max(max(img)))),[600 600],-70,'symmetric'));
    end
    img=imopen(imerode(img,se),se); 
    img=imfill(img,'holes');
    imgSeries{z}=img;

	objProps  = regionprops(img,'centroid','area');
    areas=cat(1,objProps.Area);
    objProps=objProps(areas>100);  % for nuclei 63x
    centroidsForComparison{z}=cat(1,objProps.Centroid);
    allAreas{z}=cat(1,objProps.Area);
    
    if z==1
        imgDims=size(img);
    end
end
medArea=nanmedian(cellfun(@nanmedian,allAreas(zsteps)));
for z = zsteps 
   % read in image, preprocess, make working copy
    imgCopy=imgSeries{z};
    objProps  = regionprops(imgCopy,'centroid','area','solidity','PixelList');
    areas=cat(1,objProps.Area);
    objProps=objProps(areas>100);  % for nuclei 63x
%     objProps=objProps(areas>20);  % for vglut2 
    areas=cat(1,objProps.Area);
    objProps=objProps(areas<50*medArea);
	 % select objects to segment
    objList=getObjectsToSegment(objProps,medArea);
    segmentedInds=nan(length(objList),1);
    segmentedCentroids=nan(5*length(objList),2);
    segmentedPixFlags=nan(5*length(objList),1);
    segmentedAreas=nan(5*length(objList),1);
    counter=1;
   % segment each
    imgInfo = struct('img',{[]},'nZSteps',[],'currentZ',[]);
    imgInfo.img=imgCopy;
    imgInfo.medArea=medArea; %imgInfo.medArea=median(areas(areas>200));
    imgInfo.dims=imgDims;
    for object=1:length(objList)
        objectInd=objList(object);
        pixels=single(objProps(objectInd).PixelList);
        firstEstK=round(objProps(objectInd).Area/imgInfo.medArea);
        clusterInfo=clusterPixels(pixels,firstEstK,objProps(objectInd).Centroid,imgInfo,centroidsForComparison,z);  % output: clusterInfo ->   .k   .centroids  .area   .pixflags
        if clusterInfo.k>1
            segmentedInds(object)=objectInd;
            segmentedCentroids(counter:counter+clusterInfo.k-1,:)=clusterInfo.centroids;
            segmentedPixFlags(counter:counter+clusterInfo.k-1)=clusterInfo.pixflags; 
            segmentedAreas(counter:counter+clusterInfo.k-1)=clusterInfo.areas; 
            counter=counter+clusterInfo.k;
        end
    end
    
   
    segmentedInds=segmentedInds(~isnan(segmentedInds));
    segmentedCentroids=segmentedCentroids(1:counter-1,:);
    segmentedPixFlags=segmentedPixFlags(1:counter-1);
    segmentedAreas=segmentedAreas(1:counter-1);
   
    centroids=cat(1,objProps.Centroid);
    areas=cat(1,objProps.Area);
    originalPixFlags=zeros(size(centroids,1),1);
    for object=1:size(centroids,1)
        pixels=objProps(object).PixelList;
        if any(pixels(:,1)==imgDims(1))
            originalPixFlags(object)=1;
        elseif any(pixels(:,2)==imgDims(2))
            originalPixFlags(object)=1;
        end
    end
     
    centroidIndsToKeep=setdiff((1:length(centroids))',segmentedInds);
    newCentroids=vertcat(centroids(centroidIndsToKeep,:),segmentedCentroids);
    centroidsForComparison{z}=newCentroids;
    newPixFlags=vertcat(originalPixFlags(centroidIndsToKeep,:),segmentedPixFlags);
    newAreas=vertcat(areas(centroidIndsToKeep,:),segmentedAreas);
    imagesOut{z}=imgCopy;
    if printflag
        figure();imshow(imgCopy);
        hold on;
        plot(centroids(:,1),centroids(:,2),'r.')  % red points that don't get plotted over in blue were replaced with set of new centroids
        plot(newCentroids(:,1),newCentroids(:,2),'b.')
        plot(centroids(objList,1),centroids(objList,2),'m.')
        plot(segmentedCentroids(:,1),segmentedCentroids(:,2),'g.')  % replaced centroids
    end
    
    centroidsOut.original{z}=centroids;
    pixflagsOut.original{z}=originalPixFlags;
    areasOut.original{z}=areas;
    centroidsOut.segmented{z}=newCentroids;
    pixflagsOut.segmented{z}=newPixFlags;
    areasOut.segmented{z}=newAreas;
    centroidsOut.replaced{z}=centroids(segmentedInds,:);
    pixflagsOut.replaced{z}=originalPixFlags(segmentedInds,:);
    areasOut.replaced{z}=areas(segmentedInds,:);
    toc
end
% toc

function [objListOut]=getObjectsToSegment(objPropsIn,medArea)

areas=cat(1,objPropsIn.Area);
areaCutoff=medArea; % areaCutoff=median(areas);
indsA=find(areas>areaCutoff);
objListOut=indsA;

function clusterInfoOut=clusterPixels(pixels,firstEstK,centroidIn,imgInfoIn,centroidsForComparisonIn,currentZ)   % output: clusterInfo ->   .k   .centroids   .area   .pixflags

vectorDims=size(pixels);
if vectorDims(1)<vectorDims(2)      % separate X and Y positions by column
    pixels=pixels';
end

clusterInfoOut = struct('k',[],'centroids',[]);

idx=struct([]);centroids=struct([]);clusterpixels=struct([]);objAreas=struct([]);hullAreaRatios=struct([]);
idx{1}=[];centroids{1}=[]; clusterpixels{1}=pixels;

imgFromPix=makeImFromPix(pixels,imgInfoIn.img);
objInfo=regionprops(imgFromPix,'Area','ConvexArea');
objAreas{1}=objInfo.Area;
hullAreaRatios{1}=objInfo.Area/objInfo.ConvexArea;
edgepoints=getEdgePoints(pixels,1);
elbowPoints=getElbowPoints(pixels,edgepoints);

if ~isempty(elbowPoints)

    kTestSeries=[firstEstK (firstEstK+1) (firstEstK-1)];
    kTestSeries=kTestSeries(kTestSeries>1);
    
    centroidMatchScore=nan(1,length(kTestSeries)+1);
    [centroidMatchScore(1),nearestCentroid]=getCentroidMatchScore(centroidIn,centroidsForComparisonIn,currentZ);  

    for kind=1:length(kTestSeries)
        k=kTestSeries(kind);
        [idx{k},centroids{k}] = kmeans(pixels,k, 'start','cluster','replicates',5,'EmptyAction','singleton');
        hullAreaRatios{k}=nan(1,k);
        objAreas{k}=nan(1,k);
        for cluster=1:k
            clusterpixels{k}{cluster}=pixels(idx{k}==cluster,:);
            imgFromPix=makeImFromPix(clusterpixels{k}{cluster},imgInfoIn.img);
            objInfo=regionprops(imgFromPix,'Area','ConvexArea');
            objAreas{k}(cluster)=objInfo.Area;
            hullAreaRatios{k}(cluster)=objInfo.Area/objInfo.ConvexArea;   
            centroidMatchScore(kind+1)=getCentroidMatchScore(centroids{k},centroidsForComparisonIn,currentZ,nearestCentroid);  
        end
    end
    
    elbowscore=nan(1,length(kTestSeries)+1);
    elbowscore(1)=1;
    circleEquivalentScore=nan(1,length(kTestSeries)+1);
    circleEquivalentScore(1)=overlappingCircleEquivalent(1,pixels,centroidIn,objAreas{1});   
    for kind=1:length(kTestSeries)
        k=kTestSeries(kind);
        circleEquivalentScore(kind+1)=overlappingCircleEquivalent(k,pixels,centroids{k},objAreas{k}); %circleEquivalentScore(k)=overlappingCircleEquivalent(k,pixels,centroids{k},objAreas{k},covMs{k},s{k},v{k});

        if k>2

            diffmat=centroids{k};
            centroidsmat=[];
            for counter=1:k-1
                targetpoint=diffmat(1,:);
                newdiffmat=setdiff(diffmat,targetpoint,'rows');
                mattocat=horzcat(repmat(targetpoint,size(newdiffmat,1),1),newdiffmat);
                centroidsmat=vertcat(centroidsmat,mattocat);
                diffmat=newdiffmat;
            end
            centroidsmat=horzcat(centroidsmat,nan(sum(1:k-1),1));
            centroidsmat(:,5)=sqrt((centroidsmat(:,1)-centroidsmat(:,3)).^2+(centroidsmat(:,2)-centroidsmat(:,4)).^2);
            [~,longestedge]=max(centroidsmat(:,5));
            centroidsmat=setdiff(centroidsmat,centroidsmat(longestedge,:),'rows');
            elbowscores=nan(sum(1:k-1)-1,1);

        else    % k==2
            centroidsmat=horzcat(centroids{k}(1,:),centroids{k}(2,:),nan(1));
            elbowscores=nan(1);
        end           

        for row=1:size(centroidsmat,1)
            elbowPointDists=nan(size(elbowPoints,1),1);
            testVal=nan(size(elbowPoints,1),1);
            centroidDist=sqrt((centroidsmat(row,1)-centroidsmat(row,3))^2+(centroidsmat(row,2)-centroidsmat(row,4))^2);
            for elbow=1:size(elbowPoints,1)
                elbowPointDists(elbow)= abs(det([(centroidsmat(row,1:2)-centroidsmat(row,3:4))' (centroidsmat(row,3:4)-elbowPoints(elbow,1:2))']))/sqrt((centroidsmat(row,1)-centroidsmat(row,3))^2+(centroidsmat(row,2)-centroidsmat(row,4))^2);        
                centroidElbowDist1=sqrt((centroidsmat(row,1)-elbowPoints(elbow,1))^2+(centroidsmat(row,2)-elbowPoints(elbow,2))^2);
                centroidElbowDist2=sqrt((centroidsmat(row,3)-elbowPoints(elbow,1))^2+(centroidsmat(row,4)-elbowPoints(elbow,2))^2);
                Ya = elbowPointDists(elbow)/tan(asin(elbowPointDists(elbow)/centroidElbowDist1));
                Yb = elbowPointDists(elbow)/tan(asin(elbowPointDists(elbow)/centroidElbowDist2));
                if elbowPointDists(elbow)<centroidDist
                    testVal(elbow)=abs(centroidDist-Ya-Yb); %testVal(elbow)=abs(centroidDist-Ya-Yb)<0.00001;        
                else
                    testVal(elbow)=1;
                end
            end      
            testVal=sort(testVal,'ascend');
            bestElbowTestVal=testVal(1); 
            elbowscores(row)=bestElbowTestVal<0.0001;
        end

        elbowscore(kind+1)=floor(min(elbowscores));  %elbowscore(k)=sum(elbowscores(row));

    end

    circleEquivalentScore=circleEquivalentScore./max(circleEquivalentScore);
    areaScore=abs(objAreas{1}./(imgInfoIn.medArea*([1 kTestSeries]))-1);
    areaScore=areaScore./max(areaScore);
    centroidMatchScore=centroidMatchScore./max(centroidMatchScore);

    switch length(kTestSeries)
        case 3
%             score=[1 1 1 1]./(circleEquivalentScore).^2.*elbowscore./(areaScore).^1./centroidMatchScore.^1;
            score=[hullAreaRatios{1} prod(hullAreaRatios{kTestSeries(1)}) prod(hullAreaRatios{kTestSeries(2)}) prod(hullAreaRatios{kTestSeries(3)})].^1./(circleEquivalentScore).^2.*elbowscore./(areaScore).^1./centroidMatchScore.^1;
        case 2
%             score=[1 1 1]./(circleEquivalentScore).^2.*elbowscore./(areaScore).^1./centroidMatchScore.^1;
            score=[hullAreaRatios{1} prod(hullAreaRatios{kTestSeries(1)}) prod(hullAreaRatios{kTestSeries(2)})].^1./(circleEquivalentScore).^2.*elbowscore./(areaScore).^1./centroidMatchScore.^1;
        case 1
%             score=[1 1]./(circleEquivalentScore).^2.*elbowscore./(areaScore).^1./centroidMatchScore.^1;
            score=[hullAreaRatios{1} prod(hullAreaRatios{kTestSeries(1)})].^1./(circleEquivalentScore).^2.*elbowscore./(areaScore).^1./centroidMatchScore.^1;
        case 0
            score=1;
    end
    if firstEstK>3
        score(1)=-inf;
    end
    [~,bestKInd]=max(score);
    kTestSeriesFull=horzcat(1,kTestSeries);
    clusterInfoOut.k=kTestSeriesFull(bestKInd);
%     end 
else
    clusterInfoOut.k=1;
end
    
clusterInfoOut.centroids=centroids{clusterInfoOut.k};
clusterInfoOut.areas=objAreas{clusterInfoOut.k};

pixFlags=zeros(clusterInfoOut.k,1);
for clusternum=1:length(pixFlags)
    if clusterInfoOut.k==1
        pixelX=clusterpixels{clusterInfoOut.k}(:,1);pixelY=clusterpixels{clusterInfoOut.k}(:,2);
    else
        pixelX=clusterpixels{clusterInfoOut.k}{clusternum}(:,1);pixelY=clusterpixels{clusterInfoOut.k}{clusternum}(:,2);
    end
    if any(pixelX==imgInfoIn.dims(1))
        pixFlags(clusternum)=1;
    elseif any(pixelY==imgInfoIn.dims(2))
        pixFlags(clusternum)=1;
    end
end
clusterInfoOut.pixflags=pixFlags;


function edgepointlist=getEdgePoints(points,orderflag)   

edgepointlist=[];
divisionpoints=[];

vectorDims=size(points);
if vectorDims(1)<vectorDims(2)      % separate X and Y positions by column
    points=points';
end


for indx=1:length(points)            % get edge points (not ordered)
    dists=sqrt((points(:,1)-points(indx,1)).^2+(points(:,2)-points(indx,2)).^2);
    adjacentAll=find(dists<=sqrt(2));
    adjacentByOne=find(abs(dists-1)<=0.00001);
    if length(adjacentAll)<9 && ~isempty(adjacentByOne)
        edgepointlist=vertcat(edgepointlist,[points(indx,:) length(adjacentAll) length(adjacentByOne)]);
    end
end

if orderflag
    for indx=1:length(edgepointlist)       % search for meeting points of objects that will cause ordering algorithm to fail
        dists=sqrt((edgepointlist(:,1)-edgepointlist(indx,1)).^2+(edgepointlist(:,2)-edgepointlist(indx,2)).^2);
        adjacent=find(abs(dists-1)<0.00001);
        if length(adjacent)>2 && edgepointlist(indx,3)<8
            divisionpoints=vertcat(divisionpoints,edgepointlist(indx,:));
        end        
    end 

    % order list by adjacency so can use polyarea
    orderedEdgePointList=single(nan(length(edgepointlist)+1,4));
    divisionpointscounter=1;
    
    if isempty(divisionpoints)
        currentpoint=edgepointlist(1,:);
        backtrack=edgepointlist(2,:);
    else       
        currentpoint=divisionpoints(divisionpointscounter,:);     % start at a division point if there is one so completes that object before moving to next
        backtrack=currentpoint;
    end    
    orderedEdgePointList(1,:)=currentpoint;
    
    
    for ind=2:length(edgepointlist)
        dists=sqrt((edgepointlist(:,1)-currentpoint(1,1)).^2+(edgepointlist(:,2)-currentpoint(1,2)).^2);
        adjacentInds=find(dists<=sqrt(2));
        if length(adjacentInds)>1  % there is an adjacent point, up to sqrt2 away 
            adjacentPoints=edgepointlist(adjacentInds,:);
            unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
            if ~isempty(unusedAdjacentPoints)  % there is an adjacent point, up to sqrt2 away, that has not been used
                adjacentInds=find(abs(dists-1)<0.00001);
                if ~isempty(adjacentInds)  % there is an adjacent point, up to 1 away
                    adjacentPoints=edgepointlist(adjacentInds,:);
                    unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                    if ~isempty(unusedAdjacentPoints)  % there is an adjacent point, up to 1 away, that has not been used
                        [~,sortind]=sort(unusedAdjacentPoints(:,3));
                        unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                        currentpoint=unusedAdjacentPoints(1,:);
                        if size(unusedAdjacentPoints,1)>1
                            backtrack=unusedAdjacentPoints(2:end,:);
                        end
                    else % the unused adjacent point is sqrt2 away
                        adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                        adjacentPoints=edgepointlist(adjacentInds,:);
                        unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                        [~,sortind]=sort(unusedAdjacentPoints(:,3));
                        unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                        currentpoint=unusedAdjacentPoints(1,:);
                        if size(unusedAdjacentPoints,1)>1
                            backtrack=unusedAdjacentPoints(2:end,:);
                        end
                    end
                else % the unused adjacent point is sqrt2 away
                    adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                    adjacentPoints=edgepointlist(adjacentInds,:);
                    unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                    [~,sortind]=sort(unusedAdjacentPoints(:,3));
                    unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                    currentpoint=unusedAdjacentPoints(1,:);
                    if size(unusedAdjacentPoints,1)>1
                        backtrack=unusedAdjacentPoints(2:end,:);
                    end
                end                
            elseif  ~isempty(setdiff(backtrack,orderedEdgePointList,'rows')) % any adjacent points were used already, but there are still unused points in the backtrack variable
                unusedAdjacentPoints=setdiff(backtrack,orderedEdgePointList,'rows');
                [~,sortind]=sort(unusedAdjacentPoints(:,3));
                unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                currentpoint=unusedAdjacentPoints(1,:);
                if size(unusedAdjacentPoints,1)>1
                    backtrack=unusedAdjacentPoints(2:end,:);
                end
                        
            elseif size(divisionpoints,1)>=divisionpointscounter    % all adjacent points and backtrack points have been used, but there is a division point to try
                currentpoint=divisionpoints(divisionpointscounter,:);
                dists=sqrt((edgepointlist(:,1)-currentpoint(1,1)).^2+(edgepointlist(:,2)-currentpoint(1,2)).^2);
                adjacentInds=find(dists<=sqrt(2)&dists>0);  %%%%% changed 12/15/14
                adjacentPoints=edgepointlist(adjacentInds,:);
                unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                if ~isempty(unusedAdjacentPoints)  % this division point will work because there is an adjacent point, up to sqrt2 away, that has not been used
                    adjacentInds=find(abs(dists-1)<0.00001);
                    if ~isempty(adjacentInds)  % there is an adjacent point, up to 1 away
                        adjacentPoints=edgepointlist(adjacentInds,:);
                        unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                        if ~isempty(unusedAdjacentPoints)  % there is an adjacent point, up to 1 away, that has not been used
                            [~,sortind]=sort(unusedAdjacentPoints(:,3));
                            unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                            currentpoint=unusedAdjacentPoints(1,:);
                            if size(unusedAdjacentPoints,1)>1
                                backtrack=unusedAdjacentPoints(2:end,:);
                            end
                        else % the unused adjacent point is sqrt2 away
                            adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                            adjacentPoints=edgepointlist(adjacentInds,:);
                            unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                            [~,sortind]=sort(unusedAdjacentPoints(:,3));
                            unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                            currentpoint=unusedAdjacentPoints(1,:);
                            if size(unusedAdjacentPoints,1)>1
                                backtrack=unusedAdjacentPoints(2:end,:);
                            end
                        end
                    else % the unused adjacent point is sqrt2 away
                        adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                        adjacentPoints=edgepointlist(adjacentInds,:);
                        unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                        [~,sortind]=sort(unusedAdjacentPoints(:,3));
                        unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                        currentpoint=unusedAdjacentPoints(1,:);
                        if size(unusedAdjacentPoints,1)>1
                            backtrack=unusedAdjacentPoints(2:end,:);
                        end
                    end
                elseif size(divisionpoints,1)>=divisionpointscounter+1 % if there is another division point to try
                    divisionpointscounter=divisionpointscounter+1;
                    currentpoint=divisionpoints(divisionpointscounter,:);
                    dists=sqrt((edgepointlist(:,1)-currentpoint(1,1)).^2+(edgepointlist(:,2)-currentpoint(1,2)).^2);
                    adjacentInds=find(dists<=sqrt(2)&dists>0);
                    adjacentPoints=edgepointlist(adjacentInds,:);
                    unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                    if ~isempty(unusedAdjacentPoints)
                        adjacentInds=find(abs(dists-1)<0.00001);
                        if ~isempty(adjacentInds)  % there is an adjacent point, up to 1 away
                            adjacentPoints=edgepointlist(adjacentInds,:);
                            unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                            if ~isempty(unusedAdjacentPoints)  % there is an adjacent point, up to 1 away, that has not been used
                                [~,sortind]=sort(unusedAdjacentPoints(:,3));
                                unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                                currentpoint=unusedAdjacentPoints(1,:);
                                if size(unusedAdjacentPoints,1)>1
                                    backtrack=unusedAdjacentPoints(2:end,:);
                                end
                            else % the unused adjacent point is sqrt2 away
                                adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                                adjacentPoints=edgepointlist(adjacentInds,:);
                                unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                                [~,sortind]=sort(unusedAdjacentPoints(:,3));
                                unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                                currentpoint=unusedAdjacentPoints(1,:);
                                if size(unusedAdjacentPoints,1)>1
                                    backtrack=unusedAdjacentPoints(2:end,:);
                                end
                            end
                        else % the unused adjacent point is sqrt2 away
                            adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                            adjacentPoints=edgepointlist(adjacentInds,:);
                            unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                            [~,sortind]=sort(unusedAdjacentPoints(:,3));
                            unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                            currentpoint=unusedAdjacentPoints(1,:);
                            if size(unusedAdjacentPoints,1)>1
                                backtrack=unusedAdjacentPoints(2:end,:);
                            end
                        end
                    else
                        unusedPoints=setdiff(edgepointlist,orderedEdgePointList,'rows');
                        currentpoint=unusedPoints(1,:);
                        if size(unusedPoints,1)>1
                            backtrack=unusedPoints(2:end,:);
                        end
                    end
                else % the unused adjacent point is sqrt2 away
                    unusedPoints=setdiff(edgepointlist,orderedEdgePointList,'rows');
                    currentpoint=unusedPoints(1,:);
                    if size(unusedPoints,1)>1
                        backtrack=unusedPoints(2:end,:);
                    end
                end                    
            end
                        
        elseif  ~isempty(setdiff(backtrack,orderedEdgePointList,'rows')) % any adjacent points were used already, but there are still unused points in the backtrack variable
            unusedAdjacentPoints=setdiff(backtrack,orderedEdgePointList,'rows');
            [~,sortind]=sort(unusedAdjacentPoints(:,3));
            unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
            currentpoint=unusedAdjacentPoints(1,:);
            if size(unusedAdjacentPoints,1)>1
                backtrack=unusedAdjacentPoints(2:end,:);
            end

        elseif size(divisionpoints,1)>=divisionpointscounter    % all adjacent points and backtrack points have been used, but there is a division point to try
            currentpoint=divisionpoints(divisionpointscounter,:);
            dists=sqrt((edgepointlist(:,1)-currentpoint(1,1)).^2+(edgepointlist(:,2)-currentpoint(1,2)).^2);
            adjacentInds=find(dists<=sqrt(2)&dists>0);
            adjacentPoints=edgepointlist(adjacentInds,:);
            unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
            if ~isempty(unusedAdjacentPoints)  % this division point will work because there is an adjacent point, up to sqrt2 away, that has not been used
                adjacentInds=find(abs(dists-1)<0.00001);
                if ~isempty(adjacentInds)  % there is an adjacent point, up to 1 away
                    adjacentPoints=edgepointlist(adjacentInds,:);
                    unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                    if ~isempty(unusedAdjacentPoints)  % there is an adjacent point, up to 1 away, that has not been used
                        [~,sortind]=sort(unusedAdjacentPoints(:,3));
                        unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                        currentpoint=unusedAdjacentPoints(1,:);
                        if size(unusedAdjacentPoints,1)>1
                            backtrack=unusedAdjacentPoints(2:end,:);
                        end
                    else % the unused adjacent point is sqrt2 away
                        adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                        adjacentPoints=edgepointlist(adjacentInds,:);
                        unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                        [~,sortind]=sort(unusedAdjacentPoints(:,3));
                        unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                        currentpoint=unusedAdjacentPoints(1,:);
                        if size(unusedAdjacentPoints,1)>1
                            backtrack=unusedAdjacentPoints(2:end,:);
                        end
                    end
                else % the unused adjacent point is sqrt2 away
                    adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                    adjacentPoints=edgepointlist(adjacentInds,:);
                    unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                    [~,sortind]=sort(unusedAdjacentPoints(:,3));
                    unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                    currentpoint=unusedAdjacentPoints(1,:);
                    if size(unusedAdjacentPoints,1)>1
                        backtrack=unusedAdjacentPoints(2:end,:);
                    end
                end
            elseif size(divisionpoints,1)>=divisionpointscounter+1 % if there is another division point to try
                divisionpointscounter=divisionpointscounter+1;
                currentpoint=divisionpoints(divisionpointscounter,:);
                dists=sqrt((edgepointlist(:,1)-currentpoint(1,1)).^2+(edgepointlist(:,2)-currentpoint(1,2)).^2);
                adjacentInds=find(dists<=sqrt(2)&dists>0);
                adjacentPoints=edgepointlist(adjacentInds,:);
                unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                if ~isempty(unusedAdjacentPoints)
                    adjacentInds=find(abs(dists-1)<0.00001);
                    if ~isempty(adjacentInds)  % there is an adjacent point, up to 1 away
                        adjacentPoints=edgepointlist(adjacentInds,:);
                        unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                        if ~isempty(unusedAdjacentPoints)  % there is an adjacent point, up to 1 away, that has not been used
                            [~,sortind]=sort(unusedAdjacentPoints(:,3));
                            unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                            currentpoint=unusedAdjacentPoints(1,:);
                            if size(unusedAdjacentPoints,1)>1
                                backtrack=unusedAdjacentPoints(2:end,:);
                            end
                        else % the unused adjacent point is sqrt2 away
                            adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                            adjacentPoints=edgepointlist(adjacentInds,:);
                            unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                            [~,sortind]=sort(unusedAdjacentPoints(:,3));
                            unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                            currentpoint=unusedAdjacentPoints(1,:);
                            if size(unusedAdjacentPoints,1)>1
                                backtrack=unusedAdjacentPoints(2:end,:);
                            end
                        end
                    else % the unused adjacent point is sqrt2 away
                        adjacentInds=find(abs(dists-sqrt(2))<0.00001);
                        adjacentPoints=edgepointlist(adjacentInds,:);
                        unusedAdjacentPoints=setdiff(adjacentPoints,orderedEdgePointList,'rows');
                        [~,sortind]=sort(unusedAdjacentPoints(:,3));
                        unusedAdjacentPoints=unusedAdjacentPoints(sortind,:);
                        currentpoint=unusedAdjacentPoints(1,:);
                        if size(unusedAdjacentPoints,1)>1
                            backtrack=unusedAdjacentPoints(2:end,:);
                        end
                    end
                else
                    unusedPoints=setdiff(edgepointlist,orderedEdgePointList,'rows');
                    currentpoint=unusedPoints(1,:);
                    if size(unusedPoints,1)>1
                        backtrack=unusedPoints(2:end,:);
                    end
                end
            else
                unusedPoints=setdiff(edgepointlist,orderedEdgePointList,'rows');
                currentpoint=unusedPoints(1,:);
                if size(unusedPoints,1)>1
                    backtrack=unusedPoints(2:end,:);
                end
            end

        end
        orderedEdgePointList(ind,:)=currentpoint;
    end

    orderedEdgePointList(end,:)=orderedEdgePointList(1,:);
    edgepointlist=double(orderedEdgePointList(:,1:2));

    if ~isempty(divisionpoints)
        divisionpoints=divisionpoints(:,1:2);
    end
end

function elbowPoints=getElbowPoints(points,edgePoints)

vectorDims=size(points);
if vectorDims(1)<vectorDims(2)      % separate X and Y positions by column
    points=points';
end

if nargin<2
    edgePoints=getEdgePoints(points,1);
end
hullinds=convhull(edgePoints(:,1),edgePoints(:,2),'simplify',true);     % making from edgePoints instead of points so that indices match with CEPinds
hullPoints=edgePoints(hullinds,:);
[concaveEdgePoints,CEPinds]=setdiff(edgePoints,hullPoints,'rows');

minClosestHullInds=nan(length(concaveEdgePoints),1);  

distancesToHull=nan(length(concaveEdgePoints),1);  

for CEP=1:length(concaveEdgePoints)

    [~,closestHullInd]=min(abs(hullinds-CEPinds(CEP)));
    if hullinds(closestHullInd)<CEPinds(CEP)
        hullInd1=hullinds(closestHullInd);
        hullInd2=hullinds(closestHullInd+1);
    elseif closestHullInd==1
        hullInd1=hullinds(length(hullinds)-1);   % last hull point/index is repeat of the first
        hullInd2=hullinds(length(hullinds));
    else
        hullInd1=hullinds(closestHullInd-1);
        hullInd2=hullinds(closestHullInd);
    end
    
    distancesToHull(CEP)=abs(det([(edgePoints(hullInd2,:)-edgePoints(hullInd1,:))' (edgePoints(hullInd1,:)-edgePoints(CEPinds(CEP),:))'])/sqrt((edgePoints(hullInd2,1)-edgePoints(hullInd1,1))^2+(edgePoints(hullInd2,2)-edgePoints(hullInd1,2))^2));        
    minClosestHullInds(CEP)=min([hullInd1 hullInd2]);
    
end

indices=(1:length(concaveEdgePoints))';
concaveEdgePoints=[concaveEdgePoints distancesToHull minClosestHullInds indices];
elbowPoints=[];

cutoff=median(distancesToHull)+1;
% cutoff=mean(distancesToHull)+std(distancesToHull);
if cutoff<2
    cutoff=2;
elseif cutoff>5
    cutoff=5;
end


equaldepthpoints=[];

for point=1:size(concaveEdgePoints,1)

    adjacentCEPs=concaveEdgePoints(concaveEdgePoints(:,4)==concaveEdgePoints(point,4),:);

    adjacentCEPs=setdiff(adjacentCEPs,concaveEdgePoints(point,:),'rows');
    if ~isempty(adjacentCEPs)
        deepestAdjacentCEPdepth=max(adjacentCEPs(:,3));
        if concaveEdgePoints(point,3)>cutoff
            if deepestAdjacentCEPdepth<concaveEdgePoints(point,3)  
                elbowPoints=vertcat(elbowPoints,concaveEdgePoints(point,:));
            elseif all(abs(adjacentCEPs(:,3)-concaveEdgePoints(point,3))<0.000001) && concaveEdgePoints(point,3)>cutoff
                elbowPoints=vertcat(elbowPoints,concaveEdgePoints(point,:));
            elseif deepestAdjacentCEPdepth==concaveEdgePoints(point,3)
                adjacentCEPs=concaveEdgePoints(concaveEdgePoints(:,4)==concaveEdgePoints(point,4),:);
                tempequaldepthpoints=adjacentCEPs(adjacentCEPs(:,3)==deepestAdjacentCEPdepth,:);
                [~,tempmin]=min(tempequaldepthpoints(:,5));
                [~,tempmax]=max(tempequaldepthpoints(:,5));
                equaldepthpoints=vertcat(equaldepthpoints,tempequaldepthpoints(tempmin,:),tempequaldepthpoints(tempmax,:));
            end
        end
    end          
end

if ~isempty(equaldepthpoints)
    equaldepthpoints=unique(equaldepthpoints,'rows');
    elbowPoints=vertcat(elbowPoints,equaldepthpoints);
end
if ~isempty(elbowPoints)
    elbowPoints=elbowPoints(:,1:4);
end

function equivalentScore=overlappingCircleEquivalent(k,pixels,centroids,objAreas)

errorsOutsideCircle=[];nonObjPixInCircle=[];

for objectNum=1:k
    % make edges of equivalent circle
    radius=sqrt(objAreas(objectNum)/pi);
    theta=linspace(0,2*pi,1000);
    circleXs=centroids(objectNum,1)+radius*cos(theta);
    circleYs=centroids(objectNum,2)+radius*sin(theta);
    % get object pixels outside circle
    objIndsInCircle=inpolygon(pixels(:,1),pixels(:,2),circleXs,circleYs);
    objPixInCircle=pixels(objIndsInCircle,:);
    objPixOutsideCircle=setdiff(pixels,objPixInCircle,'rows');
    % get pixels inside circle that aren't in object
    [x,y]=meshgrid(floor(min(circleXs))-1:ceil(max(circleXs))+1,floor(min(circleYs)):ceil(max(circleYs)));
    pixarray=[reshape(x,size(x,2)*size(x,1),1) reshape(y,size(y,2)*size(y,1),1)];
    pixInThisCircle=pixarray(inpolygon(pixarray(:,1),pixarray(:,2),circleXs,circleYs),:);
    nonObjPixInThisCircle=setdiff(pixInThisCircle,pixels,'rows');
    nonObjPixInThisCircleDist=abs(sqrt((nonObjPixInThisCircle(:,1)-centroids(objectNum,1)).^2+(nonObjPixInThisCircle(:,2)-centroids(objectNum,2)).^2)-radius);
    nonObjPixInCircle=vertcat(nonObjPixInCircle,horzcat(nonObjPixInThisCircle,nonObjPixInThisCircleDist));
    % for error pixels (outside) get distance to circle (error)
    if ~isempty(errorsOutsideCircle)
        % for objPixOutsideCircle already part of errorsOutsideCircle, increment count of num circles pixel is outside of
        [~,ind1,ind2]=intersect(errorsOutsideCircle(:,1:2),objPixOutsideCircle,'rows');
        errorsOutsideCircle(ind1,4)=errorsOutsideCircle(ind1,4)+1;
        errorsOutsideCircle(ind1,3)=min(horzcat(errorsOutsideCircle(ind1,3),abs(sqrt((objPixOutsideCircle(ind2,1)-centroids(objectNum,1)).^2+(objPixOutsideCircle(ind2,2)-centroids(objectNum,2)).^2)-radius)),[],2);
        % add new objPixOutsideCircle to errorsOutsideCircle and indicate they've been added once
        objPixOutsideCircle=setdiff(objPixOutsideCircle,errorsOutsideCircle(:,1:2),'rows');
        errorsOutsideCircle=vertcat(errorsOutsideCircle,horzcat(objPixOutsideCircle,abs(sqrt((objPixOutsideCircle(:,1)-centroids(objectNum,1)).^2+(objPixOutsideCircle(:,2)-centroids(objectNum,2)).^2)-radius),ones(size(objPixOutsideCircle,1),1)));        
    else 
        errorsOutsideCircle=vertcat(errorsOutsideCircle,horzcat(objPixOutsideCircle,abs(sqrt((objPixOutsideCircle(:,1)-centroids(objectNum,1)).^2+(objPixOutsideCircle(:,2)-centroids(objectNum,2)).^2)-radius),ones(size(objPixOutsideCircle,1),1)));
    end
end

% eliminate repeated non-object pixels in circles
[~,sortInds]=sort(nonObjPixInCircle(:,3),1,'ascend');  
nonObjPixInCircle=nonObjPixInCircle(sortInds,:);  % inside error pixels sorted by distance to circle
[~,uniqueInds]=unique(nonObjPixInCircle(:,1:2),'rows','first');
nonObjPixInCircle=nonObjPixInCircle(uniqueInds,:);  % eliminate repeated pixels with larger distances
% get set of pixels not inside any circle
errorsOutsideCircle=errorsOutsideCircle(errorsOutsideCircle(:,4)==k,1:3);

errorDists=vertcat(nonObjPixInCircle(:,3),errorsOutsideCircle(:,3));
% sum of squared distances => equivalentScore
equivalentScore=sum(errorDists.^2);

function [matchScore,nearestNonseg] = getCentroidMatchScore(centroidsToMatch,centroidsForComparisonIn,z,nearestNonseg)


nonSeg=size(centroidsToMatch,1)<2;
nonemptyZ=find(~cellfun(@isempty,centroidsForComparisonIn));
global zstepsize
zAbove=z-2/zstepsize;
zBelow=z+2/zstepsize;
compareAbove=0;compareBelow=0;
if zAbove>0
    if ismember(zAbove,nonemptyZ)
        compareAbove=1;
    end
end
if zBelow>0
    if ismember(zBelow,nonemptyZ)
        compareBelow=1;
    end
end

compareAB=horzcat(compareAbove,compareBelow);

if any(compareAB)
    if nargin<4
        nearestNonseg=nan(2,2);
    end
    centroidsToMatch=horzcat(centroidsToMatch,nan(size(centroidsToMatch,1),2));
    
    for centroidNum=1:size(centroidsToMatch,1)
        if compareAbove        
            centroidsAbove=centroidsForComparisonIn{zAbove};
            if ~nonSeg
                centroidsAbove=setdiff(centroidsAbove,nearestNonseg(1,:),'rows');        
            end
            XYdists=sqrt((centroidsAbove(:,1)-centroidsToMatch(centroidNum,1)).^2+(centroidsAbove(:,2)-centroidsToMatch(centroidNum,2)).^2);
            if ~isempty(XYdists)
                [centroidsToMatch(centroidNum,3),minind]=min(XYdists);
            end
            if nonSeg
                nearestNonseg(1,:)=centroidsAbove(minind,1:2);
            end
        end
        if compareBelow
            centroidsBelow=centroidsForComparisonIn{zBelow};
            if ~nonSeg
                centroidsBelow=setdiff(centroidsBelow,nearestNonseg(2,:),'rows');        
            end
            XYdists=sqrt((centroidsBelow(:,1)-centroidsToMatch(centroidNum,1)).^2+(centroidsBelow(:,2)-centroidsToMatch(centroidNum,2)).^2);
            if ~isempty(XYdists)
                [centroidsToMatch(centroidNum,4),minind]=min(XYdists);
            end
            if nonSeg
                nearestNonseg(2,:)=centroidsBelow(minind,1:2);
            end
        end        
    end
    matchScore=mean(min(centroidsToMatch(:,3:4),[],2));    
else
    matchScore=nan(1);
end

function imgOUT=makeImFromPix(pixels,imgIN)

imgOUT=imgIN*0;
imgOUT(sub2ind(size(imgOUT),round(pixels(:,2)),round(pixels(:,1))))=1;



