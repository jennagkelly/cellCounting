function [colocMatOUT,channelImgSeries,unprocessedImgSeries] = evalMultChannels(centroidsIN,fileinputnames,zsteps,channelTypes,fillCutoff)


% Channel Type 1: soma filled by stain, even staining (NeuN)
% Channel Type 2: membrane ring, signal sometimes spotty (Kv3.1b)
% Channel Type 3: PV
% Channel Type 4: GABA

% centroidsIN is after evaluating across z planes: columns are x pos, y pos, zpos ***

% Coloc Type 0: unlabeled (DAPI only)
% Coloc Type 1: Channel 1 only
% Coloc Type 2: Channel 2 only
% Coloc Type 3: Channel 3 only
% Coloc Type 4: Channels 1&2 
% Coloc Type 5: Channels 1&3
% Coloc Type 6: Channels 2&3
% Coloc Type 7: Channels 1,2,&3 (all)

numChannels=length(fileinputnames);
colocMat=horzcat(centroidsIN,nan(size(centroidsIN,1),numChannels+1));   % nan init for each channel pos/neg, and summary
channelImgSeries=struct([]);
unprocessedImgSeries=struct([]);

for channel=1:numChannels        
%     filename=horzcat(fileinputname,'ch0',num2str(channel),'.tif');
    filename=horzcat(fileinputnames{channel},'.tif');
    
    for z = zsteps
        img=imread(filename,z);
        unprocessedImgSeries{channel}{z}=img;
        channelType=channelTypes(channel);
        switch channelType
            case 1
                imgAbove=imread(filename,z-1);
                imgBelow=imread(filename,z+1);
                img=processChannelType1(img,imgAbove,imgBelow);

            case 2
                imgAbove=imread(filename,z-1);
                imgAbove2=imread(filename,z-2);
                imgBelow=imread(filename,z+1);
                imgBelow2=imread(filename,z+2);
                thisZcentroids=centroidsIN(round(centroidsIN(:,4))==z,:);
                img=processChannelType2(img,imgAbove,imgBelow,imgAbove2,imgBelow2,40);
                
            case 3
                imgAbove=imread(filename,z-1);
                imgBelow=imread(filename,z+1);
                imgAbove2=imread(filename,z-2);
                imgBelow2=imread(filename,z+2);
                img=processChannelType3(img,imgAbove,imgBelow,imgAbove2,imgBelow2);
                
            case 4
                imgAbove=imread(filename,z-1);
                imgBelow=imread(filename,z+1);
                imgAbove2=imread(filename,z-2);
                imgBelow2=imread(filename,z+2);
                img=processChannelType4(img,imgAbove,imgBelow,imgAbove2,imgBelow2);
        end
        channelImgSeries{channel}{z}=img;
    end
    for centroidNum=1:size(centroidsIN,1)
        centroidZ=round(centroidsIN(centroidNum,4));
        imgProps=regionprops(channelImgSeries{channel}{centroidZ},'PixelList');
        pixels=cat(1,imgProps.PixelList);
        evalRadius = sqrt(centroidsIN(centroidNum,7)/pi)*2/3;     
        if ~isempty(pixels)
            propFilled=evalCirclePix(pixels,centroidsIN,evalRadius,centroidNum);  
        else
            propFilled=0;
        end
        if propFilled>=fillCutoff
            colocMat(centroidNum,channel+size(centroidsIN,2))=1;
        else
            colocMat(centroidNum,channel+size(centroidsIN,2))=0;
        end    
    end
end

for row=1:size(colocMat,1)
    singleChannelAssignments=colocMat(row,size(centroidsIN,2)+1:size(centroidsIN,2)+numChannels);
    if sum(singleChannelAssignments)<1
        colocMat(row,size(centroidsIN,2)+numChannels+1)=0;
    elseif sum(singleChannelAssignments)<2
        colocMat(row,size(centroidsIN,2)+numChannels+1)=find(singleChannelAssignments);
    elseif sum(singleChannelAssignments)<3
        posChans=sum(find(singleChannelAssignments)); % 1+2=3 (type4), 1+3=4 (type5), 2+3=5(type6)
        colocMat(row,size(centroidsIN,2)+numChannels+1)=posChans+1;
    elseif sum(singleChannelAssignments)==3   % change to an else after debugging (during debugging, some errors will cause this to stay NaN)
        colocMat(row,size(centroidsIN,2)+numChannels+1)=7;
    end
end
colocMatOUT=colocMat;

end


function imgOUT=processChannelType1(imgIN,imgAbove,imgBelow)  % solid: NeuN
    imgAvg1=(double(imgIN)+double(imgAbove)+double(imgBelow))/3;
    padsize = 50;
    imgAvg2=padarray(imgAvg1,[padsize padsize],'symmetric');
    imgOUT=double(bradley(imgAvg2,[200 200],20,'symmetric')); 
        
    se=strel('disk',3);se2=strel('disk',2);imgOUT=imopen(imdilate(imdilate(imdilate(imfill(imerode(imerode(imerode(imgOUT,se),se),se),'holes'),se),se),se),se);
    imgOUT = imgOUT(padsize+1:size(imgIN,1)+padsize,padsize+1:size(imgIN,2)+padsize);
end

function imgOUT=processChannelType2(imgIN,imgAbove,imgBelow,imgAbove2,imgBelow2,distTolerance)   % rings: Kv3
    
    se=strel('ball',3,3);se1=strel('ball',200,200);
    imgIN=medfilt2(double(imgIN)-imopen(double(imgIN),se1),[5 5],'symmetric');
    imgAbove=medfilt2(double(imgAbove)-imopen(double(imgAbove),se1),[5 5],'symmetric');
    imgAbove2=medfilt2(double(imgAbove2)-imopen(double(imgAbove2),se1),[5 5],'symmetric');
    imgBelow2=medfilt2(double(imgBelow2)-imopen(double(imgBelow2),se1),[5 5],'symmetric');
    imgBelow=medfilt2(double(imgBelow)-imopen(double(imgBelow),se1),[5 5],'symmetric');
    imgAvg1=(imgIN+imgAbove+imgBelow+imgAbove2+imgBelow2)/5;
    
    padsize = 50;
    imgAvg1b=padarray(imgAvg1,[padsize padsize],'symmetric');
    imgAvg2=imgAvg1b-min(min(imgAvg1b));
    imgAvg2=imgAvg2/max(max(imgAvg2));

    [fx,fy]=gradient(imgAvg2);
    gradIm=abs(fx+fy);
    [circles,~]=imfindcircles(gradIm,[round(4*1024/243) round(8.9*1024/243)],'sensitivity',0.97);
    mask=double(imgAvg2*0);
    if ~isempty(circles)
        mask(sub2ind(size(mask),round(circles(:,2)),round(circles(:,1))))=1;
        kernel=fspecial('gaussian',[20*4 20*4],20);
        maskc=conv2(mask,kernel,'same');
        maskc=(maskc/max(max(maskc))+1)/2;  
        imgAvg2=imgAvg2.*maskc;
    else
        disp('no circles')
        imgAvg2=imgAvg2;
    end




    imgAvg2=medfilt2(imgAvg2); 
    imgAvg3=bwareaopen(imfill(bradley(imgAvg2,[25 25],-5,'symmetric').*bradley(imgAvg2,[50 50],-15,'symmetric').*bradley(imgAvg2,[400 400],-20,'symmetric').*bradley(imgAvg2,[800 800],-15,'symmetric'),'holes'),200);
    
    se1=strel('disk',2);se=strel('disk',3);mask=imerode(imerode(imfill(bwareaopen(imdilate(imdilate(imgAvg3,se1),se1),1000),'holes'),se1),se);

    se=strel('disk',4);
    imgE=imerode(mask,se);  
    se=strel('disk',3);imgOUT=bwareaopen(imclose(imgE,se),300);
    imgOUT = imgOUT(padsize+1:size(imgIN,1)+padsize,padsize+1:size(imgIN,2)+padsize);
end

function imgOUT=processChannelType3(imgIN,imgAbove,imgBelow,imgAbove2,imgBelow2)  % solid: PV
    imgAvg1=(double(imgIN)+double(imgAbove)+double(imgBelow)+double(imgAbove2)+double(imgBelow2))/5;
    imgAvg2=medfilt2(imgAvg1,[10 10]);
    se=strel('disk',2);
    imgOUT=bwareaopen(imopen(imfill(bwareaopen(double(bradley(imgAvg2,[200 200],-100,'symmetric').*bradley(imgAvg2,[800 800],-100,'symmetric')),50),'holes'),se),300); 
end

function imgOUT=processChannelType4(imgIN,imgAbove,imgBelow,imgAbove2,imgBelow2)  % solid: GABA
    padsize = 50;
    se=strel('ball',200,200);
    imgIN=medfilt2(double(imgIN)-imopen(double(imgIN),se),[10 10],'symmetric');
    imgAbove=medfilt2(double(imgAbove)-imopen(double(imgAbove),se),[10 10],'symmetric');
    imgAbove2=medfilt2(double(imgAbove2)-imopen(double(imgAbove2),se),[10 10],'symmetric');
    imgBelow2=medfilt2(double(imgBelow2)-imopen(double(imgBelow2),se),[10 10],'symmetric');
    imgBelow=medfilt2(double(imgBelow)-imopen(double(imgBelow),se),[10 10],'symmetric');
    imgAvg1a=(imgIN+imgAbove+imgBelow+imgAbove2+imgBelow2)/5;
    imgAvg1a=padarray(imgAvg1a,[padsize padsize],'symmetric');
    imgAvg1=imgAvg1a-min(min(imgAvg1a));
    imgAvg1=imgAvg1/max(max(imgAvg1));    

    [fx,fy]=gradient(imgAvg1);
    gradIm=abs(fx+fy);
    [circles1,~]=imfindcircles(gradIm,[round(4*1024/243) round(7.99*1024/243)],'sensitivity',0.96);
    [circles2,~]=imfindcircles(gradIm,[round(8*1024/243) round(13*1024/243)],'sensitivity',0.96);
    circles=vertcat(circles1,circles2);
    mask=double(imgAvg1*0);
    if ~isempty(circles)
        mask(sub2ind(size(mask),round(circles(:,2)),round(circles(:,1))))=1;
        kernel=fspecial('gaussian',[20*4 20*4],20);
        maskc=conv2(mask,kernel,'same');
        maskc=(maskc/max(max(maskc))+1)/2;  
        imgAvg1=imgAvg1.*maskc;
    else
        disp('no circles')
        imgAvg1=imgAvg1;
    end
    imgAvg1=double(bradley(imgAvg1,[50 50],-15,'symmetric').*bradley(imgAvg1,[400 400],-30,'symmetric').*bradley(imgAvg1,[800 800],-20,'symmetric')); 
    imgAvg1=imfill(imgAvg1,'holes');
    imgAvg1=bwareaopen(imgAvg1,300);
    se=strel('disk',2);imgOUT=imdilate(imdilate(imdilate(imerode(imerode(imerode(imgAvg1,se),se),se),se),se),se);%imgOUT=imopen(imdilate(imdilate(imdilate(imfill(imerode(imerode(imerode(imgAvg1,se),se),se),'holes'),se),se),se),se);
    imgAvg1a=imcomplement(imgAvg1a);
    imgAvg1a=imgAvg1a-min(min(imgAvg1a));
    imgAvg1a=imgAvg1a/max(max(imgAvg1a));
    imgAvg1a=imdilate(bradley(imgAvg1a,[400 400],5,'symmetric').*bradley(imgAvg1a,[800 800],5,'symmetric'),se);    
    imgOUT=imgOUT-imgAvg1a;
    imgOUT(imgOUT<1)=0;
    imgOUT = imgOUT(padsize+1:size(imgIN,1)+padsize,padsize+1:size(imgIN,2)+padsize);
end

function propFilled=evalCirclePix(pixels,centroids,radius,centroidNum)
% make edges of equivalent circle
theta=linspace(0,2*pi,1000);
circleXs=centroids(centroidNum,2)+radius*cos(theta);
circleYs=centroids(centroidNum,3)+radius*sin(theta);

% get pixels inside circle that aren't in image
[x,y]=meshgrid(floor(min(circleXs))-1:ceil(max(circleXs))+1,floor(min(circleYs)):ceil(max(circleYs)));
pixarray=[reshape(x,size(x,2)*size(x,1),1) reshape(y,size(y,2)*size(y,1),1)];
pixInThisCircle=pixarray(inpolygon(pixarray(:,1),pixarray(:,2),circleXs,circleYs),:);
nonObjPixInCircle=setdiff(pixInThisCircle,pixels,'rows');

% % eliminate repeated non-object pixels in circles
% [~,sortInds]=sort(nonObjPixInCircle(:,3),1,'ascend');  
% nonObjPixInCircle=nonObjPixInCircle(sortInds,:);  % inside error pixels sorted by distance to circle
% [~,uniqueInds]=unique(nonObjPixInCircle(:,1:2),'rows','first');
% nonObjPixInCircle=nonObjPixInCircle(uniqueInds,:);  % eliminate repeated pixels with larger distances

% get proportion of circle filled
propFilled=(size(pixInThisCircle,1)-size(nonObjPixInCircle,1))/size(pixInThisCircle,1);
end