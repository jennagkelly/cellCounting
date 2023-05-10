function plot3dCentroids(centroidsGroupedIN,imagesIN,saveflag)

colcount=size(centroidsGroupedIN.threeD,2);
centroidPlotList=horzcat(centroidsGroupedIN.threeD,nan(size(centroidsGroupedIN.threeD,1),3));
% colcount=size(centroidsGroupedIN.stereo,2);
% centroidPlotList=horzcat(centroidsGroupedIN.stereo,nan(size(centroidsGroupedIN.stereo,1),3));
plotZ=centroidsGroupedIN.nonemptyZ;


for ind=1:size(centroidPlotList,1)
    centroidPlotList(ind,colcount+1)=round(centroidPlotList(ind,4));
    centroidPlotList(ind,colcount+2)=centroidPlotList(ind,colcount+1)-1;
    centroidPlotList(ind,colcount+3)=centroidPlotList(ind,colcount+1)+1;
    if centroidPlotList(ind,colcount+1)<plotZ(1)
        centroidPlotList(ind,colcount+1)=plotZ(1);
        centroidPlotList(ind,colcount+2)=plotZ(1);
    elseif centroidPlotList(ind,colcount+2)<plotZ(1)
        centroidPlotList(ind,colcount+2)=plotZ(1);
    end
end

cd 'C:\jgk\cellCounting\plottemp'

for z=plotZ
    figure('visible','off');imshow(imagesIN{z});
    hold on;
%     plotGreen=vertcat(centroidPlotList(centroidPlotList(:,colcount+2)==z,:),centroidPlotList(centroidPlotList(:,colcount+3)==z,:));
    plotRed=centroidPlotList(centroidPlotList(:,colcount+1)==z,:);
%     plot(plotGreen(:,2),plotGreen(:,3),'g.')  % first column is centroid/object number, second column x, third column y
    plot(plotRed(:,2),plotRed(:,3),'r.')  
    if saveflag
        saveas(gcf,horzcat('img',num2str(z),'.jpg'),'jpg')
    end
end
