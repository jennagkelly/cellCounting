function plot3dMultiChannel(colocMatIN,imagesIN,saveflag,plotZ,plotcol)

colcount=size(colocMatIN,2);
centroidPlotList=horzcat(colocMatIN,nan(size(colocMatIN,1),1));

for ind=1:size(centroidPlotList,1)
    centroidPlotList(ind,colcount+1)=round(centroidPlotList(ind,4));
end

cd 'C:\jgk\cellCounting\plottemp'

for z=plotZ
    if isa(imagesIN{z},'logical')
        figure('visible','off');imshow(imagesIN{z},[]);
        iptsetpref('ImshowBorder','tight')
    else
        figure('visible','off');imshow(imadjust(imagesIN{z},stretchlim(imagesIN{z}),[]));
        iptsetpref('ImshowBorder','tight')
    end
    hold on;
    thisZList=centroidPlotList(centroidPlotList(:,colcount+1)==z,:);
    plotGreen=thisZList(thisZList(:,plotcol)==0,:);
    plotRed=thisZList(thisZList(:,plotcol)==1,:);
    if saveflag
        plot(plotGreen(:,2),plotGreen(:,3),'g.','MarkerSize',2)  
        plot(plotRed(:,2),plotRed(:,3),'r.','MarkerSize',2) 
        set(gcf,'PaperPosition',[0 0 1 1])
        print(horzcat('MCimg',num2str(z),'.jpg'),'-djpeg','-r1024')
%         saveas(gcf,horzcat('MCimg',num2str(z),'.jpg'),'jpg')
    else
        plot(plotGreen(:,2),plotGreen(:,3),'g.')  
        plot(plotRed(:,2),plotRed(:,3),'r.') 
    end
end
