function markers=imageJCellCounterDataToMLBStruct(directory,filename,markerTypes)
% % by JK 1/12/2015
% % takes Cell Counter plugin's saved marker file (XML doc) and extracts
% %     marker XYZ position data, with marker type (1-8) specified
% % example call: markerMat=imageJCellCounterDataToMLBStruct(pwd,filename,[1 5]);

cd(directory)
cellCounterDataFull=xmlread(filename);
cellCounterMarkerFileNode=cellCounterDataFull.getDocumentElement;
ccmfChildNodes=cellCounterMarkerFileNode.getChildNodes;
markerDataNodes=ccmfChildNodes.item(3).getChildNodes;

numMarkersAll=0;
for type=markerTypes
    thisTypeNodes=markerDataNodes.item(type*2+1).getChildNodes;
%     thisTypeNodes.item(1).getTextContent   % check type: should return correct markerTypes val
    numMarkersAll=numMarkersAll+thisTypeNodes.getElementsByTagName('Marker').getLength;
end
markers=nan(numMarkersAll,4);  % init

markerTally=1;
for type=markerTypes
    thisTypeNodes=markerDataNodes.item(type*2+1).getChildNodes;
    markerType=str2double(char(thisTypeNodes.item(1).getTextContent));   % check type: should return correct markerTypes val
    numThisTypeMarkers=thisTypeNodes.getElementsByTagName('Marker').getLength;
    for thisTypeMarker=0:numThisTypeMarkers-1
        markers(markerTally,1)=str2double(char(thisTypeNodes.getElementsByTagName('MarkerX').item(thisTypeMarker).getTextContent));
        markers(markerTally,2)=str2double(char(thisTypeNodes.getElementsByTagName('MarkerY').item(thisTypeMarker).getTextContent));
        markers(markerTally,3)=str2double(char(thisTypeNodes.getElementsByTagName('MarkerZ').item(thisTypeMarker).getTextContent));
        markers(markerTally,4)=markerType;   
        markerTally=markerTally+1;
    end    
end   


