function batch3DYAdjustMC_AB

cd C:\jgk\cellCounting\batch\


clear variables
overlappix=[0 101 96 104 118 105 0];pialypos=402;wmPos=4796;zstepsize=0.5;
filenameintro='m652b542n148_V2_DKPG_box2_';
rangesa=struct(['box' num2str(1)],[20 23],['box' num2str(2)],[22 25],['box' num2str(3)],[24 27],...
    ['box' num2str(4)],[23 26],['box' num2str(5)],[24 27],['box' num2str(6)],[25 28]);
rangesb=struct(['box' num2str(1)],[75 81],['box' num2str(2)],[74 80],['box' num2str(3)],[76 82],...
    ['box' num2str(4)],[76 82],['box' num2str(5)],[79 85],['box' num2str(6)],[77 83]);
channeltags={'_GABA','_Kv3','_PV'}; channeltypes=[4 2 3];
% and change filename to save at bottom

pixadjust=1024*((1:length(overlappix))-1)-cumsum(overlappix(1:length(overlappix)))-pialypos;
sernumstr=filenameintro(end-1:end);

colocMata=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
colocMatb=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
chanImsa=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
chanImsb=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
imsa=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
imsb=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
centroids3da=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
centroids3db=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
colocMat=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[],['box' num2str(7)],[]);
colocMatCol=[];

for boxnum=1:length(overlappix)-1
    load(horzcat(filenameintro,num2str(boxnum),'_seg2d.mat'))
    centroids3da.(['box' num2str(boxnum)])=get3DCentroids(eval(horzcat('centroids',sernumstr,num2str(boxnum))),eval(horzcat('areas',sernumstr,num2str(boxnum))),eval(horzcat('pixflags',sernumstr,num2str(boxnum))),zstepsize,1,rangesa.(['box' num2str(boxnum)]));
    [colocMata.(['box' num2str(boxnum)]),chanImsa.(['box' num2str(boxnum)]),imsa.(['box' num2str(boxnum)])]=evalMultChannels(centroids3da.(['box' num2str(boxnum)]).stereo,{horzcat(filenameintro,num2str(boxnum),channeltags{1}),horzcat(filenameintro,num2str(boxnum),channeltags{2}),horzcat(filenameintro,num2str(boxnum),channeltags{3})},rangesa.(['box' num2str(boxnum)])(1):rangesa.(['box' num2str(boxnum)])(2),channeltypes,0.9);
    colocMata.((['box' num2str(boxnum)]))=colocMata.(['box' num2str(boxnum)])(colocMata.(['box' num2str(boxnum)])(:,3)<(1024-overlappix(boxnum+1)),:);
    
    centroids3db.(['box' num2str(boxnum)])=get3DCentroids(eval(horzcat('centroids',sernumstr,num2str(boxnum))),eval(horzcat('areas',sernumstr,num2str(boxnum))),eval(horzcat('pixflags',sernumstr,num2str(boxnum))),zstepsize,1,rangesb.(['box' num2str(boxnum)]));
    [colocMatb.(['box' num2str(boxnum)]),chanImsb.(['box' num2str(boxnum)]),imsb.(['box' num2str(boxnum)])]=evalMultChannels(centroids3db.(['box' num2str(boxnum)]).stereo,{horzcat(filenameintro,num2str(boxnum),channeltags{1}),horzcat(filenameintro,num2str(boxnum),channeltags{2}),horzcat(filenameintro,num2str(boxnum),channeltags{3})},rangesb.(['box' num2str(boxnum)])(1):rangesb.(['box' num2str(boxnum)])(2),channeltypes,0.9);
    colocMatb.(['box' num2str(boxnum)])=colocMatb.(['box' num2str(boxnum)])(colocMatb.(['box' num2str(boxnum)])(:,3)<(1024-overlappix(boxnum+1)),:);
    
    colocMat.(['box' num2str(boxnum)])=vertcat(colocMata.(['box' num2str(boxnum)]),colocMatb.(['box' num2str(boxnum)]));
    colocMat.(['box' num2str(boxnum)])(:,3)=colocMat.(['box' num2str(boxnum)])(:,3)+pixadjust(boxnum);
    colocMatCol=vertcat(colocMatCol,colocMat.(['box' num2str(boxnum)]));    
    colocMat.(['box' num2str(boxnum)])(:,3)=colocMat.(['box' num2str(boxnum)])(:,3)-pixadjust(boxnum);
    
end
save m652b542n148_V2_DKPG_ser2_12192019.mat


