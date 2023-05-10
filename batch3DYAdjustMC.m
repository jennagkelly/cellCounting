function batch3DYAdjustMC

clear variables
overlappix=[0 84 77 89 66 61 88 60 0];pialypos=84;wmPos=6943;zstepsize=0.5;
filenameintro='M657b741n74a_MT_DKNG_box1_';
ranges=struct(['box' num2str(1)],[10 85],['box' num2str(2)],[12 87],['box' num2str(3)],[9 84],...
    ['box' num2str(4)],[7 82],['box' num2str(5)],[10 85],['box' num2str(6)],[9 84],...
    ['box' num2str(7)],[12 87],['box' num2str(8)],[8 83]);
channeltags={'_NeuN_2'}; channeltypes=1;
% and change filename to save at bottom

pixadjust=1024*((1:length(overlappix))-1)-cumsum(overlappix(1:length(overlappix)))-pialypos;
sernumstr=filenameintro(end-1:end);

colocMat=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[]);
chanIms=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[]);
ims=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[]);
centroids3d=struct(['box' num2str(1)],[],['box' num2str(2)],[],['box' num2str(3)],[],['box' num2str(4)],[],['box' num2str(5)],[],['box' num2str(6)],[]);
colocMatCol=[];

for boxnum=1:length(overlappix)-1
    load(horzcat(filenameintro,num2str(boxnum),'_seg2d.mat'))
    centroids3d.(['box' num2str(boxnum)])=get3DCentroids(eval(horzcat('centroids',sernumstr,num2str(boxnum))),eval(horzcat('areas',sernumstr,num2str(boxnum))),eval(horzcat('pixflags',sernumstr,num2str(boxnum))),zstepsize,1,ranges.(['box' num2str(boxnum)]));
    [colocMat.(['box' num2str(boxnum)]),chanIms.(['box' num2str(boxnum)]),ims.(['box' num2str(boxnum)])]=evalMultChannels(centroids3d.(['box' num2str(boxnum)]).stereo,{horzcat(filenameintro,num2str(boxnum),channeltags{1})},ranges.(['box' num2str(boxnum)])(1):ranges.(['box' num2str(boxnum)])(2),channeltypes,0.9);
    colocMat.((['box' num2str(boxnum)]))=colocMat.(['box' num2str(boxnum)])(colocMat.(['box' num2str(boxnum)])(:,3)<(1024-overlappix(boxnum+1)),:);
    colocMat.(['box' num2str(boxnum)])(:,3)=colocMat.(['box' num2str(boxnum)])(:,3)+pixadjust(boxnum);
      
    colocMatCol=vertcat(colocMatCol,colocMat.(['box' num2str(boxnum)]));  
    colocMat.(['box' num2str(boxnum)])(:,3)=colocMat.(['box' num2str(boxnum)])(:,3)-pixadjust(boxnum);
   
end

save M657b741n74a_MT_DKNG_ser1_DN_12192018.mat

