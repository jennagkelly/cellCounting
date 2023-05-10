function applyManCorrectionsa(wkspctofix,chnls)
directory='C:\jgk\cellCounting\batch';
savedir='C:\jgk\cellCounting\batch';
cd(directory)
load(wkspctofix)

numChannels=length(chnls);
colocMatCol=[];
for boxnum=1:length(overlappix)-1
    for chnl=1:length(chnls)
        filename=strcat('CellCounter_auto_',filenameintro,num2str(boxnum),'a_',chnls{chnl},'.xml');
        markers=imageJCellCounterDataToMLBStruct(directory,filename,1);
        for marker=1:size(markers,1)
            zmatch=colocMat.(['box' num2str(boxnum)])(round(colocMat.(['box' num2str(boxnum)])(:,4))==rangesa.(['box' num2str(boxnum)])(1)+markers(marker,3)-1,:);
            dist=sqrt((zmatch(:,2)-markers(marker,1)).^2+(zmatch(:,3)-markers(marker,2)).^2);
            [~,mindist]=min(dist);
            colocMat.(['box' num2str(boxnum)])(colocMat.(['box' num2str(boxnum)])(:,1)==zmatch(mindist,1),size(colocMat.(['box' num2str(boxnum)]),2)-1-numChannels+chnl)=...
                abs(colocMat.(['box' num2str(boxnum)])(colocMat.(['box' num2str(boxnum)])(:,1)==zmatch(mindist,1),size(colocMat.(['box' num2str(boxnum)]),2)-1-numChannels+chnl)-1);
        end
    end
    colocMat.(['box' num2str(boxnum)])(:,3)=colocMat.(['box' num2str(boxnum)])(:,3)+pixadjust(boxnum);
    colocMatCol=vertcat(colocMatCol,colocMat.(['box' num2str(boxnum)]));
    colocMat.(['box' num2str(boxnum)])(:,3)=colocMat.(['box' num2str(boxnum)])(:,3)-pixadjust(boxnum);
end


for row=1:size(colocMatCol,1)
    singleChannelAssignments=colocMatCol(row,size(colocMatCol,2)-numChannels:size(colocMatCol,2)-numChannels+2);
    if sum(singleChannelAssignments)<1
        colocMatCol(row,size(colocMatCol,2))=0;
    elseif sum(singleChannelAssignments)<2
        colocMatCol(row,size(colocMatCol,2))=find(singleChannelAssignments);
    elseif sum(singleChannelAssignments)<3
        posChans=sum(find(singleChannelAssignments)); % 1+2=3 (type4), 1+3=4 (type5), 2+3=5(type6)
        colocMatCol(row,size(colocMatCol,2))=posChans+1;
    elseif sum(singleChannelAssignments)==3   
        colocMatCol(row,size(colocMatCol,2))=7;
    end
end

newfilename=strcat(wkspctofix(1:length(wkspctofix)-4),'CORR.mat');
clear wkspctofix directory
cd(savedir)
save(newfilename)