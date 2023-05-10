function batch2dSegment

cd C:\jgk\cellCounting\batch\
clear variables

% [images1_1,centroids1_1,areas1_1,pixflags1_1] = segmentObjectsByClustering('m652b542n148_V2_DKPG_box2_1_DAPI.tif',16:85,16,0.5,0);
% save m652b542n148_V2_DKPG_box2_1_seg2d.mat
% clear variables
% [images1_2,centroids1_2,areas1_2,pixflags1_2] = segmentObjectsByClustering('m652b542n148_V2_DKPG_box2_2_DAPI.tif',18:84,18,0.5,1);
% save m652b542n148_V2_DKPG_box2_2_seg2d.mat
% clear variables
[images1_3,centroids1_3,areas1_3,pixflags1_3] = segmentObjectsByClustering('m652b542n148_V2_DKPG_box2_3_DAPI.tif',20:86,20,0.5,0);
save m652b542n148_V2_DKPG_box2_3_seg2d.mat
clear variables
[images1_4,centroids1_4,areas1_4,pixflags1_4] = segmentObjectsByClustering('m652b542n148_V2_DKPG_box2_4_DAPI.tif',19:86,19,0.5,0);
save m652b542n148_V2_DKPG_box2_4_seg2d.mat
clear variables
[images1_5,centroids1_5,areas1_5,pixflags1_5] = segmentObjectsByClustering('m652b542n148_V2_DKPG_box2_5_DAPI.tif',20:89,20,0.5,0);
save m652b542n148_V2_DKPG_box2_5_seg2d.mat
clear variables
[images1_6,centroids1_6,areas1_6,pixflags1_6] = segmentObjectsByClustering('m652b542n148_V2_DKPG_box2_6_DAPI.tif',21:87,21,0.5,0);
save m652b542n148_V2_DKPG_box2_6_seg2d.mat
clear variables


