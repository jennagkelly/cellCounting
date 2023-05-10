# cellCounting
Matlab implementation of "Quantification of neuronal density across cortical 
depth using automated 3D analysis of confocal image stacks" by JG Kelly and 
MJ Hawken (2017). Code by JG Kelly. 
The algorithm expects a z-stack of confocal images through immunolabeled 
tissue counterstained with DAPI and has been used with z-step sizes of 0.5-
1 um. Each wavelength-separated channel (DAPI and individual markers) should
be saved as its own tif image stack. Analysis runs in three stages: 

1. segmentObjectsByClustering preprocesses and segments the DAPI channel in
each (2D) image plane to get a preliminary estimate of the centroid of every
nuclear profile in each plane. 
Inputs are (1) the image stack filename, (2) the vector of z planes to 
segment, (3) the image plane farthest from the objective that has signal 
(used to set a threshold intensity level), (4) the z step size in um, and 
(5) true/false to plot the output or not. 
It is recommended that input 2 include planes above and below the exclusion
boundaries used for counting - the extra planes will be excluded in the 
subsequent step but are useful for the automated analysis.
Example call:
[images,centroids,areas,pixflags] = ...
	segmentObjectsByClustering('DAPIfilename.tif',4:58,4,0.5,1);

2. get3DCentroids uses the 2D segmentations from step #1 and groups 
centroids into clusters likely to be part of a single spherical object 
(nucleus). Computes a centroid in XYZ for each grouped object and assigns
object numbers.
Inputs 1-3 are outputs 2-4 from step #1. Remaining inputs are (4) the size
of one z step in microns, (5) a true/false flag to apply stereological 
exclusion boundaries in the XY plane, and optional (6) a 2-vector with the z 
planes limiting the analysis region (second specified z plane is excluded).
Example call:
centroids3d = get3DCentroids(centroids,areas,pixflags,0.5,1,[9,53]);
 
3. evalMultChannels uses the centroids from step #2 and evaluates whether
there is label at that location in images of a fluorescence tagged marker
collected in tandem with the DAPI channel. Can be used for up to 3 non-DAPI
channels at once and assigns each cell (centroid) a category based on its 
joint labeling pattern across channels:
Type 0: unlabeled (DAPI only)
Type 1: Channel 1 only
Type 2: Channel 2 only
Type 3: Channel 3 only
Type 4: Channels 1&2 
Type 5: Channels 1&3
Type 6: Channels 2&3
Type 7: Channels 1,2,&3 (all)
Inputs are (1) the "stereo" field of the output from step 2, (2) a cell 
array containing strings with the tif image stacks corresponding to each 
non-DAPI channel to be analyzed, (3) a vector of z planes to be analyzed 
(should match the range of input 6 in step #2 but list all intermediate 
planes), (4) a vector of channel types for each channel (see below), and 
(5) a criterion level for the fraction of pixels in the local region that 
should be "labeled" in order to classify a cell as labeled. A good starting
value is 0.9 but this may need to be adjusted based on performance.
Four channel type options (input 4) are included that differ in the 
preprocessing steps used to generate a binary image for evaluation. Type 1
was used to evaluate the NeuN channel in Kelly and Hawken (2017). Types 2, 
3, and 4 were used to evaluate the Kv3.1, PV and GABA channels, 
respectively, in Kelly et al. (2019).
The first output will be a matrix with one row per 3D centroid from step #2.
Columns 1-7 are inherited from step #2. The final n+1 columns (n=number of 
non-DAPI channels between 1 and 3) will contain positive (1)/negative (0) 
labels for each channel (per cell) and the joint labeling category (0-7, as 
above) in the last column.
Example call, for a single immuno channel similar to NeuN:
[coloc,chanIms,ims] = evalMultChannels(centroids3d.stereo,...
	{'channelfilename'},9:53,1,0.9);
Example call, for two immuno channels of different types:
[coloc,chanIms,ims] = evalMultChannels(centroids3d.stereo,...
	{'channel1filename','channel2filename'},9:53,[4,2],0.9);


References:
Kelly JG, Hawken MJ. 2017. Quantification of neuronal density across 
	cortical depth using automated 3D analysis of confocal image stacks.
	Brain Structure and Function 222, 3333-3353.
Kelly JG, García-Marín V, Rudy B, Hawken MJ. 2019. Densities and laminar 
	distributions of Kv3.1b-, PV-, GABA-, and SMI-32-immunoreactive 
	neurons in macaque area V1. Cerebral Cortex 29, 1921-1937.
