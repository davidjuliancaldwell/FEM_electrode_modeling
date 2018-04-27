% This is a script to viusalize a segmented image with CSF, Gray, and White
% Matter in 3-dimensions 
%
% David J. Caldwell 
%
% Requires: 
% CSF - matrix of n x 3 points
% GrayMatter - matrix of n x 3 points
% WhiteMatter - matrix of n x 3 points 
%
% x increases from front (anterior) to back (posterior), dx=1 mm
% y increases from top of head (dorsal) to bottom (basal), dy=1 mm
% z increases from left to the right. Z is slice number, dz=0.9 mm
% xy is MR slice, [256x256]

%% clear workspace, load data 
clear all;clc

% load in the different layers
load LL_CSF.mat;
load LL_GrayMatter.mat;
load LL_WhiteMatter.mat;

% define the coordinate space of x, y and z
x = [1:256];
y= [1:256];
z = [1:256];

%% Generate matrices for data exploration and plotting

% BrainReshape - this is a 256 x 256 x 256 matrix, with values of:
% 1 @ csf locations,
% 3 @ gray matter locations
% 5 @ white matter locations
%
% CSF reshape - this is only the CSF
% WhiteMatterReshape - this is only the white matter
% GrayMatterReshape - this is only the gray matter

[BrainReshape,CSFReshape,WhiteMatterReshape,GrayMatterReshape] =  generate_reshaped_brain(CSF,...
    GrayMatter,WhiteMatter);

%% Plot either the CSF, Gray Matter, or White Matter exterior

% threshold of:
% 0.5 - for CSF 
% 2 - for gray
% 4 -  for white
% These were set up in the generate_reshaped_brain.m function, as shown
% above

% ask user for the layer of interest
[threshold,titleName] = layer_choice;

% plot the layer of interest
[fig,p] = plot_layer_interest(BrainReshape,x,y,z,threshold,titleName);

%% Plot different cross sections of the brain 

% ask user for the layer of interest
[threshold,titleName] = layer_choice;

% plot the layer of interest, and also ask the user for the dimensions of
% the cross section. The first element in the vector sets the lower slice
% of the image, while the second sets the upper slice. 
% e.g., Y lims of [55, NaN] puts a cut at 55 mm from the top, and
% visualizes this slice, with the rest of the image shown as a 3D volume 

[fig,hiso,hcap] = plot_brain_cross_section(BrainReshape,threshold);

%% Interactively
% 1) select points on the layer of interest
% 2) find all points on that layer within a certain distance
% 3) find the normal vector from the nearest vertex 
% 4) Export data to workspace to have electrode locations 
% 5) NOTE: Clinical ECoG electrodes are 2.3 mm^2 in diameter, with 10 mm
% spacing between them 

% ask user for the layer of interest
% probably should pick gray matter for both of these 

[threshold,titleName] = layer_choice_electrodes;
[layerChoice] = layer_choice_interactive(CSF,GrayMatter,WhiteMatter,1);

% figure out the distances
% distanceWithin - node distance to select. e.g 10 mm
% distanceFixed - bracket range on either side of distanceWithin
% e.g 0.5 would yield 9.5 mm to 10.5, as in a ring of vertices
[distanceWithin,distanceRange] = distance_choice;

%% extract the electrodes given the options above 

% choose the second layer. e.g. if gray above, probably pick CSF to find
% electrodes on CSF/Gray boundary 
[secondLayer] = layer_choice_interactive(CSF,GrayMatter,WhiteMatter,2);

[electrodeLocFirst,closestVertexFirst,normalVecFirst] = extract_electrodes(BrainReshape,...
    layerChoice,secondLayer,threshold,titleName,distanceWithin,x,y,z);

%% find another electrode with specific distance 

[electrodeLocSecond,closestVertexSecond,normalVecSecond] = find_additional_electrode(BrainReshape,electrodeLocFirst,...
    layerChoice,secondLayer,threshold,titleName,distanceWithin,distanceRange,x,y,z);

%% plot both electrodes

[fig,p] = plot_layer_interest(BrainReshape,x,y,z,threshold,titleName);
hold on
scatter3(electrodeLocFirst(:,1),electrodeLocFirst(:,2),electrodeLocFirst(:,3),100,'filled')
quiver3(electrodeLocFirst(:,1),electrodeLocFirst(:,2),electrodeLocFirst(:,3), ...
     normalVecFirst(:,2),normalVecFirst(:,1),normalVecFirst(:,3),'LineWidth',5,'MaxHeadSize',40,'AutoScaleFactor',10,'color','g');
 
scatter3(electrodeLocSecond(:,1),electrodeLocSecond(:,2),electrodeLocSecond(:,3),100,'filled')
quiver3(electrodeLocSecond(:,1),electrodeLocSecond(:,2),electrodeLocSecond(:,3), ...
     normalVecSecond(:,2),normalVecSecond(:,1),normalVecSecond(:,3),'LineWidth',5,'MaxHeadSize',40,'AutoScaleFactor',10,'color','g');

% print out distance to make sure they are pretty close
distance2points = pdist2(electrodeLocFirst,electrodeLocSecond)
%% look at 2D cross sections for each of these 
electrodeTotal = ([electrodeLocFirst;electrodeLocSecond]);
analyze_electrodes_2d(BrainReshape,CSF,GrayMatter,WhiteMatter,electrodeTotal)

%% repermute to get in the appropriate position for GrayMatter,WhiteMatter, and CSF .mat files
permuteOrder = [2 1 3];
electrodes = electrodeTotal(:,permuteOrder);

% verify
min(pdist2(electrodes(1,:),CSF)) % check if close to CSF
min(pdist2(electrodes(1,:),GrayMatter)) % check if close to Gray Matter

%% if want to move electrodes along normal line additionally 

distMove = move_electrodes; % move 1 mm
electrodeLocFirstMoved = round(electrodeLocFirst + distMove*permute(normalVecFirst,[2 1 3])');
electrodeLocSecondMoved = round(electrodeLocSecond + distMove*permute(normalVecSecond,[2 1 3])');
figure
[fig,p] = plot_layer_interest(BrainReshape,x,y,z,threshold,titleName);
hold on
scatter3(electrodeLocFirstMoved(:,1),electrodeLocFirstMoved(:,2),electrodeLocFirstMoved(:,3),100,'filled')
quiver3(electrodeLocFirstMoved(:,1),electrodeLocFirstMoved(:,2),electrodeLocFirstMoved(:,3), ...
     normalVecFirst(:,2),normalVecFirst(:,1),normalVecFirst(:,3),'LineWidth',5,'MaxHeadSize',40,'AutoScaleFactor',10,'color','g');
 
scatter3(electrodeLocSecondMoved(:,1),electrodeLocSecondMoved(:,2),electrodeLocSecondMoved(:,3),100,'filled')
quiver3(electrodeLocSecondMoved(:,1),electrodeLocSecondMoved(:,2),electrodeLocSecondMoved(:,3), ...
     normalVecSecond(:,2),normalVecSecond(:,1),normalVecSecond(:,3),'LineWidth',5,'MaxHeadSize',40,'AutoScaleFactor',10,'color','g');

% print out distance to make sure they are pretty close
distance2points = pdist2(electrodeLocFirstMoved,electrodeLocSecondMoved)

%% look at 2D cross sections for each of these 
electrodeTotalMoved = ([electrodeLocFirstMoved;electrodeLocSecondMoved]);
analyze_electrodes_2d(BrainReshape,CSF,GrayMatter,WhiteMatter,electrodeTotalMoved)

%% repermute to get in the appropriate position for GrayMatter,WhiteMatter, and CSF .mat files
permuteOrder = [2 1 3];
electrodesMoved = electrodeTotalMoved(:,permuteOrder);

% verify
min(pdist2(electrodesMoved(1,:),CSF)) % check if close to CSF
min(pdist2(electrodesMoved(1,:),GrayMatter)) % check if close to Gray Matter


