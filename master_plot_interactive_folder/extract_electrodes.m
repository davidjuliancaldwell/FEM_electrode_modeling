function [electrodeLoc,closestVertex,normalVec] = extract_electrodes(BrainReshape,layerChoice,secondLayer,threshold,titleName,distanceWithin,x,y,z)


[figInteractive,pInteractive] = plot_layer_interest(BrainReshape,x,y,z,threshold,titleName);

dcm_obj = datacursormode(figInteractive);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on')

msgbox('Click brain to display a data tip, then press return.')
% Wait while the user does this.
pause

c_info = getCursorInfo(dcm_obj);
electrodeLocCursor = round(c_info.Position);

% must permute order because columns need to be flipped for 3D plots 
permuteOrder = [2 1 3];
%permuteOrder = [1 2 3];
coordsPerm = electrodeLocCursor(:,permuteOrder);
%
i = 1;
[idx,dist] = rangesearch(layerChoice,coordsPerm(i,:),distanceWithin);
firstLayerCoords = [layerChoice(idx{:},1),layerChoice(idx{:},2),layerChoice(idx{:},3)];
[idxInterface,distInterface] = rangesearch(firstLayerCoords,secondLayer,1.5);
intersectInds = find(~cellfun('isempty', idxInterface));
secondLayerPts = [secondLayer(intersectInds,1),secondLayer(intersectInds,2),secondLayer(intersectInds,3)];

DT_A= delaunayTriangulation(layerChoice(idx{:},1),layerChoice(idx{:},2),layerChoice(idx{:},3)); 
[F,P] = freeBoundary(DT_A);
[K,V] = convexHull(DT_A); 
TR = triangulation(F,P);
P_1 = incenter(TR);
F_1 = faceNormal(TR);

hold on  
ID = nearestNeighbor(TR,coordsPerm(i,:));
normalVec = vertexNormal(TR,ID);
closestVertex = TR.Points(ID,:);

% [figInteractive,pInteractive] = plot_layer_interest(BrainReshape,x,y,z,threshold,titleName);
% hold on
% scatter3(coordsPerm(:,2),coordsPerm(:,1),coordsPerm(:,3),100,'filled')
% scatter3(layerChoice(idx{:},2),layerChoice(idx{:},1),layerChoice(idx{:},3),10,'filled')
% quiver3(closestVertex(:,2),closestVertex(:,1),closestVertex(:,3), ...
%      normalVec(:,2),normalVec(:,1),normalVec(:,3),'LineWidth',5,'MaxHeadSize',40,'AutoScaleFactor',10,'color','g');

[newPt] = round(move_electrodes_along_normal(coordsPerm,normalVec,secondLayerPts));

electrodeLoc = newPt(permuteOrder);

[figInteractive,pInteractive] = plot_layer_interest(BrainReshape,x,y,z,threshold,titleName);
hold on
scatter3(newPt(:,2),newPt(:,1),newPt(:,3),100,'filled')
scatter3(layerChoice(idx{:},2),layerChoice(idx{:},1),layerChoice(idx{:},3),10,'filled')
quiver3(closestVertex(:,2),closestVertex(:,1),closestVertex(:,3), ...
     normalVec(:,2),normalVec(:,1),normalVec(:,3),'LineWidth',5,'MaxHeadSize',40,'AutoScaleFactor',10,'color','g');

end

