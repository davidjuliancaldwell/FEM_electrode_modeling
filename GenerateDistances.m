%% 4-25-2018 - script to generate distances for Ceon's modeling
% if looking at a 3d plot, the x and y coordinates are reversed!

generate_gray
%%





fixedD = 5;
coords1 = [138 66.13 40.13];
coords2 = [150 66 41.13];
permuteOrder = [2 1 3];
gcf
hold on
coords = [coords1; coords2];
scatter3(coords(:,1),coords(:,2),coords(:,3),100,'filled')
coordsPerm = coords(:,permuteOrder);

difference = coords1 - coords2;
distance = sqrt(difference * difference')
%%
i = 1
[idx,dist] = rangesearch(GrayMatter,coordsPerm(i,:),fixedD)
DT_A= delaunayTriangulation(GrayMatter(idx{:},1),GrayMatter(idx{:},2),GrayMatter(idx{:},3)); 
figure
[F,P] = freeBoundary(DT_A);
trisurf(F,P(:,1),P(:,2),P(:,3), ...
       'FaceColor','cyan','FaceAlpha',0.8);

[K,V] = convexHull(DT_A); 
figure;
%trisurf(tri,DT_A.Points(:,1),DT_A.Points(:,2),DT_A.Points(:,3), 'FaceColor', 'cyan', 'faceAlpha', 0.8); 
trisurf(K,DT_A.Points(:,1),DT_A.Points(:,2),DT_A.Points(:,3)); 

hold on
scatter3(coordsPerm(i,1),coordsPerm(i,2),coordsPerm(i,3),'filled')

TR = triangulation(F,P);
P_1 = incenter(TR);
F_1 = faceNormal(TR);
hold on  
quiver3(P_1(:,1),P_1(:,2),P_1(:,3), ...
     F_1(:,1),F_1(:,2),F_1(:,3),0.5,'color','r');
 
ID = nearestNeighbor(TR,coordsPerm(i,:));
Vertex = vertexNormal(TR,ID);
Point = TR.Points(ID,:);
scatter3(Point(i,1),Point(i,2),Point(i,3),'filled')

%%
figure
lightScale = 0.8;

figure;
set(gca, 'ZDir','reverse')
p = patch(isosurface(X,Y,Z,BrainReshape,1.25));
isonormals(x,y,z,BrainReshape,p)
if verLessThan('matlab','8.5')
    set(p,'FaceColor',[0.8 0.8 0.8]);
    set(p,'EdgeColor','none');
else
    p.FaceColor = [0.8 0.8 0.8];
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.7;
end
view(3);
%axis tight
box off
axis off
daspect([1 1 0.9])
view(90,90);
l = camlight('headlight');
set(l,'color',[lightScale lightScale lightScale]*.3);
view(90,-90);
l = camlight('headlight');
set(l,'color',[lightScale lightScale lightScale]*.3);
view(90,0);
l = camlight('headlight');
set(l,'color',[lightScale lightScale lightScale]);
view(-90,0);
l = camlight('headlight');
set(l,'color',[lightScale lightScale lightScale]);

lighting gouraud
material([.3 .8 .1 10 1]);
hold on
scatter3(coords(:,1),coords(:,2),coords(:,3),100,'filled')
scatter3(GrayMatter(idx{:},2),GrayMatter(idx{:},1),GrayMatter(idx{:},3),50,'filled')
quiver3(Point(:,2),Point(:,1),Point(:,3), ...
     Vertex(:,2),Vertex(:,1),Vertex(:,3),'LineWidth',5,'MaxHeadSize',40,'AutoScaleFactor',10,'color','g');



%%
i = 1
[idx,dist] = rangesearch(CSF,coords(i,:),fixedD)
gcf
DT_A= delaunayTriangulation(CSF(idx{:},1),CSF(idx{:},2),CSF(idx{:},3));

[K,V] = convexHull(DT_A); 
%figure;
%trisurf(tri,DT_A.Points(:,1),DT_A.Points(:,2),DT_A.Points(:,3), 'FaceColor', 'cyan', 'faceAlpha', 0.8); 
trisurf(K,DT_A.Points(:,1),DT_A.Points(:,2),DT_A.Points(:,3), 'FaceColor', 'cyan', 'faceAlpha', 0.8); 

hold on
scatter3(coords(i,1),coords(i,2),coords(i,3),20,'filled')


    
%%  Find values for projection
%   This function is called in ProjectElectrodes() to find the minimum
%   distance on a surface "gs" that is closest for each electrode in "els"

function out_els = p_zoom(els, gs)
%function out_els=p_zoom(els, gs, index);
% this function finds the minimum distance on a surface "gs" (pts x 3) that is
% closest for each electrode in "els" (N x 3)
% index indicates the number of point for calculation of norm, 0 if global
% checkdistance = 1 to indicate that electrodes within 3 mm of the surface
% are not projected

%use a local set of electrodes to determine the orthogonal direction of a
%given electrode? -- enter "0" for global, and number to use otherwise
%   Created by:
%   D. Hermes & K.J. Miller
%   Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%   Version 1.1.0, released 26-11-2009

checkdistance_dist = 3; % 3 mm

for k = 1:size(els,1)
    
    if find(ignoreTrodes == k,1,'first')>0
        continue;
    end

    npls = [gs(:,1)-els(k,1) gs(:,2)-els(k,2) gs(:,3)-els(k,3)]; %x,y,z lengths
    % calculate distance
    npls_dist = sqrt(sum((npls).^2,2));
    % check whether distance is < 3 mm
    distancesm2 = 0;
    if npls_dist(find(npls_dist == min(npls_dist)),:) < checkdistance_dist
        %disp(['distance < 3 mm electrode ' int2str(k) ]);
        distancesm2 = 1;
    end
    
    if checkdistance == 1 && distancesm2 == 1 % electrode too close to surface to project
        out_ind(k) = find(npls_dist == min(npls_dist),1); %find minimum distance
    elseif checkdistance == 2
        out_ind(k) = find(npls_dist == min(npls_dist),1); %find minimum distance
    else
        npls_unit = npls./repmat((sum(npls.^2,2).^.5),1,3); % normalize npls to get unit vector
        npdot = (npls_unit*nm); %distance along eigenvector direction (dot product)
        % only take gs within 2.5 cm distance
        npdot(npls_dist>25) = 0;
        %npdotrev=(npls_unit*-nm); % distance in reverse eigenvector direction
        [a b] = find(abs(npdot) == max(max(abs(npdot))),1);
        out_ind(k) = a;%find minimum distance, max dot product
        %out_ind_rev(k)=find(npdotrev==max(npdotrev),1); %find minimum distance, max dot product
    end
end
out_ind(out_ind == 0) = [];
out_els = els;
out_els(setdiff(1:size(els,1),ignoreTrodes),:) = gs(out_ind,:);
% out_els_rev=gs(out_ind_rev,:);
% % check distance to new electrodes
% out_els_dist=sqrt(sum((els-out_els).^2,2));
% out_els_rev_dist=sqrt(sum((els-out_els_rev).^2,2));

% plot on surface to check
figure
plot3(els(:,1),els(:,2),els(:,3),'r.','MarkerSize',20);
hold on;
plot3(out_els(:,1),out_els(:,2),out_els(:,3),'g.','MarkerSize',20);
plot3(gs(:,1),gs(:,2),gs(:,3),'k.','MarkerSize',1);
axis equal;

end
