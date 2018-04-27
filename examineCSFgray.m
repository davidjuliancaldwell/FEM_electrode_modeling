close all;clear all;clc 

load('LL_CSF.mat')
load('LL_GrayMatter.mat')

%%
image(


%Qucik Program to demo the use of findPointNormals

%generate a set of 3d points
x = repmat(1:49,49,1);
y = x';
z = peaks;
points = [x(:),y(:),z(:)];

%find the normals and curvature
[normals,curvature] = findPointNormals(CSF,[],[0,0,10],true);

%plot normals and colour the surface by the curvature
hold off;
surf(CSF(:,1),CSF(:,2),CSF(:,3)),reshape(curvature,256,256));
hold on;
quiver3(points(:,1),points(:,2),points(:,3),...
    normals(:,1),normals(:,2),normals(:,3),'r');
axis equal;

%%



figure
colormap(map)
Ds = smooth3(D);
hiso = patch(isosurface(Ds,5),...
   'FaceColor',[1,.75,.65],...
   'EdgeColor','none');
   isonormals(Ds,hiso)