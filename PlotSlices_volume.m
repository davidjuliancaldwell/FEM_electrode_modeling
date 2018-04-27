%Plot CSF Gray and White matter for LL slices
%x increases from front (anterior) to back (posterior), dx=1 mm
%y increases from top of head (dorsal) to bottom (basal), dy=1 mm
% z increases from left to the right. Z is slice number, dz=0.9 mm
% xy is MR slice, [256x256]
% 4-25-2018 - David J. Caldwell 

clear all;clc
load LL_CSF.mat;
load LL_GrayMatter.mat;
load LL_WhiteMatter.mat;

% permute the order of x and y to make plotting easier
colOrder = [3 2 1];
CSF(:,colOrder) = flipud(CSF);
GrayMatter(:,colOrder) = flipud(GrayMatter);
WhiteMatter(:,colOrder) = flipud(WhiteMatter);

%%%%% first we will reshape it
jVec = [34:158]; % this iterates through the slices
x = [1:256];
y= [1:256];
z = [1:256];

% make a meshgrid to visualize the cube
[X,Y,Z] = meshgrid(x,y,z);

% make individual  matrices 
CSFReshape = zeros(size(X));
GrayMatterReshape = zeros(size(Y));
WhiteMatterReshape = zeros(size(Z));
BrainReshape = zeros(size(X));

% iterate through the j vector
for j=jVec
    
    % CSF
    % take a slice
    SliceNumber=j;
    
    % take z dimension of CSF
    J1=CSF(:,3); 
    J2=find(J1==SliceNumber);
    CSFSlice = CSF(J2,:);
    
    %convert to linear indices for indexing
    linearInd = sub2ind([256,256,256],CSFSlice(:,1),CSFSlice(:,2),CSFSlice(:,3));
    CSFReshape(linearInd) = 1;
    % assign CSF values of 1 
    BrainReshape(linearInd) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRAY MATTER
    J1=GrayMatter(:,3);
    J2=find(J1==SliceNumber);
    GrayMatterSlice = (GrayMatter(J2,:));
    linearInd = sub2ind([256,256,256],GrayMatterSlice(:,1),GrayMatterSlice(:,2),GrayMatterSlice(:,3));
    
    GrayMatterReshape(linearInd) = 1;
    BrainReshape(linearInd) = 3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WHITE MATTER
    J1=WhiteMatter(:,3);
    J2=find(J1==SliceNumber);
    WhiteMatterSlice = (WhiteMatter(J2,:));
    linearInd = sub2ind([256,256,256],WhiteMatterSlice(:,1),WhiteMatterSlice(:,2),WhiteMatterSlice(:,3));
    
    WhiteMatterReshape(linearInd) = 1;
    BrainReshape(linearInd) = 5;
    
end;

%
%BrainReshape = permute(BrainReshape,[3 1 2]);

%% Now plot the Cortex highlighting either the gray or white or CSF
% Highlight Gray Matter
% set light scale for visualization 
lightScale = 0.8;
figure;
set(gca, 'ZDir','reverse')
p = patch(isosurface(X,Y,Z,BrainReshape,1.25));
% get normals to the surface - this helps with lighting
isonormals(x,y,z,BrainReshape,p)

% check MATLAB version for compatability 
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

% daspect should take care of the 1 mm, 1mm, 0.9 mm scales for visualizing
% the x, y, and z 
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

% lighting additions 
lighting gouraud
material([.3 .8 .1 10 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Highlight White Matter
figure
set(gca, 'ZDir','reverse')
p = patch(isosurface(X,Y,Z,BrainReshape,3));
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
axis off
colormap(gray(50))
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


%% Now plot the cortex and a slice through it 

% ui box for input
prompt = {'X limit?','Y limit?','Z limit?'};
dlg_title = 'Subvolumes';
num_lines = 1;
defaultans = {'[NaN NaN]','[NaN NaN]','[55 NaN]'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

xLim = str2num(answer{1});
yLim = str2num(answer{2});
zLim = str2num(answer{3});

% get limits of the volume we are interested in
limits = [xLim(1) xLim(2) yLim(1) yLim(2) zLim(1) zLim(2)];

% extract a subset of the volume data
[x, y, z, D] = subvolume(BrainReshape, limits);          

 % isosurface for the outside of the volume
[fo,vo] = isosurface(x,y,z,D,1.7);          

% isocaps for the end caps of the volume
[fe,ve,ce] = isocaps(x,y,z,D,1.7);               

figure
set(gca, 'ZDir','reverse')

% draw the outside of the volume
hiso = patch('Faces', fo, 'Vertices', vo);

 % draw the end caps of the volume
hcap = patch('Faces', fe, 'Vertices', ve, ...   
    'FaceVertexCData', ce);

% check MATLAB version for compatability 

if verLessThan('matlab','8.5')
    set(hiso,'FaceColor',[0.7 0.7 0.7]);
    set(hiso,'EdgeColor','none');
    set(hcap,'FaceColor','interp');
    set(hcap,'EdgeColor','none');
else
    hiso.FaceColor = [0.7 0.7 0.7];
    %hiso.FaceAlpha = 0.7;
    hiso.EdgeColor = 'none';
    
    hcap.FaceColor = 'interp';
    hcap.EdgeColor = 'none';

end

colormap(gray(50))
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


