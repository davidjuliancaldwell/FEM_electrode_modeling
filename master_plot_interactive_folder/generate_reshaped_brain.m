function [BrainReshape,CSFReshape,WhiteMatterReshape,GrayMatterReshape] =  generate_reshaped_brain(CSF,GrayMatter,WhiteMatter)

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

end

