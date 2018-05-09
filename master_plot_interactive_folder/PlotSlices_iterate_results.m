%Plot CSF Gray and White matter for LL slices
%x increases from front (anterior) to back (posterior), dx=1 mm
%y increases from top of head (dorsal) to bottom (basal), dy=1 mm
% z increases from left to the right. Z is slice number, dz=0.9 mm
% xy is MR slice, [256x256]


% load in the different layers
load('FuxMagTest1.mat')
load('VoltTest1.mat')

%% go through z - right to left it looks like
figure;
for j=1:192
    subplot(2,1,1)
    imagesc(VoltMatrix(:,:,j));
    title('Volt Matrix')
    colorbar()
    
    subplot(2,1,2)
    imagesc(FluxMatrix(:,:,j));
    title('Flux Matrix')
    colorbar()
    
    pause(0.5)
    clf
end


%%
% go from front to back 
figure;
for j=1:256
    subplot(2,1,1)
    imagesc(squeeze(VoltMatrix(:,j,:)));
    title('Volt Matrix')
    colorbar()
    
    subplot(2,1,2)
    imagesc(squeeze(FluxMatrix(:,j,:)));
    title('Flux Matrix')
    colorbar()
    
    pause(0.5)
    clf
end


%%
% go from top to bottom
figure;
for j=1:256
    subplot(2,1,1)
    imagesc(squeeze(VoltMatrix(j,:,:)));
    title('Volt Matrix')
    colorbar()
    
    subplot(2,1,2)
    imagesc(squeeze(FluxMatrix(j,:,:)));
    title('Flux Matrix')
    colorbar()
    
    pause(0.5)
    clf
end

