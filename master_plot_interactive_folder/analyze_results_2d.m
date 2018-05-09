function [] = analyze_results_2d(resultsMatrix)

%Plot CSF Gray and White matter for LL slices
%x increases from front (anterior) to back (posterior), dx=1 mm
%y increases from top of head (dorsal) to bottom (basal), dy=1 mm
% z increases from left to the right. Z is slice number, dz=0.9 mm
% xy is MR slice, [256x256]


colormapUse = gray(6);
colormapUse = [colormapUse(end,:); colormapUse(2:end-1,:)];

for i = 1:size(electrodes,1)
    figure;
    colormap(colormapUse)
    hold on
    imagesc(squeeze(resultsMatrix(xSlice,:,:)))

    set(gca,'YDir','normal')
    title([' axial'])
    
    figure;
    colormap(colormapUse)
    hold on
    imagesc(squeeze(resultsMatrix(:,ySlice,:)))
    set(gca,'YDir','reverse')
    title([' coronal'])
    
    figure;
    colormap(colormapUse)
    hold on
    imagesc(squeeze(resultsMatrix(:,:,zSlice)))
    set(gca,'YDir','reverse')
    title([' sagital'])
end

end