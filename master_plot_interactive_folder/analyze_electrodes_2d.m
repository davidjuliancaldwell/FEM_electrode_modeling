function [] = analyze_electrodes_2d(BrainReshape,CSF,GrayMatter,WhiteMatter,electrodes)

%Plot CSF Gray and White matter for LL slices
%x increases from front (anterior) to back (posterior), dx=1 mm
%y increases from top of head (dorsal) to bottom (basal), dy=1 mm
% z increases from left to the right. Z is slice number, dz=0.9 mm
% xy is MR slice, [256x256]

permuteOrder = [2 1 3];
electrodes = electrodes(:,permuteOrder);

xSlice = electrodes(1,1);
ySlice = electrodes(1,2);
zSlice = electrodes(1,3);

colormapUse = gray(6);
colormapUse = [colormapUse(end,:); colormapUse(2:end-1,:)];

for i = 1:size(electrodes,1)
    figure;
    colormap(colormapUse)
    hold on
    imagesc(squeeze(BrainReshape(xSlice,:,:)))

    scatter(electrodes(i,3),electrodes(i,2),'Filled')
    set(gca,'YDir','normal')
    title(['electrode ' num2str(i) ' axial'])
    
    figure;
    colormap(colormapUse)
    hold on
    imagesc(squeeze(BrainReshape(:,ySlice,:)))
    scatter(electrodes(i,3),electrodes(i,1),'Filled')
    set(gca,'YDir','reverse')
    title(['electrode ' num2str(i) ' coronal'])
    
    figure;
    colormap(colormapUse)
    hold on
    imagesc(squeeze(BrainReshape(:,:,zSlice)))
    scatter(electrodes(i,2),electrodes(i,1),'filled')
    set(gca,'YDir','reverse')
    title(['electrode ' num2str(i) ' sagital'])
end

end