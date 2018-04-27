%Plot CSF Gray and White matter for LL slices
%x increases from front (anterior) to back (posterior), dx=1 mm
%y increases from top of head (dorsal) to bottom (basal), dy=1 mm
% z increases from left to the right. Z is slice number, dz=0.9 mm
% xy is MR slice, [256x256]

load LL_CSF.mat;
load LL_GrayMatter.mat;
load LL_WhiteMatter.mat;

for j=36:85;
SliceNumber=j;
J1=CSF(:,3);
J2=find(J1==SliceNumber);
CSFSlice=(CSF(J2,:));

J1=GrayMatter(:,3);
J2=find(J1==SliceNumber);
GrayMatterSlice=(GrayMatter(J2,:));

J1=WhiteMatter(:,3);
J2=find(J1==SliceNumber);
WhiteMatterSlice=(WhiteMatter(J2,:));

subplot(5,10,j-35);plot(CSFSlice(:,1), CSFSlice(:,2), 'G.');
hold on
subplot(5,10,j-35);plot(GrayMatterSlice(:,1), GrayMatterSlice(:,2), 'R.');
subplot(5,10,j-35);plot(WhiteMatterSlice(:,1), WhiteMatterSlice(:,2), 'B.');
end;