function [fig,p] = plot_layer_interest(BrainReshape,x,y,z,threshold,titleName)

[X,Y,Z] = meshgrid(x,y,z);

lightScale = 0.8;
fig = figure;
set(gca, 'ZDir','reverse')
p = patch(isosurface(X,Y,Z,BrainReshape,threshold));
% get normals to the surface - this helps with lighting
isonormals(x,y,z,BrainReshape,p)

% check MATLAB version for compatability 
if verLessThan('matlab','8.5')
    set(p,'FaceColor',[0.8 0.8 0.8]);
    set(p,'EdgeColor','none');
    set(p,'FaceAlpha',1);
else
    p.FaceColor = [0.8 0.8 0.8];
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.8;
end
view(3);
%axis tight
box off
axis off
title(titleName)
% daspect should take care of the 1 mm, 1mm, 0.9 mm scales for visualizing
% the x, y, and z 
daspect([1 1 0.9])
view(90,90);
l = camlight('headlight');
set(l,'color',[lightScale lightScale lightScale]*0.8);
view(90,-90);
l = camlight('headlight');
set(l,'color',[lightScale lightScale lightScale]*0.8);
view(90,0);
l = camlight('headlight');
set(l,'color',[lightScale lightScale lightScale]*.3);
view(-90,0);
l = camlight('headlight');
set(l,'color',[lightScale lightScale lightScale]*.3);

az = 0;
el = -90;
view(az, el);

% lighting additions 
lighting gouraud
material([.3 .8 .1 10 1]);

end

