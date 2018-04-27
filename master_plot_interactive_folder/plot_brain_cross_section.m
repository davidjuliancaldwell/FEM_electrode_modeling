function [fig,hiso,hcap] = plot_brain_cross_section(BrainReshape,threshold)

% ui box for input
prompt = {'X limit?','Y limit?','Z limit?'};
dlg_title = 'Subvolumes';
num_lines = 1;
defaultans = {'[NaN NaN]','[55 NaN]','[NaN NaN]'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

xLim = str2num(answer{1});
yLim = str2num(answer{2});
zLim = str2num(answer{3});

% get limits of the volume we are interested in
limits = [xLim(1) xLim(2) yLim(1) yLim(2) zLim(1) zLim(2)];

% extract a subset of the volume data
[x, y, z, D] = subvolume(BrainReshape, limits);          

 % isosurface for the outside of the volume
[fo,vo] = isosurface(x,y,z,D,threshold);          

% isocaps for the end caps of the volume
[fe,ve,ce] = isocaps(x,y,z,D,threshold);               

fig = figure;
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
lightScale = 0.8;
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

az = 0;
el = -90;
view(az, el);

lighting gouraud
material([.3 .8 .1 10 1]);


end

