function [distanceWithin,distanceFixed] = distance_choice()

% ui box for input
prompt = {'distance of points to search within for? (mm)','distance range? (mm)'};
dlg_title = 'point search';
num_lines = 1;
defaultans = {'10','0.25'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

distanceWithin = str2num(answer{1});
distanceFixed = str2num(answer{2});

end