function [moveDist] = move_electrodes()

% ui box for input
prompt = {'how many mm to move electrodes along normal vector ? (mm)'};
dlg_title = 'vector move';
num_lines = 1;
defaultans = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

moveDist = str2num(answer{1});

end