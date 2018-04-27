function [threshold,titleName,specificLayer] = layer_choice_electrodes()

% threshold of
% 0.9 - for CSF
% 2.9 - for gray
% 4.9 -  for white
% These were set up in the generate_reshaped_brain.m function

answer = questdlg('Which layer would you like to pick electrodes on?', ...
    'layer pick electrodes', ...
    'CSF','Gray Matter','White Matter','Gray Matter');

% Handle response
switch answer
    case 'CSF'
        threshold = 0.9;
        specificLayer = 0.5;
        titleName = 'CSF';
    case 'Gray Matter'
        threshold = 2.9;
        specificLayer = 3;
        titleName = 'Gray Matter';
    case 'White Matter'
        threshold = 4.9;
        specificLayer = 5;
        titleName = 'White Matter';
end

end