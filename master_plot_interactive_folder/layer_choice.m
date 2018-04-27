function [threshold,titleName,specificLayer] = layer_choice()

% threshold of
% 0.5 - for CSF
% 2 - for gray
% 4 -  for white
% These were set up in the generate_reshaped_brain.m function

answer = questdlg('Which layer would you like to visualize and see a cross section on?', ...
    'layer visualize', ...
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