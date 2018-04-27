function [layerChoice] = layer_choice_interactive(CSF,GrayMatter,WhiteMatter,iteration)

% Handle response
switch iteration
    
    case 1
        
        answer = questdlg('Pick the first layer to find the electrodes on', ...
            'interface', ...
            'CSF','Gray Matter','White Matter','Gray Matter');
        
        switch answer
            case 'CSF'
                layerChoice = CSF;
                
            case 'Gray Matter'
                layerChoice = GrayMatter;
                
            case 'White Matter'
                layerChoice = WhiteMatter;
        end
        
    case 2
        
        answer = questdlg('Pick a separate layer to find the boundary interface on', ...
            'interface', ...
            'CSF','Gray Matter','White Matter','CSF');
        
        switch answer
            case 'CSF'
                layerChoice = CSF;
                
            case 'Gray Matter'
                layerChoice = GrayMatter;
                
            case 'White Matter'
                layerChoice = WhiteMatter;
        end
end

end