function CloseADMESH(varargin)

choice = questdlg('Are you sure you want to close ADMESH?','ADMESH',...
    'Yes','No','No');
drawnow; pause(.005);

switch choice
    
    case 'Yes'
        warning_status = warning;
        warning('off');
        rmpath(genpath('libraries'));
        warning(warning_status);
        delete(varargin{1})

    case 'No'
        
        
end



end