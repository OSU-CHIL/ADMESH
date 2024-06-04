function CloseADMESH(varargin)

msg = 'Are you sure you want to close ADMESH?';
choice = uiconfirm(app.UIFigure,msg,'ADMESH',...
    'Options',{'Yes','No'},'DefaultOption',2,'Icon','Warning');
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