function AssignBC(app,varargin)

%------------------------------------------------------------------------------
% Get line handle
%------------------------------------------------------------------------------
fig    = varargin{1};
handle = gco;

%------------------------------------------------------------------------------
% Get user data
%------------------------------------------------------------------------------
UserData = get(handle,'UserData');

switch UserData{1}
    
    case 'External Boundary'
                
        % Set the button down function
        set(handle, 'ButtonDownFcn',@SelectNodeString)
        
        % Perform set up
        SelectNodeString(handle,fig);
        
    case 'Internal Boundary'
        
end

    
    
    

%------------------------------------------------------------------------------
% Get GUI data
%------------------------------------------------------------------------------
%guiFig = findobj('Tag','ADmesh Figure'); guiH = % guidata(guiFig);


%------------------------------------------------------------------------------
% Ask user what type of boundary condition they wold like to assign
%------------------------------------------------------------------------------

% Construct a questdlg with three options
% choice = questdlg('Would you like a dessert?', ...
% 	'Dessert Menu', ...
% 	'Ice cream','Cake','No thank you','No thank you');
% % Handle response
% switch choice
%     case 'Ice cream'
%         disp([choice ' coming right up.'])
%         dessert = 1;
%     case 'Cake'
%         disp([choice ' coming right up.'])
%         dessert = 2;
%     case 'No thank you'
%         disp('I''ll bring you your check.')
%         dessert = 0;
% end

end