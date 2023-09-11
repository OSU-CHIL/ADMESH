function bpoly = DrawBoundingPolygon(fig,pH,linecolor)

gui = % guidata(fig);

%--------------------------------------------------------------------------
% Update title
%--------------------------------------------------------------------------
% gui.sb.setText('Draw a bounding polygon around the section of the domain to be extracted')

%------------------------------------------------------------------------------------
% Set window, button and ket functions
%------------------------------------------------------------------------------------
set(fig, 'windowbuttonmotionfcn', @MouseMovement);   % Track mouse movement
set(fig, 'windowkeypressfcn'    , @DeletePoint)
set(fig, 'windowbuttondownfcn'  , @MouseClick);      % Track mouse clicks
set(fig, 'Pointer'              , 'crosshair');      % Set pointer type for user

%------------------------------------------------------------------------------------
% Initially global variables, the list of points is empty.
%------------------------------------------------------------------------------------
X = []; Y = []; xInd = 0; yInd = 0; pHandle = nan;

%------------------------------------------------------------------------------------
% Wait until the user chooses to continue or cancels
%------------------------------------------------------------------------------------
uiwait;

%------------------------------------------------------------------------------------
% Subroutines************************************************************************
%------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------
% Mouse Movement callback
%------------------------------------------------------------------------------------
    function MouseClick(varargin)
        
        if strcmp(get(fig,'SelectionType'),'alt') % Right click
            
            % Transpose X & Y to column vectors
            X = X';
            Y = Y';
            
            % Reset figure functions
            if ishandle(pHandle)
                delete(pHandle)
            end
            set(fig, 'Pointer','arrow');
            set(fig, 'windowbuttonmotionfcn', @CoordDisplay);   % Track mouse movement
            set(fig, 'windowbuttondownfcn'  , '');      % Track mouse clicks
            set(fig, 'windowkeypressfcn'          , '');
            
            bpoly = [X,Y];

            % gui.sb.setText('Ready')
            
            drawnow;
            uiresume;
            
                        
        else % Left click
            
            gui.sb.ProgressBar.setVisible(true)
            gui.sb.ProgressBar.setIndeterminate(true);
            % gui.sb.setText('Adding points...')
            
            % Add selection to vector
            X(end+1) = xInd;
            Y(end+1) = yInd;
            
            % Plot user selection
            if ishandle(pHandle)
                set(pHandle,'Xdata',X,'Ydata',Y)
                uistack(pHandle,'top')
                drawnow
            else
                pHandle = line(...
                    'Parent',pH,...
                    'Xdata',X,...
                    'Ydata',Y,...
                    'LineWidth',2,...
                    'Marker','o',...
                    'Color',linecolor,...
                    'MarkerFaceColor',linecolor);
                uistack(pHandle,'top')
                drawnow
            end
            
            gui.sb.ProgressBar.setVisible(false)
            gui.sb.ProgressBar.setIndeterminate(false);
            % gui.sb.setText('Draw a bounding polygon around the section of the domain to be extracted')
            
        end

    end

%------------------------------------------------------------------------------------
% Mouse Movement callback
%------------------------------------------------------------------------------------
    function MouseMovement(varargin)
        
        % Get current point
        pt = get(pH, 'CurrentPoint');
        xInd = pt(1, 1);
        yInd = pt(1, 2);
        
        % check if its within axes limits
        xLim = get(pH, 'XLim');
        yLim = get(pH, 'YLim');
        
        if xInd < xLim(1) || xInd > xLim(2) || yInd < yLim(1) || yInd > yLim(2)
            return;
        end
        
        gui.st.setText(['( ' num2str(xInd,'%.3f') ' , ' num2str(yInd,'%.3f') ' )'])
                
        drawnow
        
    end

%------------------------------------------------------------------------------------
% Delete Point callback
%------------------------------------------------------------------------------------
    function DeletePoint(varargin)
       
        [~,eventData] = deal(varargin{:});
     
        % If d is pressed, delete last point in [X,Y]
        if strcmpi(eventData.Key, 'd') && ~isempty(X)
            
            gui.sb.ProgressBar.setVisible(true)
            gui.sb.ProgressBar.setIndeterminate(true);
            % gui.sb.setText('Deleting previous selection...')
            
            % Delete last point in [X,Y]
            X(end) = [];
            Y(end) = [];
            
            % Re-plot line
            set(pHandle,'Xdata',X,'Ydata',Y)
            
            gui.sb.ProgressBar.setVisible(false)
            gui.sb.ProgressBar.setIndeterminate(false);
            % gui.sb.setText('Draw a bounding polygon around the section of the domain to be extracted')
            
            drawnow

            
        end
        
        
    end

end