function PTS = DrawOpenOceanBoundary(PTS,gui,pH,sb)

%--------------------------------------------------------------------------
% Plot Edge structure end points
%--------------------------------------------------------------------------
% Plot beginning and end points
axes(pH)
hp1 = line(PTS.Poly(1).x(1),PTS.Poly(1).y(1),'color','r','marker','+');
hp2 = line(PTS.Poly(1).x(end),PTS.Poly(1).y(end),'color','b','marker','+');

sb.setText(['Select boundary points going from the' ...
    ' start point to the end point. Right click to continue'])

leg = legend([hp2 hp1], 'Start Point', 'End Point');

drawnow

%------------------------------------------------------------------------------------
% Set window, button and ket functions
%------------------------------------------------------------------------------------
set(gui.Window, 'windowbuttonmotionfcn', @MouseMovement);   % Track mouse movement
set(gui.Window, 'windowkeypressfcn'    , @DeletePoint)
set(gui.Window, 'windowbuttondownfcn'  , @MouseClick);      % Track mouse clicks
set(gui.Window, 'Pointer'              , 'crosshair');      % Set pointer type for user

%------------------------------------------------------------------------------------
% Initially global variables, the list of points is empty.
%------------------------------------------------------------------------------------
X = []; Y = []; xInd = 0; yInd = 0; pHandle = nan; st = gui.st;

%------------------------------------------------------------------------------------
% Start (x,y) in PTS
%------------------------------------------------------------------------------------
sx = PTS.Poly(1).x(end); sy = PTS.Poly(1).y(end);
ex = PTS.Poly(1).x(1); ey = PTS.Poly(1).y(1);

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
        
        
        if strcmp(get(gui.Window,'SelectionType'),'alt') % Right click
            
            % Transpose X & Y to column vectors
            X = X';
            Y = Y';
            
            % Reset figure functions
            delete(leg)
            set(gui.Window, 'Pointer','arrow');
            set(gui.Window, 'windowbuttonmotionfcn', @CoordDisplay);   % Track mouse movement
            set(gui.Window, 'windowbuttondownfcn'  , '');      % Track mouse clicks
            set(gui.Window, 'windowkeypressfcn'          , '');
            
            title('');
            
            % Create IBtype
            PTS.Constraints(1,1) = struct('num',[],'type',[],'xy',[],'data',[]);
            PTS.Constraints(1).type = 'Open Ocean';
            PTS.Constraints(1).num  = -1;
            PTS.Constraints(1).xy = [[PTS.Poly(1).x(end); X; PTS.Poly(1).x(1) ],...
                                     [PTS.Poly(1).y(end); Y; PTS.Poly(1).y(1) ] ];
                        
            % Comlplete edge structure, PTS
            PTS.Poly(1).x = [PTS.Poly(1).x; X; PTS.Poly(1).x(1) ];
            PTS.Poly(1).y = [PTS.Poly(1).y; Y; PTS.Poly(1).y(1) ];

            uiresume;
            
                        
        else % Left click
            
            X(end+1) = xInd;
            Y(end+1) = yInd;
            
            if ishandle(pHandle)
                delete(pHandle)
            end
            
            axes(pH)
            pHandle = line([sx X ex],[sy Y ey],'color','b','Marker','.');
            
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
            st.setText('')
            return;
        end
        
        % update figure title
        st.setText(['( ' num2str(xInd,'%.3f') ' , ' num2str(yInd,'%.3f') ' )'])
         
        drawnow
        
    end

%------------------------------------------------------------------------------------
% Delete Point callback
%------------------------------------------------------------------------------------
    function DeletePoint(varargin)
       
        [~,eventData] = deal(varargin{:});
     
        % If d is pressed, delete last point in [X,Y]
        if strcmpi(eventData.Key, 'd') && ~isempty(X)
         
            delete(pHandle); % Delete plot
            
            % Delete last point in [X,Y]
            X(end) = [];
            Y(end) = [];
            
            axes(pH)
            pHandle = line([sx X ex],[sy Y ey],'color','b','Marker','.');

        end
        
        
    end

end