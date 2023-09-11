function PTS = CloseInternalBoundary(PTS,guiFig,pH,i)

%--------------------------------------------------------------------------
% Plot Edge structure end points
%--------------------------------------------------------------------------
% Plot beginning and end points
hp1 = plot(pH, PTS.Poly(i).x(1),PTS.Poly(i).y(1),'r+');
hp2 = plot(pH, PTS.Poly(i).x(end),PTS.Poly(i).y(end),'b+');

title(pH, ['Select boundary points going from the' ...
    ' start point to the end point. Right click to continue']);
leg = legend([hp2 hp1], 'Start Point', 'End Point');

drawnow

%------------------------------------------------------------------------------------
% Set window, button and ket functions
%------------------------------------------------------------------------------------
set(guiFig, 'windowbuttonmotionfcn', @MouseMovement);   % Track mouse movement
set(guiFig, 'windowkeypressfcn'    , @DeletePoint)
set(guiFig, 'windowbuttondownfcn'  , @MouseClick);      % Track mouse clicks
set(guiFig, 'Pointer'              , 'crosshair');      % Set pointer type for user

%------------------------------------------------------------------------------------
% Initially global variables, the list of points is empty.
%------------------------------------------------------------------------------------
X = []; Y = []; xInd = 0; yInd = 0; pHandle = nan;

%------------------------------------------------------------------------------------
% Start (x,y) in PTS
%------------------------------------------------------------------------------------
sx = PTS.Poly(i).x(end); sy = PTS.Poly(i).y(end);
ex = PTS.Poly(i).x(1); ey = PTS.Poly(i).y(1);

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
        
        
        if strcmp(get(guiFig,'SelectionType'),'alt') % Right click
            
            % Transpose X & Y to column vectors
            X = X';
            Y = Y';
            
            % Reset figure functions
            delete(leg)
            set(guiFig, 'Pointer','arrow');
            set(guiFig, 'windowbuttonmotionfcn', @PlotHoverCallback);   % Track mouse movement
            set(guiFig, 'windowbuttondownfcn'  , '');      % Track mouse clicks
            set(guiFig, 'windowkeypressfcn'          , '');
            
            title('');
                                    
            % Comlplete edge structure, PTS
            PTS.Poly(i).x = [PTS.Poly(i).x; X; PTS.Poly(i).x(1) ];
            PTS.Poly(i).y = [PTS.Poly(i).y; Y; PTS.Poly(i).y(1) ];

            uiresume;
            
                        
        else % Left click
            
            X(end+1) = xInd;
            Y(end+1) = yInd;
            
            if ishandle(pHandle)
                delete(pHandle)
            end
            
            pHandle = plot(pH,[sx X ex],[sy Y ey],'g.-');
            
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
            xlabel(pH,'');
            return;
        end
        
        % update figure title
        xlabel(...
            pH,['X = ' num2str(xInd,'%.3f') ', Y = ' num2str(yInd,'%.3f')],...
            'FontName', get(0,'defaultTextFontName'),'color','k');
        
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
            
            pHandle = plot(pH,[sx X ex],[sy Y ey],'g.-');

        end
        
        
    end

end