function SelectNodeString(varargin)

%---------------------------------------------------------------------
% Initialize persistent variables
%---------------------------------------------------------------------
persistent Handle gui fig pH nodeCount Snode Mnode Enode h xd yd L

%------------------------------------------------------------------------------
% Check input
%------------------------------------------------------------------------------
if nargin == 3
    % Action input
    [~,~,action] = deal(varargin{:});
else
    % Set up input
    [Handle,fig] = deal(varargin{:});
    gui = % guidata(fig);
    action = 'Set Up';
end


switch action
    
    case 'Set Up'
        
        % Assign ADMESH plot window handle
        pH        = findobj('Tag','ADMESH Plot Window');
        
        % Initialize node count variable
        nodeCount = 0;
        
        % Change line style on the polygon
        set(Handle, 'marker','.','markersize',15)
        
        % Assign x data
        xd = get(Handle,'XData')';
        
        % Assign y data
        yd = get(Handle,'YData')';
        
        % Assign length of data
        L = length(xd);
        
        set(Handle, 'ButtonDownFcn',{@SelectNodeString,'Select Node'})
        
    case 'Select Node'
        
        if nodeCount == 0 % Select first point
            
            % Get the current point on the screen
            currPt = get(pH,'CurrentPoint');
            
            x = currPt(1,1); y = currPt(1,2);
            
            % Compute distance & choose points within min(d)
            [~,Snode] = min((x - xd).^2 + (y - yd).^2);
            
            if isempty(Snode);return;end
            
            % Shift coordinate list so Snode is in the middle of the list
            i2 = floor(L/2); % half to the vector length
            
            % Shift to middle
            if Snode < i2 % shift up
                
                shift = i2 - Snode;
                
                xd = circshift(xd,shift);
                yd = circshift(yd,shift);
                
                Snode = i2;
                
            elseif Snode > i2% shift down
                
                shift = Snode - i2;
                
                xd = circshift(xd,-shift);
                yd = circshift(yd,-shift);
                
                Snode = i2;
                
            end
            
            % Plot the point
            h(1) = plot(pH,xd(Snode),yd(Snode),'b.','markersize',25);
            
            % Increment node count
            nodeCount = 1;
            
        elseif nodeCount == 1 % Select middle point
            
            % Get the current point on the screen
            currPt = get(pH,'CurrentPoint');
            
            x = currPt(1,1); y = currPt(1,2);
            
            % Compute distance & choose points within min(d)
            [~,Mnode] = min((x - xd).^2 + (y - yd).^2);
            
            if isempty(Mnode);return;end
            
            % Plot the point
            h(2) = plot(pH,xd(Mnode),yd(Mnode),'b.','markersize',25);
            
            % Increment node count
            nodeCount = 2;
            
        elseif nodeCount == 2
            
            % Connect the line
            currPt = get(pH,'CurrentPoint');
            
            x = currPt(1,1); y = currPt(1,2);
            
            % Compute distance & choose points within min(d)
            [~,Enode] = min((x - xd).^2 + (y - yd).^2);
            
            % What direction are we going? Build node string
            if Mnode > Snode
                
                if Enode > Mnode
                    
                    % Build node string
                    nStr = (Snode:Enode);
                    
                else
                    
                    % Build node string
                    nStr = [Snode:L, 1:Enode];
                    
                end
                
                
            elseif Mnode < Snode
                
                if Enode < Mnode
                    
                    % Build node string
                    nStr = (Snode:-1:Enode);
                    
                else
                    
                    % Build node string
                    nStr = [Snode:-1:1, L:-1:Enode];
                    
                end
                
                
            end
            
            % Remove double points if they exist
            nStr = unique(nStr,'stable');
            
            
            % Delete points and plot the new Line
            delete(h)
            
            h = plot(pH,xd(nStr),yd(nStr),'b-','markersize',25);
                   
            gui.PTS.Constraints(end+1).xy = [xd(nStr),yd(nStr)];
            gui.PTS.Constraints(end).num = -1;
            gui.PTS.Constraints(end).type = 'Open Ocean';
            gui.PTS.Constraints(end).data = [];

            delete(h)
            
            % Clear mouse click function
            set(Handle, 'ButtonDownFcn','')
            
            % Remove markers
            set(Handle, 'marker','none')
            
            % guidata(fig,gui)
            gui = % guidata(fig);
            
            PlotConstraints(gui.PTS,gui.ViewAxes);
            
        end
        
        
end


end