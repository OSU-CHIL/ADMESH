function [output,status] = CoordinateConversion(app,input,mode,varargin)

output = input;
status = 0;

if isempty(input)    
    return;
end

if length(varargin) == 2
    cpplon = varargin{1};
    cpplat = varargin{2};
elseif isempty(varargin)
    try
        cpplon = input.cpplon;
        cpplat = input.cpplat;
    catch
        cpplon = [];
        cpplat = [];
    end
end

app.ProgressBarButton.Text = 'Converting coordinates...'; drawnow;
switch lower(mode)
    case 'forward' % Convert to XY (meters)
        if ~isempty(cpplon)
            output = Geo2Cart(input,cpplon,cpplat);
            status = 1;
        else
            [output,output.cpplon,output.cpplat] = Geo2Cart(input);
            status = 1;
        end

    case 'reverse' % Convert to lat/lon
        output = Cart2Geo(input);
        status = 1;

    case 'auto'
        if isfield(input,'cpplon') && ~isempty(input.cpplon)
            output = Geo2Cart(input,input.cpplon,input.cpplat);
            status = 1;
        else
            if isfield(input,'Poly') % PTS structure input
                x = vertcat(input.Poly.x);
                y = vertcat(input.Poly.y);
                x = x(~isnan(x));
                y = y(~isnan(y));

            elseif isfield(input,'Points') || isprop(input,'Points')% MESH or xyzFun
                x = input.Points(:,1);
                y = input.Points(:,2);

            else
                error('Unknown format of input is passed in CoordinateConversion');
            end

            if all(x < 180 & x > -180 & y < 90 & y > -90)

                msg = ['The input may be in geographic coordinate system (in degrees). '...
                    'Do you want to convert it into planar coordinate system (in meters)?'];
                choice = uiconfirm(app.UIFigure,msg,'Coordinate conversion',...
                    'Options',{'Convert coordinates','Do not convert coordinates'},'DefaultOption',1,'Icon','Warning');

                if strcmpi(choice,'Convert coordinates')
                    % Convert to XY (meters)
                    [output,output.cpplon,output.cpplat] = Geo2Cart(input);
                    status = 1;
                end

            end
        end

    otherwise
        error('The mode must be one of "forward", "reverse", and "auto".')

end



