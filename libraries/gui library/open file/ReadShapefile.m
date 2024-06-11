function [app,status] = ReadShapefile(file,app)

% Read SHP file
try
    SHP = shaperead(file);
    status = 1;
catch
    status = 0;
    return;
end

% Add SHP file to PTS variable
PTS = app.PTS;
if isfield(PTS,'Poly')
    n = length(PTS.Poly);
else
    n = 0;
end

% Add empty attributes if original PTS was empty
if isempty(PTS)
    PTS.Poly = [];
    PTS.Constraints = [];
    PTS.cpplon = [];
    PTS.cpplat = [];
end

flag = 1;
for i = 1 : length(SHP)
    % check if attribute already exists in PTS
    for j = 1 : length(PTS.Poly)
        if isequal(PTS.Poly(j).x, SHP(i).X(:)) && isequal(PTS.Poly(j).y, SHP(i).Y(:))
            flag = 0;
            break;
        end
    end

    % Add to PTS if it's new
    if flag
        PTS.Poly(n+i).x = SHP(i).X(:);
        PTS.Poly(n+i).y = SHP(i).Y(:);
    end
end

% Save to app
app.PTS = PTS;

