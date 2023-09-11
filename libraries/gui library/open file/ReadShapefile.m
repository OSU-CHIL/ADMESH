function [app,status] = ReadShapefile(file,app)

status = 1;
try
    SHP = shaperead(file);
catch
    status = 0;
end
    
if status == 0
    return;
end
for i = 1 : length(SHP)
    PTS.Poly(i).x = SHP(i).X(:);
    PTS.Poly(i).y = SHP(i).Y(:);    
end

PTS.Constraints = [];
PTS.cpplon = [];
PTS.cpplat = [];

app.PTS = PTS;