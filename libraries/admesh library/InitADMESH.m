function InitADMESH(app)

clc;

% Compile mex functions in admesh library
CompileMEX; 

% Hide abandoned tab. Keep it for backup
app.CrossSectionbetaTab.Parent = [];