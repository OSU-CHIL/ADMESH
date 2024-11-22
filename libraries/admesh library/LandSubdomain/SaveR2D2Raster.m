function SaveR2D2Raster(app)

msg = 'Save R2D2 DEM....';
progdlg = uiprogressdlg(app.UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

f_dummy = figure('Position',[-100 -100 0 0]); % Create a dummy figure so that uigetfile doesn't minimize our GUI
[file,path] = uiputfile('*.tif','Save R2D2 DEM As');
delete(f_dummy); % Delete the dummy figure
figure(app.UIFigure); % Put focus on ADMESH app

outfilename = [path,file];

infilename = app.ElevationDataFilename;

[~,R] = readgeoraster(infilename);
info = geotiffinfo(infilename);

A = flipud(-app.xyzFun_new.Values');

geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
tiffTags = struct(TileLength=1024,TileWidth=1024);

geotiffwrite(outfilename,A,R,TiffType="bigtiff", ...
    GeoKeyDirectoryTag=geoTags,TiffTags=tiffTags)