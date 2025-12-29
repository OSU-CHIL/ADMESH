function id_MA_new = ThinningMA(id_MA,UIFigure)

msg = 'Thinning medial axis...';
uiprogressdlg(UIFigure,'Title','ADMESH','Message',msg,'Indeterminate','on');

bw = full(id_MA);

%--------------------------------------------------------------------------
% Fill holes with less than 10 pixels
%--------------------------------------------------------------------------
conn = 4;
filled = imfill(bw,conn,'holes');
holes = filled & ~bw;
bigholes = bwareaopen(holes, 10,conn);
smallholes = holes & ~bigholes;
bw = bw | smallholes;
clear filled holes bigholes smallholes;

%--------------------------------------------------------------------------
% Get skeleton of mask
% no idea why.. but bwskel looks better than bwmorph...
%--------------------------------------------------------------------------
% id_MA_thinned = bwmorph(bw,'thin');
id_MA_new = bwskel(bw);
% id_MA_thinned = bwmorph(bw,'skel',inf);

%--------------------------------------------------------------------------
% Remove isolated pixels
%--------------------------------------------------------------------------
id_MA_new = bwmorph(id_MA_new,'clean',inf);

id_MA_new = sparse(id_MA_new);
