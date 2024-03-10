function id_MA_raw = ComputeMA_DDM(Vxg,Vyg,dd_ID)

id_MA_raw = sparse(size(Vxg,1),size(Vxg,2));
id_MA_raw = logical(id_MA_raw);

wbar = waitbar(0);
fwbar = @(x,y) waitbar(x/y,wbar,sprintf('Wairbar (%d/%d)',x,y));
for iDD = 1 : size(dd_ID,1)   
    Vxg1 = full(Vxg(dd_ID{iDD,1},dd_ID{iDD,2}));
    Vyg1 = full(Vyg(dd_ID{iDD,1},dd_ID{iDD,2}));

    DIV1 = (divergence(Vxg1,Vyg1));
    id_MA_raw1 = DIV1 > 0;
    id_MA_raw(dd_ID{iDD,1},dd_ID{iDD,2}) = id_MA_raw(dd_ID{iDD,1},dd_ID{iDD,2}) + id_MA_raw1;

    fwbar(iDD,size(dd_ID,1))
end
delete(wbar);