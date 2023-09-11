function plist = RemoveFromPlist(p,plist,id)

% Now remove from points list
%id = id(2:end-1);

for q = 1:length(id)
    
    id2 = (p(id(q),1) == plist(:,1) & p(id(q),2) == plist(:,2));
    
    plist(id2,:) = nan;
    
end

end