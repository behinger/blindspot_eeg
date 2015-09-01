function XDesc = bs_limo_designMat_betterDesc(Xin)
if ~iscell(Xin)
    error('input needs to be a cell')
end
XDesc = [];
for k = 1:length(Xin)
    [~,~,~,matchstring] = regexp(Xin{k},'([^,]*)');
    XDesc = [XDesc matchstring];
    

end